## Make Position DB

# Dependencies
library(arrow)
library(tidyr)
library(dplyr)
library(dbplyr)
library(duckdb)
library(duckplyr)
library(progress)

# Function
make_pos_db <- function(ch3_files, 
                        ch3_db,
                        min_call_prob = 0.9,
                        min_length = 100,
                        min_base_qual = 10) 
{
  # Setup files and db
  if (!grepl(".ch3.db$", ch3_db))
    ch3_db <- paste0(ch3_db, ".ch3.db")
  
  ch3_db <- list(db_file = ch3_db)
  class(ch3_db) <- "ch3_db"
  db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
  ch3_db$tables <- dbListTables(db_con)
  
  ch3_files <- list.files(ch3_files, pattern = "\\.ch3$", full.names = TRUE)
  
  # Check if the table already exists and delete it if it does
  if (dbExistsTable(db_con, "positions"))
    dbRemoveTable(db_con, "positions")
  
  if (dbExistsTable(db_con, "windows"))
    dbRemoveTable(db_con, "windows")
  
  # Loop through files to add to db
  # Create Progress Bar
  pb <- progress_bar$new(
    format = "[:bar] :percent [Elapsed time: :elapsed]",
    total = length(ch3_files) + 1,
    complete = "=",   
    incomplete = "-", 
    current = ">",    
    clear = FALSE,    
    width = 100)   
  
  pb$tick()
  
  for (ch3_file in ch3_files) {
    open_dataset(ch3_file) |>
      select(sample_name, chrom, ref_position, call_prob, 
             read_length, base_qual, call_code) |>
      filter(
        call_prob >= min_call_prob,
        read_length >= min_length,
        base_qual >= min_base_qual) |>
      summarize(
        .by = c(sample_name, chrom, ref_position),
        cov = n(),
        c_counts = sum(as.integer(call_code == "-"),
                       na.rm = TRUE),
        m_counts = sum(as.integer(call_code == "m"),
                       na.rm = TRUE),
        h_counts = sum(as.integer(call_code == "h"),
                       na.rm = TRUE),
        mh_counts = sum(as.integer(call_code %in% c("m", "h")),
                        na.rm = TRUE)) |>
      mutate(
        m_frac = m_counts / cov,
        h_frac = h_counts / cov,
        mh_frac = mh_counts / cov) |>
      arrange(
        sample_name, chrom, ref_position) |>
      collect() |>
      dbWriteTable(
        conn = db_con, 
        name = "positions", 
        append = TRUE)
    
    pb$tick()
  }
  
  # Close progress bar
  pb$terminate()
  
  # Finish Up
  ch3_db$tables <- dbListTables(db_con) # Update table list
  dbDisconnect(db_con, shutdown = TRUE)
  return(ch3_db)
}