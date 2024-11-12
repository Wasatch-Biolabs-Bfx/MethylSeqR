#' Create Position Database from CH3 Files
#'
#' This function processes a collection of CH3 files to create a positions database, 
#' filtering and summarizing the methylation data based on specified criteria. 
#' The resulting database is stored in a DuckDB format, allowing for efficient queries.
#'
#' @param ch3_files A string representing the directory path containing CH3 files. 
#' The function will read all files with the ".ch3" extension in this directory.
#' @param ch3_db A string representing the path where the positions database will be created. 
#' The extension ".ch3.db" will be appended if not already present.
#' @param min_call_prob A numeric value representing the minimum call probability threshold. 
#' Only calls with a probability greater than or equal to this value will be included. Default is 0.9.
#' @param min_length A numeric value representing the minimum read length. 
#' Only reads meeting or exceeding this length will be processed. Default is 100.
#' @param min_base_qual A numeric value representing the minimum base quality. 
#' Only reads with quality scores at or above this threshold will be included. Default is 10.
#'
#' @details
#' The function connects to a DuckDB database and processes each CH3 file to extract relevant methylation data. 
#' The data is filtered according to the specified thresholds and summarized by sample name, chromosome, and reference position. 
#' The final summary statistics, including counts and fractions for methylated and hydroxymethylated bases, are written to the database.
#'
#' @return A list containing the database file path and the updated list of tables in the database.
#' This can be useful for further processing or analysis.
#'
#' @import arrow
#' @import dplyr
#' @import dbplyr
#' @import duckdb
#' @import duckplyr
#' @import progress
#'
#' @examples
#' # Set up the file path for the test data located in inst/test_data/
#' ch3_files <- system.file("test_data", package = "MethylseqR")
#' ch3_db <- tempfile("example_db")
#' 
#' # Run the function with the example data
#' result <- make_pos_db(ch3_files, ch3_db)
#' print(result)
#' 
#' @export
make_pos_db <- function(ch3_files, 
                        ch3_db,
                        min_call_prob = 0.9,
                        min_length = 100,
                        min_base_qual = 10) 
{
  # Setup files and db
  if (!grepl(".ch3.db$", ch3_db))
    ch3_db <- paste0(ch3_db, ".ch3.db")
  
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  ch3_files <- list.files(ch3_files, pattern = "\\.ch3$", full.names = TRUE)
  
  # Check if the table already exists and delete it if it does
  if (dbExistsTable(db_con, "positions"))
    dbRemoveTable(db_con, "positions")
  
  if (dbExistsTable(db_con, "windows"))
    dbRemoveTable(db_con, "windows")
  
  tryCatch(
    {
      # Loop through files to add to db
      # Create Progress Bar
      cat("Building positions table...")
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
            base_qual >= min_base_qual,
            nchar(chrom) < 6) |> # remove unneccessary chromosomes
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
    }, 
    error = function(e) 
    {
      # Print custom error message
      message("An error occurred: ", e$message)
    }, 
    finally = 
      {
        # Finish Up
        database$last_table = "positions"
        .helper_closeDB(database)
        return(database)
      })
  # return(ch3_db)
}