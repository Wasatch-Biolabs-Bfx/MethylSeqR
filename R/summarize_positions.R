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
#' @param chrs A list representation of all chromosomes wanted to be included in the positional data.
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
#' @import progress
#'
#' @examples
#' # Set up the file path for the test data located in inst/test_data/
#' ch3_files <- system.file("test_data", package = "MethylSeqR")
#' ch3_db <- tempfile("example_db")
#' 
#' # Run the function with the example data
#' result <- make_pos_db(ch3_files, ch3_db)
#' print(result)
#' 
#' @export
summarize_positions <- function(ch3_db, 
                             min_cov = 1) 
{
  # Open the database connection - first check to make sure correct name is there
  if (is.character(ch3_db)) {
    if (!grepl(".ch3.db$", ch3_db)) {
      ch3_db <- paste0(ch3_db, ".ch3.db")
    }
  }
  
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  on.exit(ch3_db$tables <- dbListTables(ch3_db$db_con), add = TRUE)
  on.exit(dbDisconnect(db_con, shutdown = TRUE), add = TRUE)
  on.exit(ch3_db$db_con <<- NULL, add = TRUE)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
  
  # Set up the progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent [Elapsed: :elapsed]",
    total = 3,  # 3 major steps: summarizing, creating temp table, and creating positions table
    complete = "=",   
    incomplete = "-", 
    current = ">",    
    clear = FALSE,    
    width = 60
  )
  
  pb$tick(0)  # Initialize bar without moving it
  
  # Process data using duckplyr
  message("Building positions data from calls table...")
  
  # Process data using duckplyr
  summarized_data <- tbl(db_con, "calls") |>
    summarize(
      .by = c(sample_name, chrom, ref_position),
      cov = n(),
      c_counts = sum(as.integer(call_code == "-"), na.rm = TRUE),
      m_counts = sum(as.integer(call_code == "m"), na.rm = TRUE),
      h_counts = sum(as.integer(call_code == "h"), na.rm = TRUE),
      mh_counts = sum(as.integer(call_code %in% c("m", "h")), na.rm = TRUE)) |>
    mutate(
      m_frac = m_counts / cov,
      h_frac = h_counts / cov,
      mh_frac = mh_counts / cov) |> 
    filter(cov >= min_cov)
  
  pb$tick()  # Progress bar update after summarizing
  
  # Materialize summarized data into a temporary table before creating positions
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_summary;")
  copy_to(db_con, summarized_data, "temp_summary", temporary = TRUE)
  
  pb$tick()  # Progress bar update after summarizing
  
  # Drop the positions table if it already exists
  dbExecute(db_con, "DROP TABLE IF EXISTS positions;")
  
  # Create the final 'positions' table from the materialized data
  dbExecute(db_con, "CREATE TABLE positions AS SELECT * FROM temp_summary;")
  on.exit(dbExecute(db_con, "DROP TABLE IF EXISTS temp_summary;"))
  
  pb$tick()  # Progress bar update after summarizing
  
  message("Positions table successfully created!")
  # message("\n")
  # message("Printing preview of positions table.")
  print(head(summarized_data))
  
  # ch3_db$tables <- dbListTables(db_con)
  invisible(ch3_db)
}