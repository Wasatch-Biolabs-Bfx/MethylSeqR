#' Summarize Methylation Data using a Sliding Window
#'
#' This function summarizes methylation data from a DuckDB database by creating 
#' sliding windows over the specified genomic regions. It allows for the adjustment 
#' of window size and step size to control the granularity of the summarization.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param call_type A string indicating the type of data to summarize. Default is "positions".
#' @param window_size An integer specifying the size of the sliding window in base pairs. Default is 1000.
#' @param step_size An integer specifying the number of base pairs to step forward with each window. Default is 10.
#' @param overwrite A logical indicating whether to overwrite the existing "windows" table if it exists. Default is TRUE.
#'
#' @details
#' The function connects to a DuckDB database and removes any existing "windows" and "temp_table" tables if necessary. 
#' It creates a sequence of offsets based on the specified window and step sizes, and then it iterates through 
#' these offsets to generate sliding windows of methylation data. A progress bar is displayed during the operation.
#'
#' The function utilizes the helper function `.make_window` to perform the actual window calculation and summarization.
#' The resulting summarized data is stored in a table called "windows" within the database.
#'
#' @return The updated `ch3_db` object with the summarized windows data added to the DuckDB database.
#'
#' @import dplyr
#' @import dbplyr
#' @import duckdb
#' @import duckplyr
#' @import progress
#' 
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Calculate windows
#'  summarize_windows(ch3_db, window_size = 100, step_size = 100)
#'
#' @export
summarize_windows <- function(ch3_db,
                         call_type = "positions",
                         window_size = 1000,
                         step_size = 10,
                         overwrite = TRUE) 
{
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  if (dbExistsTable(db_con, "windows") & overwrite)
    dbRemoveTable(db_con, "windows")
  
  if (dbExistsTable(db_con, "temp_table"))
    dbRemoveTable(db_con, "temp_table")
  
  # Open table
  db_tbl <- tbl(db_con, call_type)
  
  # Calc windows in each frame
  offsets <- seq(1, window_size - 1, by = step_size)
  
  cat("Building windows table...")
  # Create Progress Bar
  pb <- progress_bar$new(
    format = "[:bar] :percent [Elapsed time: :elapsed]",
    total = length(offsets) + 1,
    complete = "=",   
    incomplete = "-", 
    current = ">",    
    clear = FALSE,    
    width = 100)
  
  # Tick progress bar to make it show up (first loop can be long)
  pb$tick(-1)
  pb$tick()
  
  
  # Conduct analysis. 
  # Creates tiled windows and then loops to create sliding window
  tryCatch(
    {
      for (offset in offsets) {
        .make_window(db_tbl, db_con, offset, window_size)
        # Advance progress bar
        pb$tick()
      }
    }, 
    error = function(e) 
    {
      if (dbExistsTable(db_con, "windows"))
        dbRemoveTable(db_con, "windows")
      message("Error: ", e$message)
    },
    finally = 
      {
        if (dbExistsTable(db_con, "temp_table"))
          dbRemoveTable(db_con, "temp_table")
        # Close progress bar
        pb$terminate()
        # Finish Up
        database$last_table = "windows"
        .helper_closeDB(database)
        return(database)
      })
}

.make_window <- function(db_tbl, db_con, offset, window_size)
{
  db_tbl |>
    mutate(
      start = ref_position - ((ref_position - offset) %% window_size),
      na.rm = TRUE) |>
    filter(
      start > 0) |>
    summarize(
      .by = c(sample_name, chrom, start),
      total_calls = sum(cov, na.rm = TRUE),
      across(ends_with("_counts"), ~sum(.x, na.rm = TRUE)),
      across(ends_with("_frac"), ~sum(.x * cov, na.rm = TRUE) / sum(cov, na.rm = TRUE))) |>
    mutate(
      end = start + window_size - 1) |>
    select(sample_name, chrom, start, end, everything()) |>
    compute(name = "temp_table", temporary = TRUE)
  
  # Create or append table
  dbExecute(db_con, 
            "CREATE TABLE IF NOT EXISTS windows AS 
              SELECT * FROM temp_table WHERE 1=0")
  
  dbExecute(db_con, 
            "INSERT INTO windows SELECT * FROM temp_table")
  
  dbRemoveTable(db_con, "temp_table")
}