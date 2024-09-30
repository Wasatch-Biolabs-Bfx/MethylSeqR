## Sliding Window

# Dependencies
library(tidyr)
library(dplyr)
library(dbplyr)
library(duckdb)
library(duckplyr)
library(progress)

# Function
summarize_windows <- function(ch3_db,
                         call_type = "positions",
                         window_size = 1000,
                         step_size = 10,
                         overwrite = TRUE) 
{
  # Open the database connection
  db_con <- helper_connectDB(ch3_db)
  
  if (dbExistsTable(db_con, "windows") & overwrite)
    dbRemoveTable(db_con, "windows")
  
  if (dbExistsTable(db_con, "temp_table"))
    dbRemoveTable(db_con, "temp_table")
  
  # Calc windows in each frame
  offsets <- seq(1, window_size - 1, by = step_size)
  
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
  pb$tick()
  
  # Conduct analysis. 
  # Creates tiled windows and then loops to create sliding window
  for (offset in offsets) {
    tbl(
      db_con, call_type) |>
      mutate(
        start = ref_position - ((ref_position - offset) %% window_size)) |>
      filter(
        start > 0) |>
      summarize(
        .by = c(sample_name, chrom, start),
        total_calls = sum(cov),
        across(ends_with("_counts"), ~sum(.x, na.rm = TRUE)),
        across(ends_with("_frac"), ~sum(.x * cov, na.rm = TRUE) / sum(cov))) |>
      mutate(
        end = start + window_size - 1) |>
      compute(name = "temp_table", temporary = TRUE)
    
    # Create or append table
    dbExecute(db_con, 
              "CREATE TABLE IF NOT EXISTS windows AS 
              SELECT * FROM temp_table WHERE 1=0")
    
    dbExecute(db_con, 
              "INSERT INTO windows SELECT * FROM temp_table")
    
    dbRemoveTable(db_con, "temp_table")
    
    # Advance progress bar
    pb$tick()
  }
  
  # Close progress bar
  pb$terminate()
  
  # Finish Up
  ch3_db$tables <- dbListTables(db_con) # Update table list
  dbDisconnect(db_con, shutdown = TRUE)
  return(ch3_db)
}