#' Summarize Methylation Positions from Calls Table
#'
#' This function processes methylation call data from a DuckDB database, summarizing 
#' coverage and call counts into a `positions` table.
#'
#' @param ch3_db A DuckDB database connection or file path (character) to the `.ch3.db` file.
#' @param mod_type A character vector specifying the modification types to include. Options are  `"c"` (unmodified cytosine),
#' `"m"` (methylation), `"h"` (hydroxymethylation), 
#'   and `"mh"` (methylated + hydroxymethylated).
#' @param chrs A character vector specifying which chromosomes to include. Default includes all autosomes, 
#'   sex chromosomes (chrX, chrY), and mitochondrial chromosome (chrM).
#' @param min_num_calls Minimum number of calls required to include a position in the summary. Default is 1.
#'
#' @return The modified `ch3_db` object with the updated `positions` table.
#' @details
#' The function:
#' - Connects to the DuckDB database.
#' - Summarizes methylation calls by `sample_name`, `chrom`, `start`, and `end`.
#' - Computes coverage and call counts (`m`, `h`, `mh`, `c`).
#' - Filters based on `min_num_calls`.
#' - Creates a `positions` table in the database.
#'
#' A progress bar is displayed during execution.
#'
#' @examples
#' # Specify the path to the database
#' ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#' 
#' # Summarize Positions
#' summarize_positions(ch3_db)
#'
#' @importFrom DBI dbExecute dbDisconnect dbListTables
#' @importFrom dplyr tbl summarize mutate filter
#' @importFrom progress progress_bar
#' 
#' @export
summarize_positions <- function(ch3_db,
                                mod_type = c("c", "m", "h", "mh"),
                                chrs = c(as.character(1:22), 
                                         paste0("chr", 1:22), "chrX", "chrY", "chrM",
                                         paste0("Chr", 1:22), "ChrX", "ChrY", "ChrM"),
                                min_num_calls = 1) 
{
  # Open the database connection - first check to make sure correct name is there
  if (is.character(ch3_db)) {
    if (!grepl(".ch3.db$", ch3_db)) {
      ch3_db <- paste0(ch3_db, ".ch3.db")
    }
  }
  
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # Finish up: purge extra tables & update table list and close the connection
  keep_tables = c("calls","positions", "regions", "windows", "mod_diff")
  on.exit(.helper_purgeTables(db_con, keep_tables), add = TRUE)
  on.exit(.helper_closeDB(database), add = TRUE)
  
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
    filter(chrom %in% chrs) |>  # Filter for selected chromosomes
    summarize(
      .by = c(sample_name, chrom, start, end),
      num_calls = n(),
      c_counts = sum(as.integer(call_code == "-"), na.rm = TRUE),
      m_counts = sum(as.integer(call_code == "m"), na.rm = TRUE),
      h_counts = sum(as.integer(call_code == "h"), na.rm = TRUE),
      mh_counts = sum(as.integer(call_code %in% c("m", "h")), na.rm = TRUE)) |>
    mutate(
      m_frac = m_counts / num_calls,
      h_frac = h_counts / num_calls,
      mh_frac = mh_counts / num_calls) |> 
    filter(num_calls >= min_num_calls)
  
  pb$tick()  # Progress bar update after summarizing
  
  # Select only requested modtype columns (always keeping num_calls)
  selected_columns <- c("sample_name", "chrom", "start", "end", "num_calls", 
                        paste0(mod_type, "_counts"), paste0(mod_type, "_frac"))
  selected_columns <- intersect(selected_columns, colnames(summarized_data))
  
  summarized_data <- summarized_data |> select(all_of(selected_columns))
  
  # Materialize summarized data into a temporary table before creating positions
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_summary;")
  copy_to(db_con, summarized_data, "temp_summary", temporary = TRUE)
  
  pb$tick()  # Progress bar update after summarizing
  
  # Drop the positions table if it already exists
  dbExecute(db_con, "DROP TABLE IF EXISTS positions;")
  
  # Create the final 'positions' table from the materialized data
  dbExecute(db_con, "CREATE TABLE positions AS SELECT * FROM temp_summary;")
  on.exit(dbExecute(db_con, "DROP TABLE IF EXISTS temp_summary;"), add = TRUE)
  
  pb$tick()  # Progress bar update after summarizing
  
  message("Positions table successfully created!")
  # message("\n")
  # message("Printing preview of positions table.")
  print(head(summarized_data))
  
  # ch3_db$tables <- dbListTables(db_con)
  invisible(ch3_db)
}