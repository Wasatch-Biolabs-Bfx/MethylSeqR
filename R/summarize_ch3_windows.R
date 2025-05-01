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
#' @param mod_type A character vector specifying the modification types to include. Options are  `"c"` (unmodified cytosine),
#' `"m"` (methylation), `"h"` (hydroxymethylation), 
#'   and `"mh"` (methylated + hydroxymethylated).
#' @param chrs A character vector specifying which chromosomes to include. Default includes all autosomes, 
#'   sex chromosomes (chrX, chrY), and mitochondrial chromosome (chrM).
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
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbExistsTable dbRemoveTable
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl
#' @importFrom glue glue glue_collapse
#' 
#' @examples
#'  # Specify the path to the database
#' ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#' 
#' # Summarize Windows
#' summarize_ch3_windows(ch3_db, window_size = 100, step_size = 100)
#'
#' @export
summarize_ch3_windows <- function(ch3_db,
                         window_size = 1000,
                         step_size = 10,
                         mod_type = c("c", "m", "h", "mh"),
                         chrs = c(as.character(1:22), 
                                  paste0("chr", 1:22), "chrX", "chrY", "chrM",
                                  paste0("Chr", 1:22), "ChrX", "ChrY", "ChrM"),
                         min_num_calls = 1,
                         overwrite = TRUE) 
{
  # Open the database connection
  database <- .ch3helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # Finish up: purge extra tables & update table list and close the connection
  on.exit(.ch3helper_purgeTables(db_con), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(.ch3helper_closeDB(database), add = TRUE)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
  
  if (dbExistsTable(db_con, "windows") & overwrite)
    dbRemoveTable(db_con, "windows")
  
  if (dbExistsTable(db_con, "temp_table"))
    dbRemoveTable(db_con, "temp_table")
  
  query <- glue("
    CREATE TABLE temp_positions AS 
    SELECT
        sample_name,
        chrom,
        start,
        \"end\",
        COUNT(*) AS num_calls,
        SUM(CASE WHEN call_code = '-' THEN 1 ELSE 0 END) AS c_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) AS m_counts,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) AS h_counts,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) AS mh_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS m_frac,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS h_frac,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS mh_frac
    FROM calls
    WHERE chrom IN ({glue::glue_collapse(glue(\"'{chrs}'\"), sep = ', ')})
    GROUP BY sample_name, chrom, start, \"end\"
    HAVING num_calls >= {min_num_calls};  -- Filter based on min_num_calls
")
  
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_positions;")  # Drop existing table
  dbExecute(db_con, "VACUUM;")  # Clean up storage
  dbExecute(db_con, query)  # Execute the query
  
  # Calc windows in each frame
  offsets <- seq(1, window_size - 1, by = step_size)
  
  cat("Building windows table...")
  
  # Conduct analysis. 
  # Creates tiled windows and then loops to create sliding window
  for (offset in offsets) {
    .make_window(db_tbl, db_con, offset, window_size)
  }

  if (dbExistsTable(db_con, "temp_table"))
    dbRemoveTable(db_con, "temp_table")
  # Close progress bar
  pb$terminate()
  
  message("Windows table successfully created!")
  print(head(tbl(db_con, "windows")))
  
  invisible(database)
}

.make_window <- function(db_tbl, db_con, offset, window_size)
{
  query <- glue::glue("
    CREATE TEMP TABLE temp_table AS
    WITH windowed AS (
      SELECT 
        sample_name,
        chrom,
        start - ((start - {offset}) % {window_size}) AS start,
        SUM(num_calls) AS num_calls,
        {paste0('SUM(', c('c_counts', 'm_counts', 'h_counts', 'mh_counts'), ') AS ', c('c_counts', 'm_counts', 'h_counts', 'mh_counts'), collapse = ', ')},
        {paste0('SUM(', c('m_counts', 'h_counts', 'mh_counts'), ' * num_calls) / NULLIF(SUM(num_calls), 0) AS ', c('m_frac', 'h_frac', 'mh_frac'), collapse = ', ')}
      FROM temp_positions
      WHERE start > 0
      GROUP BY sample_name, chrom, start
    )
    SELECT 
      sample_name,
      chrom,
      start,
      start + {window_size} - 1 AS end,
      COUNT(*) AS num_CpGs,
      SUM(num_calls) AS num_calls,
      {paste0('SUM(', c('c_counts', 'm_counts', 'h_counts', 'mh_counts'), ') AS ', c('c_counts', 'm_counts', 'h_counts', 'mh_counts'), collapse = ', ')},
      {paste0('SUM(', c('m_counts', 'h_counts', 'mh_counts'), ' * num_calls) / NULLIF(SUM(num_calls), 0) AS ', c('m_frac', 'h_frac', 'mh_frac'), collapse = ', ')}

    FROM windowed
    GROUP BY sample_name, chrom, start
  ")

  dbExecute(db_con, query)
  
  dbExecute(db_con, "CREATE TABLE IF NOT EXISTS windows AS SELECT * FROM temp_table WHERE 1=0")
  dbExecute(db_con, "INSERT INTO windows SELECT * FROM temp_table")
  dbRemoveTable(db_con, "temp_table")
}