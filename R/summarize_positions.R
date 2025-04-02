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
  
  # Connect to the database
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # purge extra tables, update table list, and then close the connection
  keep_tables = c("calls", "positions", "regions", "windows", 
                  "mod_diff_positions", "mod_diff_regions", "mod_diff_windows",
                  "collapsed_windows")
  on.exit(.helper_purgeTables(db_con, keep_tables), add = TRUE)  # Purge tables FIRST
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(.helper_closeDB(database), add = TRUE)        # Close DB LAST 
  
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
  
  # Process data using duckplyr
  cat("Building positions table...\n")
  
  query <- glue("
    CREATE TABLE positions AS 
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
    WHERE chrom IN ({glue_collapse(glue(\"'{chrs}'\"), sep = ', ')})
    GROUP BY sample_name, chrom, start, \"end\"
    HAVING num_calls >= {min_num_calls};  -- Filter based on min_num_calls
")
  
  dbExecute(db_con, "DROP TABLE IF EXISTS positions;")  # Drop existing table
  dbExecute(db_con, "VACUUM;")  # Clean up storage
  dbExecute(db_con, query)  # Execute the query
  
  message("Positions table successfully created!")
  # Print a preview of what table looks like
  print(head(tbl(db_con, "positions")))
  
  # ch3_db$tables <- dbListTables(db_con)
  invisible(ch3_db)
}