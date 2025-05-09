#' Summarize Methylation Positions from Calls Table
#'
#' This function processes methylation call data from a DuckDB database, summarizing 
#' coverage and call counts into a `positions` table.
#'
#' @param ch3_db A DuckDB database connection or file path (character) to the `.ch3.db` file.
#' @param table_name A string specifying what the user would like the name to be called in the database. Default is "positions".
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
#' summarize_ch3_positions(ch3_db)
#'
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbExistsTable dbRemoveTable
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl
#' @importFrom glue glue glue_collapse
#' 
#' @export

summarize_ch3_positions <- function(ch3_db,
                                    table_name = "positions",
                                mod_type = c("c", "m", "h", "mh"),
                                chrs = c(as.character(1:22), 
                                         paste0("chr", 1:22), "chrX", "chrY", "chrM",
                                         paste0("Chr", 1:22), "ChrX", "ChrY", "ChrM"),
                                min_num_calls = 1) 
{
  start_time <- Sys.time()
  # Connect to the database
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(ch3_db$con, "PRAGMA max_temp_directory_size='100GiB';")
  
  # Process data using duckplyr
  cat("Building positions table...\n")
  
  query <- glue("
    CREATE TABLE {table_name} AS 
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
  
  if (dbExistsTable(ch3_db$con, table_name))
    dbRemoveTable(ch3_db$con, table_name)
  
  dbExecute(ch3_db$con, "VACUUM;")  # Clean up storage
  dbExecute(ch3_db$con, query)  # Execute the query
  
  end_time <- Sys.time()
  message("Positions table successfully created as ", table_name, " in database!\n", 
          "Time elapsed: ", end_time - start_time, "\n")
  
  # Print a preview of what table looks like
  print(head(tbl(ch3_db$con, table_name)))
  
  ch3_db$current_table = table_name
  ch3_db <- .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}