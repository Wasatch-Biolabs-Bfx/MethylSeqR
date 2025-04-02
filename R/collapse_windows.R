#' Collapse Windows Based on Methylation Differences
#'
#' This function collapses significant windows in a methylation dataset by merging
#' contiguous regions that meet the specified criteria. Can only collapse windows 
#' once a differential modification analysis (mod_diff()) has been called.
#'
#' @param ch3_db A DuckDB database connection object or path to the database.
#' @param max_distance Numeric. The maximum allowable distance between consecutive
#'        significant windows for merging (default: 1000).
#' @param sig_cutoff Numeric. The significance threshold for adjusted p-values
#'        (default: 0.05).
#' @param min_diff Numeric. The minimum absolute methylation difference required
#'        for inclusion in the analysis (default: 0.5).
#' @param output_table_name Character. Name of the output table to store collapsed
#'        windows (default: "collapsed_windows").
#'
#' @return This function does not return an object; it creates or replaces the
#'         `collapsed_windows` table in the database.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters the `mod_diff_windows` to retain only significant windows where
#'         `p_adjust <= sig_cutoff` and `ABS(meth_diff) >= min_diff`.
#'   \item Assigns a new region identifier based on proximity (`max_distance`) and
#'         the direction of methylation differences.
#'   \item Collapses regions by grouping contiguous windows, computing the
#'         average methylation difference (`avg_meth_diff`), and counting the
#'         number of merged windows.
#' }
#'
#' @import DBI
#' @importFrom glue glue
#' @export
collapse_windows <- function(ch3_db, 
                             max_distance = 1000,
                             sig_cutoff = 0.05,
                             min_diff = 0.5,
                             output_table_name = "collapsed_windows") 
{
  database <- MethylSeqR:::.helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Check if "mod_diff" table exists
  if (!DBI::dbExistsTable(db_con, "mod_diff_windows")) {
    stop(glue::glue("Error: Table 'mod_diff_windows' not found in the database. 
                     Please run 'mod_diff()' on windows data first to generate it."))
  }
  
  on.exit(MethylSeqR:::.helper_purgeTables(db_con), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(MethylSeqR:::.helper_closeDB(database), add = TRUE)
  
  query <- 
    glue(
      "CREATE OR REPLACE TABLE {output_table_name} AS
    WITH FilteredWindows AS (
      SELECT *
      FROM mod_diff_windows
      WHERE p_adjust <= {sig_cutoff} AND ABS(meth_diff) >= {min_diff}
    ),
    NumberedWindows AS (
      SELECT *,
        CASE
          WHEN LAG(\"end\") OVER w IS NULL 
               OR LAG(\"end\") OVER w + {max_distance} < start 
               OR SIGN(meth_diff) != SIGN(LAG(meth_diff) OVER w) THEN 1
          ELSE 0
        END AS new_region_flag
      FROM FilteredWindows
      WINDOW w AS (PARTITION BY chrom ORDER BY start)
    ),
    RegionGroups AS (
      SELECT *,
        SUM(new_region_flag) OVER (PARTITION BY chrom ORDER BY start 
                                   ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) AS region_id
      FROM NumberedWindows
    )
    SELECT 
      chrom,
      MIN(start) AS start,
      MAX(\"end\") AS \"end\",
      AVG(meth_diff) AS avg_meth_diff,
      COUNT(*) AS num_windows
    FROM RegionGroups
    GROUP BY chrom, region_id
    ORDER BY chrom, start;")
  
  DBI::dbExecute(db_con, query)
  
  message(paste0("Windows successfully collapsed - ", output_table_name, " created!"))
  
  return(invisible(NULL))
}