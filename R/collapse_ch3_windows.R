#' Collapse Windows Based on Methylation Differences
#'
#' This function collapses significant windows in a methylation dataset by merging
#' contiguous regions that meet the specified criteria. Can only collapse windows 
#' once a differential modification analysis (calc_ch3_diff()) has been called.
#'
#' @param ch3_db A DuckDB database connection object or path to the database.
#' @param table_name Character. Name of the output table to store collapsed
#'        windows (default: "collapsed_windows").
#' @param max_distance Numeric. The maximum allowable distance between consecutive
#'        significant windows for merging (default: 1000).
#' @param sig_cutoff Numeric. The significance threshold for adjusted p-values
#'        (default: 0.05).
#' @param min_diff Numeric. The minimum absolute methylation difference required
#'        for inclusion in the analysis (default: 0.5).
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
#' @importFrom DBI dbConnect dbDisconnect dbExistsTable dbExecute
#' @importFrom duckdb duckdb
#' @importFrom glue glue
#' 
#' @export

collapse_ch3_windows <- function(ch3_db, 
                                 table_name = "collapsed_windows",
                                 max_distance = 1000,
                                 sig_cutoff = 0.05,
                                 min_diff = 0.5) 
{
  start_time <- Sys.time()
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  # Check if "mod_diff" table exists
  if (!DBI::dbExistsTable(ch3_db$con, "mod_diff_windows")) {
    stop(glue::glue("Error: Table 'mod_diff_windows' not found in the database. 
                     Please run 'mod_diff()' on windows data first to generate it.\n"))
  }
  
  cat("Collapsing windows on differential analysis results...\n")
  
  # --- Build AVG(...) list dynamically based on existing columns ---
  cols <- DBI::dbListFields(ch3_db$con, "mod_diff_windows")
  patterns <- c("_counts_control$", "_counts_case$", "_frac_control$", "_frac_case$")
  match_cols <- cols[grepl(paste(patterns, collapse = "|"), cols)]
  
  # Quote identifiers safely for DuckDB
  qi <- function(x) as.character(DBI::dbQuoteIdentifier(ch3_db$con, x))
  
  avg_exprs <- if (length(match_cols)) {
    paste0("AVG(", qi(match_cols), ") AS ", qi(paste0("avg_", match_cols)))
  } else character(0)
  
  # This block STARTS with a comma if nonempty; otherwise it's empty.
  avg_select_sql <- if (length(avg_exprs)) {
    paste0(",\n      ", paste(avg_exprs, collapse = ",\n      "))
  } else {
    ""
  }
  
  query <- glue::glue(
    "CREATE OR REPLACE TABLE {table_name} AS
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
     {qi('chrom')},
     MIN(start) AS start,
     MAX({qi('end')}) AS {qi('end')}
     {avg_select_sql}
     , AVG(meth_diff) AS avg_meth_diff
     , COUNT(*) AS num_windows
   FROM RegionGroups
   GROUP BY {qi('chrom')}, region_id
   ORDER BY {qi('chrom')}, start;"
  )
  
  
  DBI::dbExecute(ch3_db$con, query)
  
  end_time <- Sys.time()
  cat("\n")
  total_seconds <- as.numeric(end_time - start_time, units = "secs")
  
  if (total_seconds > 60) {
    total_minutes <- as.numeric(end_time - start_time, units = "mins")
    message("Windows successfully collapsed - ", table_name, " created!",
            "\nTime elapsed: ", round(total_minutes, 2), " minutes\n")
  } else {
    message("Windows successfully collapsed - ", table_name, " created!", 
            "\nTime elapsed: ", round(total_seconds, 2), " seconds\n")
  }
  
  ch3_db$current_table <- table_name
  ch3_db <- .ch3helper_cleanup(ch3_db)
  
  invisible(ch3_db)
}
