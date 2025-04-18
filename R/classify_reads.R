#' Classify Reads as Case or Control Based on Methylation Profiles
#'
#' This function classifies reads in the `reads` table of a DuckDB database
#' as either `"case"`, `"control"`, or `"unknown"` based on similarity
#' to reference methylation fractions in a `key_table` (e.g., a collapsed windows file).
#' Classification is based on how close the read's `mh_frac` is to the average
#' case or control `mh_frac`, within a user-defined `meth_diff_threshold`.
#'
#' @param ch3_db Path to a DuckDB `.db` file created by this package (e.g., from `summarize_reads()`).
#' @param key_table Path to a CSV, TSV, or BED file generated by `collapse_windows()`. Must include the columns:
#' `chrom`, `start`, `end`, `avg_mh_frac_control`, `avg_mh_frac_case`, and `avg_meth_diff`.
#' @param case Character string used to label case reads (e.g., `"case"`).
#' @param control Character string used to label control reads (e.g., `"control"`).
#' @param meth_diff_threshold Numeric value specifying the maximum difference in `mh_frac` allowed
#' to match either case or control averages. Must be less than half the minimum absolute value of
#' `avg_meth_diff` to prevent ambiguous classifications.
#'
#' @return Invisibly returns the open database connection with a new table named `classified_reads`
#' added to the database. This table includes:
#' \itemize{
#'   \item \code{sample_name} – Sample identifier
#'   \item \code{read_id} – Unique read identifier
#'   \item \code{first_cpg_pos} – First CpG position of the read
#'   \item \code{last_cpg_pos} – Last CpG position of the read
#'   \item \code{mh_frac} – Methylation fraction of the read
#'   \item \code{classification} – `"case"`, `"control"`, or `"unknown"`
#' }
#'
#' @details 
#' This function runs entirely in SQL for scalability. It performs an interval join
#' between the `reads` table and the key table on `chrom` and CpG position range,
#' then classifies each read based on proximity of `mh_frac` to either `avg_mh_frac_control` or `avg_mh_frac_case`.
#'
#' @examples
#' \dontrun{
#' classify_reads(
#'   ch3_db = "my_data.ch3.db",
#'   key_table = "key_table.csv",
#'   case = "treated",
#'   control = "untreated",
#'   meth_diff_threshold = 0.1
#' )
#' }
#'
#' @export

classify_reads <- function(ch3_db,
                           key_table,
                           case,
                           control,
                           meth_diff_threshold = 0.1) {
  
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Check if "mod_diff" table exists
  if (!DBI::dbExistsTable(db_con, "reads")) {
    stop(glue::glue("Error: Table 'reads' not found in the database. 
                     Please run 'summarize_reads()' on windows data first to generate it."))
  }
  
  on.exit(.helper_purgeTables(db_con), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(.helper_closeDB(database), add = TRUE)
  
  # Read in key_table and make sure the key_table looks like collapsed_windows...
  file_ext <- tools::file_ext(key_table)
  
  if (file_ext == "csv") {
    annotation <- read_csv(key_table,
                           show_col_types = FALSE)
  } else if (file_ext %in% c("bed", "tsv")) {
    annotation <- read_tsv(key_table,
                           show_col_types = FALSE)
  } else {
    stop("Invalid file type. Only CSV, TSV, or BED files are supported.")
  }
  
  # Make sure key_table is a collapsed window format from collapse_windows().
  required_cols <- c("chrom", "start", "end", "avg_mh_frac_control", "avg_mh_frac_case", "avg_meth_diff")
  
  if (!all(required_cols %in% colnames(annotation))) {
    stop("\nMissing one or more required columns: start, end, avg_mh_frac_control, avg_meth_diff.\nkey_table must be a collapsed_windows table from the function collapse_windows().")
  }
  
  
  # Check to make sure threshold, make sure it is less than 1/2 of min_diff... make sure threshold doesn't overlap
  if (any(meth_diff_threshold >= abs(annotation$avg_meth_diff) / 2, na.rm = TRUE)) {
    stop("`\nmeth_diff_threshold` must be less than half the absolute value of all `avg_meth_diff` values.\nThis is to avoid overlapping in classification.\n")
  }
  
  # Upload annotation as a temporary table
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_key_table;")
  dbWriteTable(db_con, "temp_key_table", annotation, temporary = TRUE)
  
  dbExecute(db_con, "DROP TABLE IF EXISTS classified_reads;")
  
  query <- glue("
  CREATE TABLE classified_reads AS
  SELECT
    r.sample_name,
    r.read_id,
    r.first_cpg_pos,
    r.last_cpg_pos,
    r.mh_frac,
    CASE
      WHEN ABS(r.mh_frac - k.avg_mh_frac_control) <= {meth_diff_threshold} THEN '{control}'
      WHEN ABS(r.mh_frac - k.avg_mh_frac_case) <= {meth_diff_threshold} THEN '{case}'
      ELSE 'unknown'
    END AS classification
  FROM
    reads r
  JOIN
    temp_key_table k
  ON
    r.chrom = k.chrom
    AND r.first_cpg_pos BETWEEN k.start AND k.end;
  ")
  
  dbExecute(db_con, query)
  
  # Drop temporary tables
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_key_table;")
  
  cat("\n")
  message("classified_reads table successfully created!")
  print(head(tbl(db_con, "classified_reads")))
  
  invisible(database)
}