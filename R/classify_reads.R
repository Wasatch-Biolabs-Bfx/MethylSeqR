classify_reads <- function(ch3_db,
                           key_table,
                           case,
                           control,
                           meth_diff_threshold = 0.1) {
  
  database <- MethylSeqR:::.helper_connectDB(ch3_db)
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
  
  query <- glue("
  CREATE TABLE classified_reads AS
  SELECT
    r.sample_name,
    r.read_id,
    r.first_cpg_pos,
    r.last_cpg_pos,
    r.mh_frac,
    CASE
      WHEN ABS(r.mh_frac - k.avg_mh_frac_control) <= 0.1 THEN '{control}'
      WHEN ABS(r.mh_frac - k.avg_mh_frac_case) <= 0.1 THEN '{case}'
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