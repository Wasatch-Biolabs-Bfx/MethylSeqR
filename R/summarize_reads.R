summarize_reads <- function(ch3_db,
                            key_table = NULL,
                            min_CGs = 5) {
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # Finish up: purge extra tables & update table list and close the connection
  keep_tables = c("calls", "positions", "regions", "windows", 
                  "mod_diff_positions", "mod_diff_regions", "mod_diff_windows",
                  "collapsed_windows", "filtered_calls", "reads")
  on.exit(.helper_purgeTables(db_con, keep_tables), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(.helper_closeDB(database), add = TRUE)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
  
  # First - filter calls for region
  # Create the filtered_calls table
  dbExecute(db_con, "DROP TABLE IF EXISTS filtered_calls;")
  
  dbExecute(db_con, "
    CREATE TABLE filtered_calls AS 
    SELECT calls.*
    FROM calls
    JOIN collapsed_windows 
    ON calls.chrom = collapsed_windows.chrom
    AND calls.start >= collapsed_windows.start
    AND calls.end <= collapsed_windows.end;
")
  
  # Second - summarize reads using the filtered calls table
  # Included: filter total_calls by the min_CGs parameter wanted!
  dbExecute(db_con, "DROP TABLE IF EXISTS reads;")
  
  query <- glue("
    CREATE TABLE reads AS 
    SELECT
        ANY_VALUE(sample_name) AS sample_name,
        read_id,
        COUNT(*) AS total_calls,
        MIN(start) AS first_cpg_pos,
        MAX(\"end\") AS last_cpg_pos,
        ANY_VALUE(read_length) AS read_length,
        SUM(CASE WHEN call_code = '-' THEN 1 ELSE 0 END) AS c_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) AS m_counts,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) AS h_counts,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) AS mh_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS m_frac,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS h_frac,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS mh_frac
    FROM filtered_calls
    GROUP BY read_id
    HAVING total_calls >= {min_CGs};
")
  
  dbExecute(db_con, query)
  
  # Then make classify_reads() function
  cat("\n")
  message("reads table successfully created!")
  print(head(tbl(db_con, "reads")))
  
  invisible(database)
  
}