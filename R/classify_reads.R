classify_reads <- function(ch3_db,
                           key_table,
                           case,
                           control) {
  
  database <- MethylSeqR:::.helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Check if "mod_diff" table exists
  if (!DBI::dbExistsTable(db_con, "reads")) {
    stop(glue::glue("Error: Table 'reads' not found in the database. 
                     Please run 'summarize_reads()' on windows data first to generate it."))
  }
  
  on.exit(MethylSeqR:::.helper_purgeTables(db_con), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(MethylSeqR:::.helper_closeDB(database), add = TRUE)
  
  # make sure the key_table looks like collapsed_windows...
  
  
}