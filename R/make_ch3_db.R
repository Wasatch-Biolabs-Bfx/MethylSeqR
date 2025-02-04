make_ch3_db <- function(ch3_files, 
                        db_name,
                        chrom = NULL, 
                        min_read_length = 50, 
                        min_call_prob = 0.9, 
                        min_base_qual = 30, 
                        flag = NULL)
{
  # Check if folder/files exist and create string for query
  if (any(!file.exists(ch3_files))) {
    stop("One or more file/folder names do not exist.")
  }
  
  ch3_files[dir.exists(ch3_files)] <- 
    paste0(ch3_files[dir.exists(ch3_files)], "/*")
  path <- paste0("'", ch3_files, "'", collapse = ", ")
  
  # Setup files and db
  if (!grepl(".ch3.db$", db_name)) {
    db_name <- paste0(db_name, ".ch3.db")
  } 
  
  ch3_db <- list(db_file = db_name)
  class(ch3_db) <- "ch3_db"
  db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
  ch3_db$tables <- dbListTables(db_con)
  
  on.exit(dbDisconnect(db_con, shutdown = TRUE))
  
  # Check if the table already exists and delete it if it does
  for (table in dbListTables(db_con)) {
    dbExecute(db_con, paste0("DROP TABLE ", table))
  }
  
  # Start building the WHERE clause
  filters <- c()
  
  # Filter chromosome if specified
  if (!is.null(chrom)) {
    filters <- c(filters, paste0("chrom = '", chrom, "'"))
  }
  
  # Filter minimum read length
  filters <- c(filters, paste0("read_length >= ", min_read_length))
  
  # Filter minimum call probability
  filters <- c(filters, paste0("call_prob >= ", min_call_prob))
  
  # Filter minimum base quality
  filters <- c(filters, paste0("base_qual >= ", min_base_qual))
  
  # Filter flag if specified
  if (!is.null(flag)) {
    filters <- c(filters, paste0("flag = ", flag))
  }
  
  # Combine all filters into a WHERE clause
  where_clause <- 
    if (length(filters) > 0) 
      paste("WHERE", paste(filters, collapse = " AND ")) else ""
  
  # Execute the query
  dbExecute(db_con, 
            paste0("CREATE TABLE calls AS SELECT * FROM read_parquet([", 
                   path, "]) ", 
                   where_clause)) 
  
  message("Database successfully created at ", ch3_db$db_file)
  ch3_db$tables <- dbListTables(db_con)
  invisible(ch3_db)
}