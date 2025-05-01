#' Create a CH3 Database from Parquet Files
#'
#' This function processes a collection of CH3 files stored in Parquet format and creates 
#' a DuckDB database containing filtered and structured methylation call data.
#'
#' @param ch3_files A character vector of file or directory paths containing CH3 Parquet files. 
#' If a directory is provided, all files within it will be processed.
#' @param db_name A string representing the path where the CH3 database will be created. 
#' The extension ".ch3.db" will be appended if not already present.
#' @param chrom An optional character string specifying a chromosome to filter the data. 
#' If NULL, all chromosomes will be included.
#' @param min_read_length A numeric value specifying the minimum read length required 
#' for inclusion in the database. Default is 50.
#' @param min_call_prob A numeric value representing the minimum call probability threshold. 
#' Only calls with a probability greater than or equal to this value will be included. Default is 0.9.
#' @param min_base_qual A numeric value representing the minimum base quality threshold. 
#' Only reads with quality scores at or above this value will be included. Default is 30.
#' @param flag An optional numeric value specifying a flag-based filter for the data. 
#' If NULL, no flag filtering is applied. 
#' @param chr_prefix A boolean value stating whether or not to keep the chr prefix 
#' in the chromosome column. Default is TRUE, and "chr" will be kept in front of every
#' chr number.
#'
#' @details
#' This function reads CH3 files stored in Parquet format and imports them into a DuckDB database. 
#' The data is filtered based on user-specified criteria, including chromosome, read length, 
#' call probability, base quality, and flag values. If a table already exists in the database, 
#' it will be dropped and recreated.
#'
#' @return A list of class `"ch3_db"`, containing:
#' \item{db_file}{The path to the created DuckDB database file.}
#' \item{tables}{A list of tables available in the database.}
#'
#'
#' @examples
#' # Example usage
#' ch3_files <- system.file("new_test_data", package = "MethylSeqR")
#' db_name <- tempfile("example_ch3")
#' 
#' make_ch3_db(ch3_files, db_name)
#'
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbListTables
#' @importFrom duckdb duckdb
#' 
#' @export

make_ch3_db <- function(ch3_files, 
                        db_name,
                        chrom = NULL, 
                        min_read_length = 50, 
                        min_call_prob = 0.9, 
                        min_base_qual = 10, 
                        flag = NULL,
                        chr_prefix = TRUE)
{
  # Check if folder/files exist and create string for query
  if (any(!file.exists(ch3_files))) {
    stop("One or more file/folder names do not exist.")
  }
  
  ch3_files[dir.exists(ch3_files)] <- paste0(ch3_files[dir.exists(ch3_files)],
                                             "/*.ch3")
  
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
  if (chr_prefix) {
    dbExecute(db_con, 
              paste0("CREATE TABLE calls AS 
                   SELECT *, 
                          CASE WHEN chrom LIKE 'chr%' THEN chrom 
                                    ELSE 'chr' || chrom 
                               END AS chrom_fixed 
                   FROM read_parquet([", 
                   path, "]) ", 
                   where_clause)) 
  } else {
    dbExecute(db_con, 
              paste0("CREATE TABLE calls AS 
                   SELECT *, 
                          CASE WHEN chrom LIKE 'chr%' THEN SUBSTRING(chrom, 4) 
                              ELSE chrom 
                         END AS chrom_fixed 
                   FROM read_parquet([", 
                   path, "]) ", 
                   where_clause)) 
    dbExecute(db_con, "ALTER TABLE calls DROP COLUMN chrom")
    dbExecute(db_con, "ALTER TABLE calls RENAME COLUMN chrom_fixed TO chrom")
  }
  
  message("Database successfully created at ", ch3_db$db_file)
  ch3_db$tables <- dbListTables(db_con)
  
  # message(paste0("Samples include: ", 
  #                unique(pull(tbl(db_con, "calls"), sample_name))))
  
  invisible(ch3_db)
}