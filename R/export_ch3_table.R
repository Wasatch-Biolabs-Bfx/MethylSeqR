#' Export Tables from the ch3 Database
#'
#' This function exports specified tables from the ch3 database to CSV files. Can export one or multiple tables as a time. It checks whether each table exists in the database 
#' before exporting, and provides informative messages for any missing tables. The output CSV files are saved at the specified path.
#'
#' @param ch3_db A string. The path to the database containing ch3 files from nanopore data.
#' @param table A character vector specifying the table to be exported from the database. Default is "positions".
#' @param out_path A string. The path to the directory where the CSV files will be saved. The file will automatically be named "{table name}.csv". 
#'
#' @details
#' The function connects to the specified database and iterates through the list of table names provided in the `tables` parameter. 
#' For each table that exists in the database, it reads the table into R and writes it as a CSV file to the location specified by `out_path`. 
#' If a table does not exist in the database, a message is printed indicating this.
#'
#' In case of any error during the execution, a custom error message is displayed. The function ensures that the database connection 
#' is closed safely using the `finally` block.
#'
#' @note The function assumes that the tables specified in `tables` exist in the database and can be accessed via the `DBI` package.
#' 
#' @return NULL. The function writes the specified tables to CSV files.
#' 
#' @importFrom DBI dbConnect dbDisconnect dbListTables
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl
#' @importFrom utils write.csv
#'
#' @export

export_ch3_table <- function(ch3_db, 
                             table = "positions", 
                             out_path) {
  # Open DB connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  on.exit({ .ch3helper_closeDB(ch3_db) }, add = TRUE)
  
  # Build output path
  if (!grepl("\\.csv$", out_path, ignore.case = TRUE)) {
    if (!grepl("/$", out_path)) out_path <- paste0(out_path, "/")
    out_path <- paste0(out_path, table, ".csv")
  }
  
  if (!DBI::dbExistsTable(ch3_db$con, table)) {
    message("Table '", table, "' does not exist.")
    return(invisible(ch3_db))
  }
  
  tbl_id  <- as.character(DBI::dbQuoteIdentifier(ch3_db$con, table))
  file_q  <- as.character(DBI::dbQuoteString(ch3_db$con, out_path))
  
  # COPY is executed entirely inside DuckDB (fast + avoids R memory)
  sql <- sprintf(
    "COPY (SELECT * FROM %s) TO %s WITH (FORMAT CSV, HEADER, DELIMITER ',', QUOTE '\"');",
    tbl_id, file_q
  )
  
  message("Exporting '", table, "' to: ", out_path)
  n <- DBI::dbExecute(ch3_db$con, sql)
  message("Rows written: ", n)
  
  invisible(ch3_db)
}
