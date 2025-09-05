#' Helper Function to Print Database Tables
#'
#' This function connects to a DuckDB database and prints the names of all tables currently in the 
#' database. It can also print the first few rows of specific tables.
#'
#' @param ch3_db A character string or an object of class ch3_db representing the DuckDB database to connect to.
#' @param tables A character vector specifying the names of tables to print. Defaults to the last table given, which 
#' will print all tables. If specific table names are provided, it will only print those (ex. "positions", "regions", "meth_diff").
#' If tables ="all", the function will print out all table sin the database...
#'
#' @details
#' The function establishes a connection to the database, retrieves the list of tables, and prints their 
#' names. If the specified table exists, it also prints the first few rows of that table. If a specified 
#' table does not exist, a message is printed to indicate this.
#'
#' @return None. This function is called for its side effects (printing information).
#' 
#' @examples 
#' # Specify the path to the database
#'  ch3_files <- system.file("test_data", package = "MethylSeqR")
#'  ch3_db <- tempfile("example_db")
#'  
#'  # Print out tables in the database
#'  make_pos_db(ch3_files, ch3_db) |> print()
#'  
#' @importFrom DBI dbListTables
#' @importFrom dplyr tbl
#'
#' @export

print.ch3_db <- function(ch3_db) {
  cat("<ch3_db object>\n")
  
  cat("Database file:\n")
  cat("  ", ch3_db$db_file, "\n")
  
  cat("Current table:\n")
  cat("  ", ifelse(is.null(ch3_db$current_table), "NULL", ch3_db$current_table), "\n")
  
  cat("Connection:\n")
  
  if (is.character(ch3_db$con) && ch3_db$con == "none") {
    cat("  NULL\n")
  } else if (!dbIsValid(ch3_db$con))
  {
    cat("  NULL\n")
  } else
  {
    cat("  Active DBI connection\n")
  }
}