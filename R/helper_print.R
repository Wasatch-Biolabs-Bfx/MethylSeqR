#' Helper Function to Print Database Tables
#'
#' This function connects to a DuckDB database and prints the names of all tables currently in the 
#' database. It can also print the first few rows of specific tables.
#'
#' @param ch3_db A character string or an object of class `ch3_db` representing the DuckDB database to connect to.
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
#' @import DBI
#' 
#' @examples 
#' # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Print out tables in the database
#'  print_db(ch3_db, tables = "all")
#'
#' @export

print_db <- function(ch3_db, tables = "") {
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  tryCatch(
    {
      # List all tables in the database
      cat("Last table added: ")
      cat(database$last_table)
      cat("\n")
      cat("\n")
      all_tables <- dbListTables(db_con)
      cat("All tables currently in database:\n")
      
      # Print each table name
      for (tbl_name in all_tables) {
        cat(tbl_name, "\n")
      }
      cat("\n")
      
      if (tables == "all") {
        tables <- all_tables
      } else if (tables == "") {
        tables <- database$last_table
      }
      
      for (tb_name in tables) {
        if (tb_name %in% all_tables) {
          # Print the table name
          cat("Printing table: ", tb_name, "\n")
          
          # Fetch the first few rows and print them
          print(head(tbl(db_con, tb_name)))
          
        } else {
          # Print a message if the table does not exist
          message(paste0("Table '", tb_name, "' does not exist in the database."))
        }
      }
    }, 
    error = function(e) {
      # Print custom error message
      message("An error occurred: ", e$message)
    }, 
    finally = {
      .helper_closeDB(database)
    }
  )
}
