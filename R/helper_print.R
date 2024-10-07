#' Helper Function to Print Database Tables
#'
#' This internal function connects to a DuckDB database and prints the names of all tables currently in the 
#' database. It can also print the first few rows of specific tables.
#'
#' @param ch3_db A character string or an object of class `ch3_db` representing the DuckDB database to connect to.
#' @param tables A character vector specifying the names of tables to print. Defaults to "all", which 
#' will print all tables. If specific table names are provided, it will only print those.
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
#' @keywords internal

.helper_print <- function(ch3_db, tables = "all") {
  # Open the database connection
  db_con <- .helper_connectDB(ch3_db)
  
  tryCatch(
    {
      # List all tables in the database
      all_tables <- dbListTables(db_con)
      cat("All tables currently in database:\n")
      
      # Print each table name
      for (tbl_name in all_tables) {
        cat(tbl_name, "\n")
      }
      
      if (tables == "all") {
        tables <- all_tables
      }
      
      for (tb_name in tables) {
        if (tb_name %in% all_tables) {
          # Print the table name
          cat("Table: ", tb_name, "\n")
          
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
      .helper_closeDB(ch3_db, db_con)
    }
  )
}
