#' Export Tables from the ch3 Database
#'
#' This function exports specified tables from the ch3 database to CSV files. It checks whether each table exists in the database 
#' before exporting, and provides informative messages for any missing tables. The output CSV files are saved at the specified path.
#'
#' @param ch3_db A string. The path to the database containing ch3 files from nanopore data.
#' @param tables A character vector specifying the tables to be exported from the database. Default is "positions".
#' @param out_path A string. The path where the CSV files will be saved.
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
#' @importFrom DBI dbListTables
#' @importFrom dplyr tbl
#' @import utils
#' 
#' @return NULL. The function writes the specified tables to CSV files.
#'
#' @export
export_tables <- function(ch3_db,
                         tables = c("positions"),
                        out_path) 
{
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  tryCatch(
    {
      all_tables <- dbListTables(db_con)
      
      for (tbl in tables) {
        if (tbl %in% all_tables) {
          # Print the table name
          message(paste0("Writing out ", tbl, " table..."))
          # write out table to path given
          write.csv(tbl(db_con, tbl), file = out_path, row.names = FALSE, quote = FALSE)
          
          
        } else {
          # Print a message if the table does not exist
          message(paste0("Table '", tbl, "' does not exist in the database."))
        }
      }
    }, 
    error = function(e) 
    {
      # Print custom error message
      message("An error occurred: ", e$message)
    }, 
    finally = 
      {
        # Finish Up
        .helper_closeDB(database)
        return(database)
      })
}