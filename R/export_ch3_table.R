#' Export Tables from the ch3 Database
#'
#' This function exports specified tables from the ch3 database to CSV files. Can export one or multiple tables as a time. It checks whether each table exists in the database 
#' before exporting, and provides informative messages for any missing tables. The output CSV files are saved at the specified path.
#'
#' @param ch3_db A string. The path to the database containing ch3 files from nanopore data.
#' @param table A character vector specifying the table to be exported from the database. Default is "positions".
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
                        out_path) 
{
  # Open the database connection
  database <- .ch3helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  on.exit(.ch3helper_closeDB(database), add = TRUE)
  
  # If the out_path does not end with ".csv", append "_<tablename>.csv"
  if (!grepl("\\.csv$", out_path)) {
    
    # If out_path is a directory, make sure it ends with a "/"
    if (grepl("/$", out_path) == FALSE) {
      out_path <- paste0(out_path, "/")
    }
    
    out_path <- paste0(out_path, "_", table, ".csv")
  }
  
  all_tables <- dbListTables(db_con)
  
  if (table %in% all_tables) {
      # Print the table name
      message(paste0("Writing out ", table, " table to ", out_path))
      # Specify the file path
      
    #  FUTURE CODE TO IMPLEMENT...
    # dbExecute(db_con, 
    #           paste("COPY (SELECT * FROM", table, 
    #                 ") TO '", out_path, 
    #                 "' WITH (FORMAT CSV, HEADER TRUE);", sep = ""))
    #   
    #   
      # write out table to path given
      write.csv(tbl(db_con, table),
                file = out_path,
                row.names = FALSE,
                quote = FALSE)
    } else {
      # Print a message if the table does not exist
      message(paste0("Table '", tbl, "' does not exist in the database."))
    }
  
  invisible(database)
}