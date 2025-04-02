#' Internal Function to Purge Unwanted Tables from Database
#'
#' This internal function connects to a DuckDB database and removes tables that are not specified in the 
#' `keep_tables` list. It retains only the tables that you want to keep in the database.
#'
#' @param ch3_db A character string or an object of class `ch3_db` representing the DuckDB database to connect to.
#'
#' @details
#' The function connects to the specified database, lists all tables, and removes those not included in 
#' the `keep_tables` vector. After purging, it prints the names of the remaining tables in the database.
#'
#' @return None. This function is called for its side effects (modifying the database).
#'
#' @import DBI
#'
#' @keywords internal
.helper_purgeTables <- function(db_con)
{
  # keep_tables involves any possible table created for the database.
  keep_tables = c("calls", "positions", "regions", "windows", 
                  "mod_diff_positions", "mod_diff_regions", "mod_diff_windows",
                  "collapsed_windows", "reads")
  tryCatch(
    {
      # List all tables in the database
      all_tables <- dbListTables(db_con)
      
      # Remove tables that are not in the 'keep_tables' list
      for (table in all_tables) {
        if (!(table %in% keep_tables)) {
          # message(paste0("Removing table: ", table))
          dbRemoveTable(db_con, table)
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
        return(db_con)
      })
}