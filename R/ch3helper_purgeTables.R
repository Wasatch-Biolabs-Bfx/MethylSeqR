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
#' @importFrom DBI dbListTables dbRemoveTable
#'
#' @keywords internal

.ch3helper_purgeTables <- function(db_con)
{
  # List all tables in the database
  all_tables <- dbListTables(db_con)
  
  # Remove tables that start with "temp"
  for (table in all_tables) {
    if (startsWith(table, "temp")) {
      dbRemoveTable(db_con, table)
    }
  }
  
  return(db_con)
}