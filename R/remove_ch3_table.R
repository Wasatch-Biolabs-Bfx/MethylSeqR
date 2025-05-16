#' Remove a table from a database connection.
#'
#' This function checks if a table exists and, if so, removes it.
#' It provides messages about the outcome.
#'
#' @param ch3_db A database object.
#' @param table_name A string specifying the name of the table to remove.
#'
#' @return Invisibly returns the database with the removed table.
#'
#' @importFrom DBI dbExistsTable dbRemoveTable
#' @importFrom glue glue
#'
#' @export
remove_ch3_table <- function(ch3_db, table_name) 
{
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  if (!is.character(table_name) || length(table_name) != 1 || nchar(trimws(table_name)) == 0) {
    stop("Error: 'table_name' must be a non-empty single string.")
  }
  
  if (dbExistsTable(ch3_db$con, table_name)) {
    dbRemoveTable(ch3_db$con, table_name)
    message(glue("Table '{table_name}' successfully removed from the database."))
    } 
  else {
    message(glue("Table '{table_name}' does not exist in the database. No action taken."))
  }
  
  ch3_db <- .ch3helper_cleanup(ch3_db)
  invisible(ch3_db)
}