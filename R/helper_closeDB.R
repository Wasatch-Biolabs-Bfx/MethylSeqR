#' Close Database Connection
#'
#' This function updates the table list in the `ch3_db` object and closes the database connection.
#'
#' @param ch3_db An object containing information about the database, including the list of tables.
#' @param db_con A database connection object. This should be a connection created by a database connection function.
#'
#' @details
#' The function updates the `tables` attribute of the `ch3_db` object with the current list of tables in the connected database 
#' before closing the connection.
#'
#' @import DBI
#'
#' @return None. This function is called for its side effects (updating the object and closing the connection).
#'
#' @export
.helper_closeDB <- function(ch3_db)
{
  # add this on every function for CONSISTENCY
  # Finish up: update table list and close the connection
  ch3_db$tables <- dbListTables(ch3_db$db_con)
  
  dbDisconnect(ch3_db$db_con, shutdown = TRUE)
  ch3_db$db_con = NULL
  
  return(ch3_db)
}