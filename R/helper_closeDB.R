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
#' @return None. This function is called for its side effects (updating the `ch3_db` object and closing the connection).
#'
#' @examples
#' # Assuming `ch3_db` is an object and `db_con` is a valid database connection
#' helper_closeDB(ch3_db, db_con)
#'
#' @export
.helper_closeDB <- function(ch3_db, db_con)
{
  # Finish up: update table list and close the connection
  ch3_db$tables <- dbListTables(db_con)
  dbDisconnect(db_con, shutdown = TRUE)
}