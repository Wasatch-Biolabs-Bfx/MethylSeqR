#' Close Database Connection
#'
#' This function updates the table list in the `ch3_db` object and closes the database connection.
#'
#' @param ch3_db An object containing information about the database, including the list of tables.
#'
#' @details
#' The function updates the `tables` attribute of the `ch3_db` object with the current list of tables in the connected database 
#' before closing the connection.
#'
#' @return None. This function is called for its side effects (updating the object and closing the connection).
#'
#' @importFrom DBI dbDisconnect dbListTables dbIsValid
#'
#' @keywords internal

.ch3helper_closeDB <- function(ch3_db)
{
  if (is.character(ch3_db$con) && ch3_db$con == "none") {
    return(NULL)
  }
  
  if (dbIsValid(ch3_db$con)) {
    dbDisconnect(ch3_db$con, shutdown = TRUE)
    ch3_db$con <- "none"
  }
  
  return(ch3_db)
}