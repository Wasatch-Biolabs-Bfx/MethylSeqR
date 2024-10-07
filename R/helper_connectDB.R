#' Connect to a Database
#'
#' This internal function establishes a connection to a DuckDB database. It can handle both a character file name 
#' or an object of class `ch3_db` to open the database.
#'
#' @param ch3_db A character string representing the file path to the DuckDB database or an object of class `ch3_db`.
#'
#' @details
#' This function checks the class of `ch3_db` and attempts to connect to the database. If `ch3_db` is a character string, 
#' it will create an object of class `ch3_db`. If `ch3_db` is already of class `ch3_db`, it will directly establish a 
#' connection to the database.
#'
#' @note This function is intended for internal use within the package.
#'
#' @import DBI duckdb
#' 
#' @return A database connection object.
#'
#' @examples
#' # db_con <- .helper_connectDB("path/to/database.db")
#' # db_con <- .helper_connectDB(ch3_db_object)
#'
#' @keywords internal
#' 
.helper_connectDB <- function(ch3_db)
{
  # Open the database connection
  if (inherits(ch3_db, "character")) { # if given a string file name
    ch3_db <- list(db_file = ch3_db)
    class(ch3_db) <- "ch3_db"
    db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
    ch3_db$tables <- dbListTables(db_con)
  } else if (inherits(ch3_db, "ch3_db")) { # if given a 
    db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
  } else {
    stop("Invalid ch3_db class. Must be character or ch3_db.")
  }
  
  return(db_con)
}