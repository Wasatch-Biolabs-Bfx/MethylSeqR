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
#' @return A database connection object.
#' 
#' @importFrom DBI dbConnect dbListTables
#' @importFrom duckdb duckdb
#'
#' @keywords internal
#' 

.ch3helper_connectDB <- function(ch3_db)
{
  # Open the database connection
  if (inherits(ch3_db, "character")) { # if given a string file name- create a new database
    # Open the database connection - first check to make sure correct name is there
    if (is.character(ch3_db)) {
      if (!grepl(".ch3.db$", ch3_db)) {
        ch3_db <- paste0(ch3_db, ".ch3.db")
      }
    }
    # make ch3_object
    database <- list(db_file = ch3_db, current_table = NULL, con = NULL)
    class(database) <- "ch3_db"
    
    #add in the connection
    database$con <- dbConnect(duckdb(database$db_file), read_only = FALSE)
    defer(.ch3helper_closeDB(database$con), parent.frame())
    
    # return database object
    return(database)
  
  } else if (inherits(ch3_db, "ch3_db")) { # if given a ch3_db OBJECT!
    ch3_db$con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
    defer(.ch3helper_closeDB(ch3_db), parent.frame())
    
    return(ch3_db)
    
  } else {
    stop("Invalid ch3_db class. Must be character or ch3_db.")
  }
  
}