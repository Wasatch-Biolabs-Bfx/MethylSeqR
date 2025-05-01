#' Collect Table from DuckDB Database as Tibble
#'
#' This function connects to a DuckDB database and collects a specified table as a tibble.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param table_name A string representing the name of the table to collect from the database.
#' 
#' @details
#' The function establishes a connection to the DuckDB database using \code{.helper_connectDB}.
#' It retrieves the specified table as a tibble. If an error occurs during table retrieval, 
#' a message with the error is displayed. The database connection is closed after retrieving 
#' the data, regardless of success or failure.
#'
#' @return A tibble containing the collected data from the specified database table. If the table retrieval fails, an empty tibble is returned.
#'
#' @examples
#' # Assuming ch3_db is a valid database object and "positions" is a table in the database
#' ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#' positions = get_ch3_table(ch3_db, "positions")
#'
#' @importFrom DBI dbConnect dbDisconnect dbListTables
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl collect
#' @importFrom tibble tibble
#' 
#' @export

get_ch3_table <- function(ch3_db, 
                      table_name)
{
  # Open the database connection
  database <- .ch3helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  on.exit(.ch3helper_closeDB(database), add = TRUE)
  
  dat <- tibble()  # Initialize an empty tibble to return if there's an error
  
  if (table_name %in% dbListTables(db_con)) {
    # write out table to path given
    dat <- tbl(db_con, table_name) |> collect()
  } else {
    # Print a message if the table does not exist
    message(paste0("Table '", table_name, "' does not exist in the database."))
  }
  
  return(dat)
}