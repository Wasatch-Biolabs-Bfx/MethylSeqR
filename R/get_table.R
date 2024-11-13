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
#' @import dplyr
#' @import duckdb
#'
#' @examples
#' # Assuming ch3_db is a valid database object and "positions" is a table in the database
#' ch3_db <- system.file("my_data.ch3.db", package = "MethylseqR")
#' positions = get_table(ch3_db, "positions")
#'
#' @export
get_table <- function(ch3_db, 
                      table_name)
{
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  tryCatch(
    {
      dat <- tbl(db_con, table_name) |>
        collect()
    }, 
    error = function(e) 
    {
      message("Error: ", e$message)
    },
    finally = 
      {
        dbDisconnect(db_con, shutdown = TRUE)
        
        return(dat)
      })
  
}