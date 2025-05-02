#' Collect Table from DuckDB Database as Tibble
#'
#' This function connects to a DuckDB database and collects a specified table as a tibble.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param table_name A string representing the name of the table to collect from the database.
#' @param max_rows The maximum amount of rows wanted for calculation. This argument can help analysis run faster when there is a lot of data.
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
#' ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#' positions = get_ch3_table(ch3_db, "positions")
#'
#' @export
get_ch3_table <- function(ch3_db, 
                      table_name,
                      max_rows = NULL)
{
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  dat <- tibble()  # Initialize an empty tibble to return if there's an error
  
  if (table_name %in% dbListTables(ch3_db$con)) {
    # Reference the table
    table_ref <- tbl(ch3_db$con, table_name)
    
    if (!is.null(max_rows)) {
      # Randomly sample max_rows rows (efficient with DuckDB)
      dat <- table_ref |> slice_sample(n = max_rows) |> collect()
    } else {
      dat <- table_ref |> collect()
    }
  } else {
    message(paste0("Table '", table_name, "' does not exist in the database."))
  }
  
  return(dat)
}