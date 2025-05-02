#' Get Column Names from a DuckDB Table
#'
#' Returns a character vector of column names from a specified table
#' in a DuckDB database file.
#'
#' @param ch3_db Path to the DuckDB database file (e.g., `"my_data.ch3.db"`).
#' @param table_name Name of the table whose column names are to be retrieved.
#'
#' @return A character vector of column names from the specified table.
#'
#' @examples
#' \dontrun{
#' get_ch3_cols("my_data.ch3.db", "windows")
#' }
#'
#' @import DBI
#' @import duckdb
#' @export
get_ch3_cols <- function(ch3_db, table_name) {
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  # Check if table exists
  if (!(table_name %in% dbListTables(ch3_db$con))) {
    stop(paste("Table", table_name, "does not exist in the database."))
  }
  
  cat(paste0("Columns in ", table_name, " table:\n"))
  print(dbListFields(ch3_db$con, table_name))
  
  ch3_db <- .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}