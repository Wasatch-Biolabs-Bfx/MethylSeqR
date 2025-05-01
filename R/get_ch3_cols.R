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
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ch3_db)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  
  DBI::dbListFields(con, table_name)
}