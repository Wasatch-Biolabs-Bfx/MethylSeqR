#' Count unique CpG sites in a CH3 database table
#'
#' Computes the number of unique CpG sites (defined by distinct pairs of
#' \code{start} and \code{end} positions) in a table of calls within a CH3
#' database. Accepts either a file path to a database or an existing CH3
#' database object.
#'
#' @param ch3_db Character string giving the path to a CH3 database file,
#'   or an object/handle accepted by \code{.ch3helper_connectDB()}.
#' @param table_name Character scalar. Name of the table to query. Defaults to
#'   \code{"calls"}.
#'
#' @return An integer (length-one) giving the number of unique CpG sites,
#'   returned \emph{invisibly}. If \code{table_name} does not exist, returns
#'   \code{NULL} invisibly after printing a message. The function also prints a
#'   formatted summary to the console as a side effect.
#'   
#' @section Errors:
#' If \code{ch3_db} is a character path and the file does not exist, an error is
#' thrown with \code{stop()} before attempting any connection.
#'
#' @examples
#' \dontrun{
#' # Using a file path
#' n <- get_ch3_cpg_count("/path/to/ch3.sqlite", table_name = "calls")
#' n
#'
#' # Using a pre-opened handle (package-internal)
#' dbh <- .ch3helper_connectDB("/path/to/ch3.sqlite")
#' get_ch3_cpg_count(dbh)              # prints summary, returns value invisibly
#' }
#'
#' @importFrom DBI dbListTables dbGetQuery
#' @export
get_ch3_cpg_count <- function(ch3_db, table_name = "calls") {
  # Check if database path exists if passed as a string
  if (is.character(ch3_db)) {
    if (!file.exists(ch3_db)) {
      stop(paste("The database file", ch3_db, "does not exist."))
    }
  }
  
  # Open DB connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  # Check if table exists
  tables <- DBI::dbListTables(ch3_db$con)
  
  if (!(table_name %in% tables)) {
    cat("Table '", table_name, "' does not exist in the database.\n", sep = "")
    ch3_db <- .ch3helper_closeDB(ch3_db)
    return(invisible(NULL))
  }
  
  # Query unique CpG (based on start + end positions)
  query <- paste0(
    "SELECT COUNT(*) AS unique_cpgs FROM (",
    "SELECT DISTINCT start, \"end\" FROM ", table_name,
    ")"
  )
  
  result <- DBI::dbGetQuery(ch3_db$con, query)
  
  # Close DB connection
  ch3_db <- .ch3helper_closeDB(ch3_db)
  
  cat("=================================================\n",
      "          Unique CpG Count in Calls Table        \n",
      "=================================================\n")
  cat("Unique CpG sites (by start+end):", result$unique_cpgs, "\n")
  
  invisible(result$unique_cpgs)
}