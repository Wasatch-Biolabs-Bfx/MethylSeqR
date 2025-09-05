#' Get Table Information from ch3 Database
#'
#' Provides a summary of a specified table within a `.ch3.db` database, including total number of records,
#' sample-level record counts (if applicable), and column names.
#'
#' @param ch3_db Path to the `.ch3.db` file or a `ch3_db` object.
#' @param table_name Name of the table to summarize.
#'
#' @examples
#' # Specify the path to the database
#' ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'
#' # Get information about the 'calls' table
#' get_ch3_tableinfo(ch3_db = ch3_db, table_name = "calls")
#'
#' @return Invisibly returns the closed `ch3_db` object after summarizing the specified table.
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbListTables dbListFields
#'
#' @export

get_ch3_tableinfo <- function(ch3_db, 
                              table_name) 
{
  if (is.character(ch3_db)) if (!file.exists(ch3_db)) {
    stop(paste("The database file", ch3_db, " does not exist."))
  }
  
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  tables <- DBI::dbListTables(ch3_db$con)
  
  # If the table does not exist
  if (!(table_name %in% tables)) {
    cat("Table '", table_name, "' does not exist in the database.\n", sep = "")
    ch3_db <- .ch3helper_closeDB(ch3_db)
    return(invisible(ch3_db))
  }
  
  # Else, print db stats
  cat("=================================================\n",
      "                 Table Statistics                \n",
      "=================================================\n",
      sep = "")
  
  cat("Table:", table_name, "\n")
  
  # Total number of records
  num_records <- DBI::dbGetQuery(ch3_db$con, paste0("SELECT COUNT(*) AS n FROM ", table_name))$n
  
  cat("\nTotal Number of Records:", num_records, "\n")
  
  # Records by sample_name (if applicable)
  col_names <- DBI::dbListFields(ch3_db$con, table_name)
  
  if ("sample_name" %in% col_names) {
    sample_counts <- DBI::dbGetQuery(
      ch3_db$con,
      paste0("SELECT sample_name, COUNT(*) AS n FROM ", table_name,
             " GROUP BY sample_name ORDER BY sample_name")
    )
    
    apply(sample_counts, 1, function(row) {
      cat(sprintf("  %s: %d\n", row[["sample_name"]], as.integer(row[["n"]])))
    })
  }
  
  # Indented column names
  cat("\nColumns:\n")
  for (col in col_names) {
    cat("  ", col, "\n")
  }
  
  ch3_db <- .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}