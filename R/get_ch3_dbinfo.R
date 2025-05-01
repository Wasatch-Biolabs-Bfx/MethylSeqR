#' Get Database Statistics
#'
#' Prints out a summary of the ch3 database, including size, tables, and unique sample names.
#'
#' @param ch3_db Path to the `.ch3.db` file or a `ch3_db` object.
#' 
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Get database statistics
#'  get_ch3_dbinfo(ch3_db = ch3_db)
#'
#' @return Invisibly returns a list of stats from the database.
#' 
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbListTables
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl distinct arrange collect pull
#' 
#' @export

get_ch3_dbinfo <- function(ch3_db) 
{
  
  if (is.character(ch3_db)) {
    # Check if the file exists
    if (!file.exists(ch3_db)) {
      stop(paste("The database file", ch3_db, "does not exist."))
    }
    
    cat("=================================================\n",
        "               Database Statistics               \n",
        "=================================================\n", 
        "Database: ", ch3_db, "\n",
        sep = "")
  } else {
    cat("=================================================\n",
        "               Database Statistics               \n",
        "=================================================\n", 
        sep = "")
  }
  
  # Open the database connection
  database <- .ch3helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...close the connection
  on.exit(.ch3helper_closeDB(database), add = TRUE)
  
  # Get DB size
  size_df <- dbGetQuery(db_con, "PRAGMA database_size")
  size <- sum(size_df$total_blocks * size_df$block_size) / 1024 / 1024
  cat(sprintf("\nDatabase Size: %.2f MB\n", size))
  
  # What tables are in the database?
  tables <- dbListTables(db_con)
  cat("\nTables in Database:\n")
  cat(paste(tables, collapse = "\n"), "\n")
  
  # Unique Sample Names (if "calls" table exists)
  if ("calls" %in% tables) {
    sample_names <- tbl(db_con, "calls") %>%
      distinct(sample_name) %>%
      arrange(sample_name) %>%
      collect() %>%
      pull(sample_name)
    
    cat("\nUnique Sample Names:\n")
    cat(paste(sample_names, collapse = "\n"), "\n")
  } else {
    cat("\nNo 'calls' table found in the database.\n")
    sample_names <- character(0)
  }
  
  # Assemble output list
  stats <- list(
    tables = tables,
    sample_names = sample_names,
    num_samples = length(sample_names)
  )
  
  return(invisible(stats))
  
  }