.ch3helper_cleanup <- function(ch3_db)
{
  # purge extra tables, update table list, and then close the connection
  ch3_db$con <- .ch3helper_purgeTables(ch3_db$con)  # Purge tables FIRST
  dbExecute(ch3_db$con, "VACUUM;")  # <-- Ensure space is reclaimed\
  ch3_db <- .ch3helper_closeDB(ch3_db)     # Close DB LAST 
  
  invisible(ch3_db)
}