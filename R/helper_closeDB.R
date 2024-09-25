helper_closeDB <- function(ch3_db, db_con)
{
  # Finish up: update table list and close the connection
  ch3_db$tables <- dbListTables(db_con)
  dbDisconnect(db_con, shutdown = TRUE)
}