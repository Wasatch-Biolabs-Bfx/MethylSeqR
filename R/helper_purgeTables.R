helper_purgeTables <- function(ch3_db, keep_tables = c("positions"))
{
  # Open the database connection
  db_con <- helper_connectDB(ch3_db)
  
  # List all tables in the database
  all_tables <- dbListTables(db_con)
  
  # Remove tables that are not in the 'keep_tables' list
  for (table in all_tables) {
    if (!(table %in% keep_tables)) {
      message(paste0("Removing table: ", table))
      dbRemoveTable(db_con, table)
    }
  }
  
   print("Remaining tables in database: ")
   print(dbListTables(db_con))
   
   helper_closeDB(ch3_db, db_con)
}