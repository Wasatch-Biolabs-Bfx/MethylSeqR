helper_purgeTables <- function(db_con, keep_tables = c("positions"))
{
  # List all tables in the database
  all_tables <- dbListTables(db_con)
  
  # Remove tables that are not in the 'keep_tables' list
  for (table in all_tables) {
    if (!(table %in% keep_tables)) {
      message(paste0("Removing table: ", table))
      dbRemoveTable(db_con, table)
    }
  }
   print(dbListTables(db_con))
   
   return(db_con)
}