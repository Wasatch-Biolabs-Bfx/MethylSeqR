helper_print <- function(ch3_db, tables = c("positions"))
{
  # Open the database connection
  db_con <- helper_connectDB(ch3_db)
  
  # List all tables in the database
  all_tables <- dbListTables(db_con)
  message(paste0("All tables currently in database: "))
  
  for (tbl_name in all_tables) {
    print(tbl_name)
  }
  
  # if (tables == "all") {
  #   tables = all_tables
  # }
  
  for (tb_name in tables) {
    if (tb_name %in% all_tables) {
      # Print the table name
      message(paste0("Table: ", tb_name))
      
      # Fetch the first few rows and print them
      print(tbl(db_con, tb_name) %>% head())

    } else {
      # Print a message if the table does not exist
      message(paste0("Table '", tb_name, "' does not exist in the database."))
    }
  }
  
  helper_closeDB(ch3_db, db_con)
}