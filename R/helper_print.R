helper_print <- function(ch3_db, tables = "all") {
  # Open the database connection
  db_con <- helper_connectDB(ch3_db)
  
  tryCatch(
    {
      # List all tables in the database
      all_tables <- dbListTables(db_con)
      cat("All tables currently in database:\n")
      
      # Print each table name
      for (tbl_name in all_tables) {
        cat(tbl_name, "\n")
      }
      
      if (tables == "all") {
        tables <- all_tables
      }
      
      for (tb_name in tables) {
        if (tb_name %in% all_tables) {
          # Print the table name
          cat("Table: ", tb_name, "\n")
          
          # Fetch the first few rows and print them
          print(head(tbl(db_con, tb_name)))
          
        } else {
          # Print a message if the table does not exist
          message(paste0("Table '", tb_name, "' does not exist in the database."))
        }
      }
    }, 
    error = function(e) {
      # Print custom error message
      message("An error occurred: ", e$message)
    }, 
    finally = {
      helper_closeDB(ch3_db, db_con)
    }
  )
}
