export_tables <- function(ch3_db,
                         tables = c("positions"),
                        out_path) 
{
  # Open the database connection
  db_con <- helper_connectDB(ch3_db)
  
  tryCatch(
    {
      all_tables <- dbListTables(db_con)
      
      for (tbl in tables) {
        if (tbl %in% all_tables) {
          # Print the table name
          message(paste0("Writing out ", tbl, " table..."))
          # write out table to path given
          write.csv(tbl(db_con, tbl), file = out_path, row.names = FALSE)
          
          
        } else {
          # Print a message if the table does not exist
          message(paste0("Table '", tbl, "' does not exist in the database."))
        }
      }
    }, 
    error = function(e) 
    {
      # Print custom error message
      message("An error occurred: ", e$message)
    }, 
    finally = 
      {
        helper_closeDB(ch3_db, db_con)
      })
}