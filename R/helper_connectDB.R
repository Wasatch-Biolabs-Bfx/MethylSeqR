helper_connectDB <- function(ch3_db)
{
  # Open the database connection
  if (inherits(ch3_db, "character")) { # if given a string file name
    ch3_db <- list(db_file = ch3_db)
    class(ch3_db) <- "ch3_db"
    db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
    ch3_db$tables <- dbListTables(db_con)
  } else if (inherits(ch3_db, "ch3_db")) { # if given a 
    db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
  } else {
    stop("Invalid ch3_db class. Must be character or ch3_db.")
  }
  
  return(db_con)
}