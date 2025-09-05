#' Execute an Expression on a DuckDB Table with Optional Materialization
#'
#' Connects to a DuckDB database, evaluates a user-supplied expression on a specified table,
#' and either collects the result into R or computes and stores it as a new table in the database.
#'
#' @param ch3_db Path to the DuckDB database file (e.g., `"my_data.ch3.db"`).
#' @param table_name Name of the table in the database to operate on.
#' @param expr A function taking one argument (`tbl_ref`) and returning a lazy `dplyr` expression.
#'   The table reference (`tbl_ref`) will be passed as a `tbl()` object connected to the database.
#' @param mode One of `"collect"` (default) or `"compute"`. If `"collect"`, the result is returned as
#'   a data frame in R. If `"compute"`, the result is stored as a new table in the database.
#' @param output_table Required if `mode = "compute"`. Name of the output table to create or overwrite
#'   with the result of `expr(tbl_ref)`.
#'
#' @return If `mode = "collect"`, returns a data frame. If `mode = "compute"`, returns `NULL` invisibly.
#'
#' @examples
#' \dontrun{
#' # Collect results of a filtered table into R
#' run_ch3_dplyr(
#'   ch3_db = "my_data.ch3.db",
#'   table_name = "methylation_data",
#'   expr = function(tbl_ref) dplyr::filter(tbl_ref, score > 0.5),
#'   mode = "collect"
#' )
#'
#' # Store the filtered result in a new table inside the database
#' run_ch3_dplyr(
#'   ch3_db = "my_data.ch3.db",
#'   table_name = "methylation_data",
#'   expr = function(tbl_ref) dplyr::filter(tbl_ref, score > 0.5),
#'   mode = "compute",
#'   output_table = "filtered_data"
#' )
#' }
#' 
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbListTables
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl collect compute
#' 
#' @export

run_ch3_dplyr <- function(
    ch3_db, 
    table_name, 
    expr, 
    mode = c("collect", "compute"), 
    output_table = NULL) 
{

  mode <- match.arg(mode)
  
  start_time <- Sys.time()
  # Connect to the database
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  tbl_ref <- tbl(ch3_db$con, table_name)
  
  result <- force(expr(tbl_ref))
  
  if (mode == "collect") {
    # Just collect the results into R
    out <- collect(result)
    
    end_time <- Sys.time()
    message("Query Finished. Time elapsed: ", end_time - start_time, "\n")
    ch3_db <- .ch3helper_closeDB(ch3_db)
    
    return(out)
  } else if (mode == "compute") {
    if (is.null(output_table)) {
      stop("output_table must be provided when mode = 'compute'.", call. = FALSE)
    }
    
    # Drop existing table
    dbExecute(ch3_db$con, paste0("DROP TABLE IF EXISTS ", output_table))
    
    # Compute into a new table
    computed <- compute(result, name = output_table, temporary = FALSE)
    
    end_time <- Sys.time()
    message("Query Finished. Time elapsed: ", end_time - start_time, "\n")
    ch3_db <- .ch3helper_closeDB(ch3_db)
    
    invisible(ch3_db)
  }
}

#' Execute a query on a CH3 Database
#'
#' Connects to a DuckDB database, evaluates a user-supplied sql query and closes the database.
#'
#' @param ch3_db Path to the DuckDB database file (e.g., `"my_data.ch3.db"`).
#' @param query An sql query supported by the duckdb database framework.
#'
#' @return a ch3_db object to allow piping
#'
#' @examples
#' \dontrun{
#' # Count the number of rows in the calls table
#' run_ch3_sql(
#'   ch3_db = "my_data.ch3.db",
#'   query = "CREATE TABLE call_count AS SELECT COUNT(*) AS num_rows FROM calls;")
#' }
#' 
#' @importFrom DBI dbExecute 
#' 
#' @export

run_ch3_sql <- function(ch3_db, 
                        query) 
{
  start_time <- Sys.time()
  # Connect to the database
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  dbExecute(ch3_db$con, query)
  
  end_time <- Sys.time()
  message("Query Finished. Time elapsed: ", end_time - start_time, "\n")
  ch3_db <- .ch3helper_closeDB(ch3_db)
  
  invisible(ch3_db)
}
