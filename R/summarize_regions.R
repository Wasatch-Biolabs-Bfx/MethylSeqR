#' Summarize Methylation Data by Regions
#'
#' This function summarizes methylation data from a DuckDB database based on specified regions 
#' defined in a BED file. It performs a join operation between the methylation data and the regions 
#' specified in the annotation file, allowing for different types of joins.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param region_bed A string representing the path to the BED file that contains the region annotations.
#' @param join_type A string indicating the type of join to perform. Options are "inner", "left", 
#' "right", or "full". Default is "inner".
#'
#' @details
#' The function reads the region annotations from the specified BED file and checks for its validity.
#' It connects to the DuckDB database, creates a summarized table of methylation data based on the specified 
#' regions, and performs the join operation according to the specified join type. A progress bar is displayed 
#' during the summarization process. The resulting data is stored in a table called "regions" within the database.
#'
#' @return The updated `ch3_db` object with the summarized regions data added to the DuckDB database.
#'
#' @import readr
#' @import dplyr
#' @import dbplyr
#' @import duckdb
#' @import duckplyr
#' @import progress
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylseqR")
#'  ch3_db <- file.path(ch3_db, "my_data.ch3.db")  # Path to the pre-existing database
#'  
#'  region_bed = system.file("Islands_hg38_ucsc.csv", package = "MethylseqR")
#'  # Summarize Regions using anotation table
#'  summarize_regions(ch3_db, region_bed)
#'
#' @export
summarize_regions <- function(ch3_db,
                              region_bed,
                              join_type = "inner")
{
  # Read annotation
  annotation <-
    read_csv(region_bed,
             col_names = c("chrom", "start", 
                           "end", 
                           "region_name"),
             show_col_types = FALSE)
  
  # make sure column names did not get included and mess up code...
  if (annotation[1,1] %in% c("chr", "Chr", "chrom", "Chrom")) {
    annotation = annotation[-1, ]
  }
  
  # check format
  if (ncol(annotation) < 3 || ncol(annotation) > 4) {
    stop("Invalid annotation format. File must have 3 or 4 columns:
         chr, start, end, region_name (optional) annotation.")
  }
  
  if (ncol(annotation == 3)) {
    annotation <-
      annotation |>
      mutate(
        region_name = paste(chrom, start, end, sep = "_"))
  }
  
  annotation <-
    annotation |>
    reframe(
      .by = c(region_name, chrom),
      ref_position = start:end)

  tryCatch(
    {
      # Open the database connection
      database <- .helper_connectDB(ch3_db)
      db_con <- database$db_con
      # Create regional data frame- offer a left, right or inner join
      
      # left join- keep reads that are outside of the annotation table
      # right join- keep reads in annotation table + regions not included in data frame
      # inner join- regions in both annotation and data frame...
      if (!join_type %in% c("inner", "right", "left", "full")) {
        stop("Invalid join type")
      }
      
      my_join <- switch(join_type,
                        "inner" = inner_join,
                        "right" = right_join,
                        "left" = left_join,
                        "full" = full_join)
      
      if (dbExistsTable(db_con, "regions"))
        dbRemoveTable(db_con, "regions")
      
      # Create regions table by unique chroms
      chroms <-
        annotation |>
        select(chrom) |>
        distinct() |>
        filter(nchar(chrom) < 6) |>
        pull()
      
      cat("Building regions table...")
      # Create Progress Bar
      pb <- progress_bar$new(
        format = "[:bar] :percent [Elapsed time: :elapsed]",
        total = length(chroms) + 1,
        complete = "=",
        incomplete = "-",
        current = ">",
        clear = FALSE,
        width = 100)
      
      pb$tick()
      
      db_tbl = tbl(db_con, "positions")
      
      for (chr in chroms) {
        # Begin summarizing by region- perform the join and aggregation
        db_tbl |>
          filter(chrom == chr) |>
          my_join(annotation, 
                  by = join_by(chrom, 
                               ref_position), 
                  copy = TRUE) |>
          summarize(
            .by = c(sample_name, region_name),
            cov = sum(cov, na.rm = TRUE),
            across(ends_with("_counts"), ~ sum(.x, na.rm = TRUE)),
            across(ends_with("_frac"), ~ sum(.x * cov, na.rm = TRUE) / sum(cov, na.rm = TRUE))) |>
          compute(name = "temp_table", temporary = TRUE)
        
        # Create or append table
        dbExecute(db_con, 
                  "CREATE TABLE IF NOT EXISTS regions AS 
              SELECT * FROM temp_table WHERE 1=0")
        
        dbExecute(db_con, 
                  "INSERT INTO regions SELECT * FROM temp_table")
        
        dbRemoveTable(db_con, "temp_table")
        
        pb$tick()
      }
      
      # Close progress bar
      pb$terminate()
    }, error = function(e)
    {
      # Print custom error message
      message("An error occurred: ", e$message)
      dbRemoveTable(db_con, "regions")
    }, 
    finally = 
      {
        # Finish up: purge extra tables & update table list and close the connection
        keep_tables = c("positions", "regions", "windows", "meth_diff")
        .helper_purgeTables(db_con, keep_tables)
        
        # Finish Up
        database$last_table = "regions"
        .helper_closeDB(database)
        return(database)
      }
  )
}