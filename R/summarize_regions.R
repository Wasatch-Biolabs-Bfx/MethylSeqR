#' Summarize Methylation Data by Regions
#'
#' This function summarizes methylation data from a DuckDB database based on specified regions 
#' defined in a BED, TSV, or CSV file. It performs a join operation between the methylation data and 
#' the regions specified in the annotation file, allowing for different types of joins.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param region_file A string representing the path to the BED or CSV file that contains the region annotations.
#' @param join_type A string indicating the type of join to perform. Options are "inner", "left", 
#' "right", or "full". Default is "inner".
#'
#' @details
#' The function reads the region annotations from the regional annotation file and checks for its validity.
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
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  region_bed = system.file("Islands_hg38_test.csv", package = "MethylSeqR")
#'  # Summarize Regions using annotation table
#'  summarize_regions(ch3_db, region_bed)
#'
#' @export
summarize_regions <- function(ch3_db,
                              region_file,
                              join_type = "inner",
                              min_cov = 1)
{
  # Determine the file type (csv, tsv, or bed)
  file_ext <- tools::file_ext(region_file)
  
  if (file_ext == "csv") {
    annotation <- read_csv(region_file,
                           col_names = c("chrom", "start", "end", "region_name"),
                           show_col_types = FALSE)
  } else if (file_ext %in% c("bed", "tsv")) {
    annotation <- read_tsv(region_file,
                           col_names = c("chrom", "start", "end", "region_name"),
                           show_col_types = FALSE)
  } else {
    stop("Invalid file type. Only CSV, TSV, or BED files are supported.")
  }
  
  # Ensure proper column names and remove header row if needed
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

  
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  on.exit(dbDisconnect(db_con, shutdown = TRUE), add = TRUE)
  on.exit(ch3_db$db_con <<- NULL, add = TRUE)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
  
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
  
  # Drop the regions table if it already exists
  dbExecute(db_con, "DROP TABLE IF EXISTS regions;")
  
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
  
  db_tbl = tbl(db_con, "calls") |>
    summarize(
      .by = c(sample_name, chrom, ref_position),
      cov = n(),
      c_counts = sum(as.integer(call_code == "-"), na.rm = TRUE),
      m_counts = sum(as.integer(call_code == "m"), na.rm = TRUE),
      h_counts = sum(as.integer(call_code == "h"), na.rm = TRUE),
      mh_counts = sum(as.integer(call_code %in% c("m", "h")), na.rm = TRUE)) |>
    mutate(
      m_frac = m_counts / cov,
      h_frac = h_counts / cov,
      mh_frac = mh_counts / cov) |> 
    filter(cov >= min_cov)
  
  for (chr in chroms) { # where error is happening
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
    
    # Drop the regions table if it already exists
    dbExecute(db_con, "DROP TABLE IF EXISTS temp_table;")
    
    pb$tick()
  }
  
  on.exit(dbRemoveTable(db_con, "regions"), add = TRUE)
  # Close progress bar
  pb$terminate()

  # Finish up: purge extra tables & update table list and close the connection
  keep_tables = c("positions", "regions", "windows", "meth_diff")
  .helper_purgeTables(db_con, keep_tables)
    
  # Finish Up
  ch3_db$tables <- dbListTables(db_con)
  database$last_table = "regions"
  .helper_closeDB(database)
  return(database)
}