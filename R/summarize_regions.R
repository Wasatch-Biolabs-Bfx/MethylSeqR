#' Summarize Methylation Data by Regions
#'
#' This function summarizes methylation data from a DuckDB database based on specified regions 
#' defined in a BED, TSV, or CSV file. It performs a join operation between the methylation data and 
#' the regions specified in the annotation file, allowing for different types of joins.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param region_file A string representing the path to the BED or CSV file that contains the region annotations.
#' @param mod_type A character vector specifying the modification types to include. Options are  `"c"` (unmodified cytosine),
#' `"m"` (methylation), `"h"` (hydroxymethylation), 
#'   and `"mh"` (methylated + hydroxymethylated).
#' @param chrs A character vector specifying which chromosomes to include. Default includes all autosomes, 
#'   sex chromosomes (chrX, chrY), and mitochondrial chromosome (chrM).
#' @param min_num_calls An integer specifying the minimum number of calls required for inclusion in the summary.
#'   Default is 1.
#'
#' @details
#' The function reads the region annotations from the regional annotation file and checks for its validity.
#' It connects to the DuckDB database, creates a summarized table of methylation data based on the specified 
#' regions, and performs the join operation according to the specified join type. A progress bar is displayed 
#' during the summarization process. The resulting data is stored in a table called "regions" within the database.
#'
#' @return The updated `ch3_db` object with the summarized regions data added to the DuckDB database.
#'
#' @importFrom tools file_ext
#' @importFrom readr read_csv read_tsv
#' @importFrom dplyr mutate summarize filter pull tbl
#' @importFrom DBI dbExecute dbWriteTable
#' @importFrom duckdb duckdb
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  region_bed = system.file("Islands_hg38_test.csv", package = "MethylSeqR")
#'  
#'  # Summarize Regions using annotation table
#'  summarize_regions(ch3_db, region_bed)
#'
#' @export
summarize_regions <- function(ch3_db,
                              region_file,
                              mod_type = c("c", "m", "h", "mh"),
                              chrs = c(as.character(1:22), 
                                      paste0("chr", 1:22), "chrX", "chrY", "chrM",
                                      paste0("Chr", 1:22), "ChrX", "ChrY", "ChrM"),
                              min_num_calls = 1)
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
  
  # check format
  if (ncol(annotation) < 3 || ncol(annotation) > 4) {
    stop("Invalid annotation format. File must have 3 or 4 columns:
         chr, start, end, region_name (optional) annotation.")
  }
  
  # Ensure proper column names and remove header row if needed
  if (annotation[1,1] %in% c("chr", "Chr", "chrom", "Chrom")) {
    annotation = annotation[-1, ]
  }
  
  if (ncol(annotation == 3)) {
    annotation <-
      annotation |>
      mutate(
        region_name = paste(chrom, start, end, sep = "_"))
  }

  
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # Finish up: purge extra tables & update table list and close the connection
  on.exit(.helper_purgeTables(db_con), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(.helper_closeDB(database), add = TRUE)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
  dbExecute(db_con, "PRAGMA memory_limit='100GiB';")
  
  
  # Create regional data frame- offer a left, right or inner join
  # left join- keep reads that are outside of the annotation table
  # right join- keep reads in annotation table + regions not included in data frame
  # inner join- regions in both annotation and data frame...
  
  # Drop the regions table if it already exists
  dbExecute(db_con, "DROP TABLE IF EXISTS regions;")
  dbExecute(db_con, "VACUUM;")  # <-- Add this to free space immediately
  
  cat("Building regions table...")
  
  db_tbl = tbl(db_con, "calls") |>
    filter(chrom %in% chrs) |>  # Filter for selected chromosomes
    summarize(
      .by = c(sample_name, chrom, start, end),
      num_calls = n(),
      c_counts = sum(as.integer(call_code == "-"), na.rm = TRUE),
      m_counts = sum(as.integer(call_code == "m"), na.rm = TRUE),
      h_counts = sum(as.integer(call_code == "h"), na.rm = TRUE),
      mh_counts = sum(as.integer(call_code %in% c("m", "h")), na.rm = TRUE)) |>
    mutate(
      m_frac = ifelse(num_calls == 0, NA_real_, m_counts / num_calls),
      h_frac = ifelse(num_calls == 0, NA_real_, h_counts / num_calls),
      mh_frac = ifelse(num_calls == 0, NA_real_, mh_counts / num_calls)) |> 
    filter(num_calls >= min_num_calls)
  
  has_missing_chr <- db_tbl |> 
    summarize(any_missing_chr = any(!grepl("^chr", chrom))) |> 
    pull(any_missing_chr)
  
  if (has_missing_chr) {
    has_chr <- any(grepl("^chr", annotation$chrom))  # TRUE if any row has "chr"
    
    if (has_chr) {
      annotation <- annotation |> 
        mutate(chrom = sub("^chr", "", chrom))  # Base R equivalent
    }
  }
  
  # Select only requested modtype columns (always keeping num_calls)
  selected_columns <- c("sample_name", "chrom", "start", "end", "num_calls", 
                        paste0(mod_type, "_counts"), paste0(mod_type, "_frac"))
  selected_columns <- intersect(selected_columns, colnames(db_tbl))
  
  db_tbl <- db_tbl |> select(all_of(selected_columns))
  
  # Upload annotation as a temporary table
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_annotation;")
  dbWriteTable(db_con, "temp_annotation", annotation, temporary = TRUE)
  
  # Upload positions (db_tbl) as a temporary table
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_positions;")
  dbWriteTable(db_con, "temp_positions", collect(db_tbl), temporary = TRUE)
  
  # Dynamically construct the SQL query for modification types
  count_columns <- ""
  frac_columns <- ""
  
  if ("c" %in% mod_type) {
    count_columns <- paste0("COALESCE(SUM(p.c_counts), 0) AS c_counts, ")
  }
  if ("m" %in% mod_type) {
    count_columns <- paste0("COALESCE(SUM(p.m_counts), 0) AS m_counts, ")
    frac_columns <- paste0(frac_columns, "COALESCE(SUM(p.m_counts * p.num_calls) / NULLIF(SUM(p.num_calls), 0), 0) AS m_frac, ")
  }
  if ("h" %in% mod_type) {
    count_columns <- paste0(count_columns, "COALESCE(SUM(p.h_counts), 0) AS h_counts, ")
    frac_columns <- paste0(frac_columns, "COALESCE(SUM(p.h_counts * p.num_calls) / NULLIF(SUM(p.num_calls), 0), 0) AS h_frac, ")
  }
  if ("mh" %in% mod_type) {
    count_columns <- paste0(count_columns, "COALESCE(SUM(p.mh_counts), 0) AS mh_counts, ")
    frac_columns <- paste0(frac_columns, "COALESCE(SUM(p.mh_counts * p.num_calls) / NULLIF(SUM(p.num_calls), 0), 0) AS mh_frac, ")
  }
  
  # Remove trailing commas
  count_columns <- substr(count_columns, 1, nchar(count_columns) - 2)
  frac_columns <- substr(frac_columns, 1, nchar(frac_columns) - 2)
  
  # Create the final query with dynamic columns
  query <- paste0("
  CREATE TABLE regions AS
  SELECT 
    p.sample_name, 
    a.region_name,
    a.chrom, 
    a.start, 
    a.end, 
    COUNT(*) AS num_CpGs, 
    NULLIF(SUM(p.num_calls), 0) AS num_calls, 
    ", count_columns, ", 
    ", frac_columns, "
  FROM temp_positions p
  JOIN temp_annotation a
    ON p.chrom = a.chrom 
    AND CAST(p.start AS DOUBLE) BETWEEN CAST(a.start AS DOUBLE) AND CAST(a.end AS DOUBLE)
  GROUP BY p.sample_name, a.region_name, a.chrom, a.start, a.end;
  ")
  
  dbExecute(db_con, query)
  
  # Drop temporary tables
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_annotation;")
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_positions;")
  
  cat("\n")
  message("Regions table successfully created!")
  print(head(tbl(db_con, "regions")))
    
  invisible(database)
}