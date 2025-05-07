#' Summarize Methylation Data by Regions
#'
#' This function summarizes methylation data from a DuckDB database based on specified regions 
#' defined in a BED, TSV, or CSV file. It performs a join operation between the methylation data and 
#' the regions specified in the annotation file, allowing for different types of joins.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param table_name A string specifying what the user would like the name to be called in the database. Default is "regions".
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
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbExistsTable dbRemoveTable dbWriteTable
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl filter summarize mutate select all_of pull
#' @importFrom glue glue glue_collapse
#' @importFrom readr read_csv read_tsv
#' @importFrom tools file_ext
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  region_bed = system.file("Islands_hg38_test.csv", package = "MethylSeqR")
#'  
#'  # Summarize Regions using annotation table
#'  summarize_ch3_regions(ch3_db, region_bed)
#'
#' @export
summarize_ch3_regions <- function(ch3_db,
                              table_name = "regions",
                              region_file,
                              mod_type = c("c", "m", "h", "mh"),
                              chrs = c(as.character(1:22), 
                                      paste0("chr", 1:22), "chrX", "chrY", "chrM",
                                      paste0("Chr", 1:22), "ChrX", "ChrY", "ChrM"),
                              min_num_calls = 1)
{
  start_time <- Sys.time()
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
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(ch3_db$con, "PRAGMA max_temp_directory_size='100GiB';")
  dbExecute(ch3_db$con, "PRAGMA memory_limit='100GiB';")
  
  # Create regional data frame- offer a left, right or inner join
  # left join- keep reads that are outside of the annotation table
  # right join- keep reads in annotation table + regions not included in data frame
  # inner join- regions in both annotation and data frame...
  
  # Drop the regions table if it already exists
  if (dbExistsTable(ch3_db$con, table_name))
    dbRemoveTable(ch3_db$con, table_name)
  
  dbExecute(ch3_db$con, "VACUUM;")  # <-- Add this to free space immediately
  
  cat("Building regions table...")
  
  db_tbl = tbl(ch3_db$con, "calls") |>
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
  dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_annotation;")
  dbWriteTable(ch3_db$con, "temp_annotation", annotation, temporary = TRUE)
  
  # Upload positions (db_tbl) as a temporary table
  dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_positions;")
  dbWriteTable(ch3_db$con, "temp_positions", collect(db_tbl), temporary = TRUE)
  
  # Dynamically construct the SQL query for modification types
  count_columns <- ""
  frac_columns <- ""
  
  if ("c" %in% mod_type) {
    count_columns <- paste0("COALESCE(SUM(p.c_counts), 0) AS c_counts, ")
  }
  if ("m" %in% mod_type) {
    count_columns <- paste0(count_columns, "COALESCE(SUM(p.m_counts), 0) AS m_counts, ")
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
  CREATE TABLE ", table_name," AS
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
  
  dbExecute(ch3_db$con, query)
  
  # Drop temporary tables
  dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_annotation;")
  dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_positions;")
  
  cat("\n")
  end_time <- Sys.time()
  message("Windows table successfully created as ", table_name, " in database!\n", 
          "Time elapsed: ", end_time - start_time, "\n")
  
  print(head(tbl(ch3_db$con, table_name)))
  
  ch3_db$current_table = table_name
  ch3_db <- .ch3helper_cleanup(ch3_db)
  invisible(ch3_db)
}