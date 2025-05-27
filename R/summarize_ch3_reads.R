#' Summarize Reads in a Database
#'
#' This function summarizes reads from a database, filtering and processing 
#' the data based on a provided key table (if given). It computes statistics 
#' on the reads such as the number of calls, CpG positions, and fractions of 
#' methylation (`m`), hemi-methylation (`h`), and total calls. The function 
#' interacts with the database to generate a `reads` table.
#'
#' @param ch3_db A character string specifying the path to the DuckDB database.
#' @param table_name A string specifying what the user would like the name to be called in the database. Default is "reads".
#' @param min_length An integer specifying the the minimum read_length.
#' @param min_CGs An integer specifying the minimum number of CG sites required for a read to be included in the summary.
#'
#' @return Invisibly returns the database object. The function also outputs a success message and the first few rows of the summarized `reads` table.
#'
#' @details
#' The function connects to the provided DuckDB database, optionally filters reads based on the key table, and then summarizes the read data. It creates a temporary table for the filtered reads (if a key table is provided) and creates a summary table called `reads` with information on the total number of calls, the positions of the first and last CG sites, and counts for different types of calls (`m`, `h`, and `-`).
#'
#' @examples
#' #Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  region_bed = system.file("Islands_hg38_test.csv", package = "MethylSeqR")
#'  
#'  # Summarize Reads
#'  summarize_ch3_reads(ch3_db, region_bed)
#'
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbExistsTable dbRemoveTable dbWriteTable
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl select
#' @importFrom glue glue
#' @importFrom tools file_ext
#' @importFrom readr read_csv read_tsv
#' @export

# Bucket Approach
summarize_ch3_reads <- function(ch3_db,
                                input_calls_table = "calls",
                                output_reads_table = "reads",
                                regions_table = NULL,
                                min_length = 100,
                                min_CGs = 5) 
{
  start_time <- Sys.time()
  
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  dbExecute(ch3_db$con, "PRAGMA max_temp_directory_size='100GiB';")
  
  # Optionally read and upload regions table if provided
  if (!is.null(regions_table)) {
    file_ext <- tools::file_ext(regions_table)
    if (file_ext == "csv") {
      annotation <- readr::read_csv(regions_table, show_col_types = FALSE)
    } else if (file_ext %in% c("bed", "tsv")) {
      annotation <- readr::read_tsv(regions_table, show_col_types = FALSE)
    } else {
      stop("Invalid regions table file type. Use CSV, TSV, or BED.")
    }
    
    required_cols <- c("chrom", "start", "end")
    if (!all(required_cols %in% colnames(annotation))) {
      stop("Regions table must include columns: chrom, start, end")
    }
    
    dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_regions_table;")
    annotation <- annotation |> dplyr::mutate(bucket = floor(start / 10000))
    DBI::dbWriteTable(ch3_db$con, "temp_regions_table", annotation, temporary = TRUE)
  }
  
  cat("Building reads table...")
  if (dbExistsTable(ch3_db$con, output_reads_table))
    dbRemoveTable(ch3_db$con, output_reads_table)
  
  query <- glue::glue("
    CREATE TABLE {output_reads_table} AS
    SELECT
        ANY_VALUE(c.sample_name) AS sample_name,
        c.read_id,
        ANY_VALUE(c.chrom) AS chrom,
        COUNT(*) AS total_calls,
        MIN(c.start) AS first_cpg_pos,
        MAX(c.end) AS last_cpg_pos,
        ANY_VALUE(c.read_length) AS read_length,
        SUM(CASE WHEN c.call_code = '-' THEN 1 ELSE 0 END) AS c_counts,
        SUM(CASE WHEN c.call_code = 'm' THEN 1 ELSE 0 END) AS m_counts,
        SUM(CASE WHEN c.call_code = 'h' THEN 1 ELSE 0 END) AS h_counts,
        SUM(CASE WHEN c.call_code IN ('m', 'h') THEN 1 ELSE 0 END) AS mh_counts,
        SUM(CASE WHEN c.call_code = 'm' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS m_frac,
        SUM(CASE WHEN c.call_code = 'h' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS h_frac,
        SUM(CASE WHEN c.call_code IN ('m', 'h') THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS mh_frac
    FROM {input_calls_table} c
    {if (!is.null(regions_table))
        \"JOIN temp_regions_table k
        ON c.chrom = k.chrom
        AND k.bucket = FLOOR(c.start / 10000)
        AND c.start >= k.start
        AND c.start <= k.end\" else \"\"}
    GROUP BY c.read_id
    HAVING total_calls >= {min_CGs}
       AND read_length >= {min_length}")
  
  dbExecute(ch3_db$con, query)
  
  cat("\n")
  end_time <- Sys.time()
  total_seconds <- as.numeric(end_time - start_time, units = "secs")
  
  if (total_seconds > 60) {
    message("Reads table successfully created as ", output_reads_table, " in database!",
            "\nTime elapsed: ", round(total_seconds / 60, 2), " minutes\n")
  } else {
    message("Reads table successfully created as ", output_reads_table, " in database!", 
            "\nTime elapsed: ", round(total_seconds, 2), " seconds\n")
  }
  
  print(head(dplyr::tbl(ch3_db$con, output_reads_table)))
  ch3_db$current_table = output_reads_table
  ch3_db <- .ch3helper_cleanup(ch3_db)
  invisible(ch3_db)
}