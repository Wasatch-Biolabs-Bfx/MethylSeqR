#' Summarize Reads in a Database
#'
#' This function summarizes reads from a database, filtering and processing 
#' the data based on a provided key table (if given). It computes statistics 
#' on the reads such as the number of calls, CpG positions, and fractions of 
#' methylation (`m`), hemi-methylation (`h`), and total calls. The function 
#' interacts with the database to generate a `reads` table.
#'
#' @param ch3_db A character string specifying the path to the DuckDB database.
#' @param key_table (Optional) A character string specifying the path to a key table (CSV, TSV, or BED) used for filtering reads. If NULL, no filtering is applied.
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
#'  summarize_reads(ch3_db, region_bed)
#'
#' @export

summarize_reads <- function(ch3_db,
                            key_table = NULL,
                            min_CGs = 5) {
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do - purge extra tables, update table list, and close the connection
  on.exit(.helper_purgeTables(db_con), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(.helper_closeDB(database), add = TRUE)
  
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
  
  # If a key_table is provided, filter reads, THEN summarize.
  if (!is.null(key_table)) {
    .filter_then_sum_reads(db_con, 
                           key_table, 
                           min_CGs)
  } else {
    # If no key table, then just regularly summarize the reads.
    .sum_reads(db_con, 
               min_CGs)
  }
  
  cat("\n")
  message("reads table successfully created!")
  print(head(tbl(db_con, "reads")))
  
  ## TO DO: make classify_reads() function
  invisible(database)
}




.filter_then_sum_reads <- function(db_con, key_table, min_CGs) {
  # Read in key_table, determine the file type (csv, tsv, or bed)
  file_ext <- tools::file_ext(key_table)
  if (file_ext == "csv") {
    annotation <- read_csv(key_table, show_col_types = FALSE)
  } else if (file_ext %in% c("bed", "tsv")) {
    annotation <- read_tsv(key_table, show_col_types = FALSE)
  } else {
    stop("Invalid file type. Only CSV, TSV, or BED files are supported.")
  }
  
  # Ensure proper column names and remove header row if needed
  if (annotation[1,1] %in% c("chr", "Chr", "chrom", "Chrom")) {
    annotation = annotation[-1, ]
  }
  
  # Define the required columns
  required_columns <- c("chrom", "start", "end")
  existing_columns <- intersect(required_columns, colnames(annotation))
  
  # If any of the required columns are missing, stop with a message
  if(length(existing_columns) < length(required_columns)) {
    stop(cat("\nError: \"chrom\", \"start\", and \"end\" column required to filter reads.\nYou may need to rename headers if not in this format.\n"))
  }
  
  annotation <- annotation |>
    select(chrom, start, end)
  
  # Upload annotation as a temporary table
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_annotation;")
  dbWriteTable(db_con, "temp_annotation", annotation, temporary = TRUE)
  
  # Create the filtered_calls table
  dbExecute(db_con, "DROP TABLE IF EXISTS filtered_calls;")
  
  dbExecute(db_con, "
    CREATE TABLE filtered_calls AS 
    SELECT calls.*
    FROM calls
    JOIN temp_annotation
    ON calls.chrom = temp_annotation.chrom
    AND calls.start >= temp_annotation.start
    AND calls.end <= temp_annotation.end;
")
  
  # Summarize reads using the filtered calls table
  # Included: filter total_calls by the min_CGs parameter wanted!
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_annotation;")
  dbExecute(db_con, "DROP TABLE IF EXISTS reads;")
  
  query <- glue("
    CREATE TABLE reads AS 
    SELECT
        ANY_VALUE(sample_name) AS sample_name,
        read_id,
        ANY_VALUE(chrom) AS chrom,
        COUNT(*) AS total_calls,
        MIN(start) AS first_cpg_pos,
        MAX(\"end\") AS last_cpg_pos,
        ANY_VALUE(read_length) AS read_length,
        SUM(CASE WHEN call_code = '-' THEN 1 ELSE 0 END) AS c_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) AS m_counts,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) AS h_counts,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) AS mh_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS m_frac,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS h_frac,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS mh_frac
    FROM filtered_calls
    GROUP BY read_id
    HAVING total_calls >= {min_CGs};
")
  
  dbExecute(db_con, query)
  dbExecute(db_con, "DROP TABLE IF EXISTS filtered_calls;")
}




.sum_reads <- function(db_con, min_CGs) {
  # Summarize reads using the filtered calls table
  # Included: filter total_calls by the min_CGs parameter wanted!
  dbExecute(db_con, "DROP TABLE IF EXISTS reads;")
  
  query <- glue("
    CREATE TABLE reads AS 
    SELECT
        ANY_VALUE(sample_name) AS sample_name,
        read_id,
        ANY_VALUE(chrom) AS chrom,
        COUNT(*) AS total_calls,
        MIN(start) AS first_cpg_pos,
        MAX(\"end\") AS last_cpg_pos,
        ANY_VALUE(read_length) AS read_length,
        SUM(CASE WHEN call_code = '-' THEN 1 ELSE 0 END) AS c_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) AS m_counts,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) AS h_counts,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) AS mh_counts,
        SUM(CASE WHEN call_code = 'm' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS m_frac,
        SUM(CASE WHEN call_code = 'h' THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS h_frac,
        SUM(CASE WHEN call_code IN ('m', 'h') THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS mh_frac
    FROM calls
    GROUP BY read_id
    HAVING total_calls >= {min_CGs};
")
  
  dbExecute(db_con, query)
  
}