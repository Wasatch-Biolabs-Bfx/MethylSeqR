#' Summarize methylation in sliding windows
#'
#' Aggregates call-level rows from \code{input_table} (typically \code{"calls"})
#' into per-sample window summaries written to \code{output_table} (default \code{"windows"}).
#' Windows are created using a \strong{tiling + offsets} strategy so that every
#' \code{start} position is assigned to one or more windows depending on the
#' \code{step_size}. For each window the function computes:
#' total CpG sites, total calls, per-class counts, and per-class fractions.
#'
#' @section How modification codes work:
#' Pass \code{mod_code} as single codes (e.g. \code{"m"}, \code{"h"}, \code{"a"})
#' or combinations with \code{"+"} (e.g. \code{"m + h"}). Labels are created by
#' removing spaces and \code{"+"} (e.g. \code{"m + h"} \u2192 \code{"mh"}). For each
#' label the table includes \code{<label>_counts} and \code{<label>_frac}. The
#' unmodified class is defined by \code{unmod_code} (default \code{"-"}), named
#' using \code{unmod_label} (default \code{"c"} \u2192 \code{c_counts}, \code{c_frac}).
#'
#' @param ch3_db Path to a \code{.ch3.db} DuckDB file or a \code{"ch3_db"} object.
#'   A connection is opened via internal helpers and closed/cleaned on return.
#' @param input_table Source table containing call-level records (default \code{"calls"}).
#'   Must contain at least: \code{sample_name}, \code{chrom}, \code{start}, \code{call_code}.
#' @param output_table Destination window table name (default \code{"windows"}).
#' @param window_size Integer window width in bases (default \code{1000}).
#' @param step_size Step, in bases, used to create staggered window offsets
#'   (default \code{10}). Offsets are \code{seq(1, window_size - 1, by = step_size)}.
#' @param chrs Character vector of chromosome filters; rows whose \code{chrom} match any
#'   value are retained. Defaults to common human aliases (1–22, \code{chrX}, \code{chrY}, \code{chrM}, …).
#' @param samples Optional character vector of \code{sample_name}s to include. If \code{NULL},
#'   all samples present in \code{input_table} are processed.
#' @param mod_code Character vector of modification specs to count (single codes or
#'   \code{"code1 + code2"} combinations). Default \code{c("m","h","m + h")}.
#' @param unmod_code Call code representing unmodified (default \code{"-"}).
#' @param unmod_label Label used to name unmodified columns (default \code{"c"}).
#' @param min_num_calls Minimum total calls required for a window to be written
#'   (default \code{1}). Windows below this threshold are skipped.
#' @param temp_dir Directory for DuckDB temporary files (default \code{tempdir()}).
#' @param threads Integer DuckDB thread count. If \code{NULL}, an internal heuristic
#'   (typically all-but-one core) is used.
#' @param memory_limit DuckDB memory limit string (e.g. \code{"16384MB"}).
#'   If \code{NULL}, an internal heuristic (~80\% of RAM) is used.
#' @param overwrite If \code{TRUE} and \code{output_table} exists, it is dropped before writing.
#'
#' @details
#' For each sample, the function first aggregates per-position counts from
#' \code{input_table} (\code{num_calls} plus dynamically generated
#' \code{<label>_counts} per \code{mod_code}/\code{unmod_code}). It then creates
#' sliding windows by assigning each position to a window start computed as:
#' \deqn{temp\_start = start - ((start - offset) \bmod window\_size).}
#' For each \code{offset} in \code{seq(1, window_size - 1, by = step_size)}, it
#' sums counts over \code{[temp_start, temp_start + window_size - 1]} and writes:
#' \itemize{
#'   \item \code{num_CpGs}: number of positions aggregated in the window
#'   \item \code{num_calls}: sum of \code{num_calls}
#'   \item \code{<label>_counts}: summed counts for each label
#'   \item \code{<label>_frac}: \code{<label>_counts / num_calls} (NULL if \code{num_calls == 0})
#' }
#' Resource pragmas (\code{temp_directory}, \code{threads}, \code{memory_limit}) are set
#' via internal heuristics unless overridden.
#'
#' @return (Invisibly) a \code{"ch3_db"} object pointing to the same DB file with
#'   \code{current_table} set to \code{output_table}. The created table has columns:
#'   \itemize{
#'     \item \code{sample_name}, \code{chrom}, \code{start}, \code{end},
#'           \code{num_CpGs}, \code{num_calls},
#'     \item for each label in \code{c(unmod_label, parsed(mod_code))}:
#'           \code{<label>_counts}, \code{<label>_frac}.
#'   }
#'
#' @examples
#' \dontrun{
#' # Default m/h windows (1kb windows, 10bp staggered offsets)
#' summarize_ch3_windows("my_db.ch3.db")
#'
#' # Custom mod codes with a novel 'a' code and stricter filtering
#' summarize_ch3_windows(
#'   ch3_db        = "my_db.ch3.db",
#'   mod_code      = c("a", "m + h"),
#'   min_num_calls = 25,
#'   window_size   = 2000,
#'   step_size     = 20
#' )
#'
#' # Limit to selected samples and chromosomes; recreate table if present
#' summarize_ch3_windows(
#'   ch3_db       = "my_db.ch3.db",
#'   samples      = c("Astrocytes","Blood_Plasma"),
#'   chrs         = c("chr1","chrX"),
#'   overwrite    = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{make_ch3_db}},
#' \code{\link{summarize_ch3_positions}},
#' \code{\link{summarize_ch3_regions}},
#' \code{\link{calc_ch3_diff}}
#'
#' @importFrom DBI dbExecute dbExistsTable dbGetQuery dbQuoteIdentifier dbRemoveTable
#' @importFrom glue glue
#' @export

summarize_ch3_windows <- function(ch3_db,
                                  input_table  = "calls",
                                  output_table = "windows",
                                  window_size = 1000,
                                  step_size = 10,
                                  chrs = c(as.character(1:22),
                                           paste0("chr", 1:22), "chrX", "chrY", "chrM",
                                           paste0("Chr", 1:22), "ChrX", "ChrY", "ChrM"),
                                  samples = NULL,             # NULL = all samples
                                  mod_code    = c("m", "h", "m + h"),
                                  unmod_code  = "-",
                                  unmod_label = "c",
                                  min_num_calls = 1,
                                  temp_dir = tempdir(),
                                  threads = NULL,             # default: all-but-one
                                  memory_limit = NULL,        # default: ~80% RAM
                                  overwrite = TRUE)
{
  start_time <- Sys.time()
  ch3_db <- MethylSeqR:::.ch3helper_connectDB(ch3_db)
  
  # Resource caps
  caps <- .auto_duckdb_resource_caps(0.80)
  thr  <- if (is.null(threads)) caps$threads else threads
  mem  <- if (is.null(memory_limit)) caps$memory_limit else memory_limit
  
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  DBI::dbExecute(ch3_db$con, sprintf("PRAGMA temp_directory='%s';", temp_dir))
  DBI::dbExecute(ch3_db$con, sprintf("PRAGMA memory_limit='%s';", mem))
  DBI::dbExecute(ch3_db$con, sprintf("PRAGMA threads=%d;", thr))
  
  in_id  <- as.character(DBI::dbQuoteIdentifier(ch3_db$con, input_table))
  out_id <- as.character(DBI::dbQuoteIdentifier(ch3_db$con, output_table))
  
  if (DBI::dbExistsTable(ch3_db$con, output_table) && overwrite)
    DBI::dbRemoveTable(ch3_db$con, output_table)
  
  # Sample list
  samp_query <- sprintf("SELECT DISTINCT sample_name FROM %s WHERE sample_name IS NOT NULL", in_id)
  all_samps <- DBI::dbGetQuery(ch3_db$con, samp_query)[,1]
  if (!is.null(samples)) {
    all_samps <- intersect(all_samps, samples)
  }
  if (length(all_samps) == 0) stop("No samples found.")
  
  specs   <- .parse_mod_specs(mod_code)
  cntsql  <- .build_pos_count_sql(unmod_code, unmod_label, specs)
  countsS <- cntsql$select_counts_pos
  labels  <- cntsql$labels_all
  
  # Prepare output schema with dynamic columns
  count_nulls <- paste(sprintf("CAST(NULL AS BIGINT) AS %s_counts", labels), collapse = ",\n      ")
  frac_nulls  <- paste(sprintf("CAST(NULL AS DOUBLE) AS %s_frac", labels),  collapse = ",\n      ")
  schema_sql <- glue::glue("
    CREATE TABLE IF NOT EXISTS {out_id} AS
    SELECT
      CAST(NULL AS VARCHAR) AS sample_name,
      CAST(NULL AS VARCHAR) AS chrom,
      CAST(NULL AS BIGINT)  AS start,
      CAST(NULL AS BIGINT)  AS \"end\",
      CAST(NULL AS BIGINT)  AS num_CpGs,
      CAST(NULL AS BIGINT)  AS num_calls,
      {count_nulls},
      {frac_nulls}
    WHERE 1=0;
  ")
  DBI::dbExecute(ch3_db$con, schema_sql)
  
  offsets <- seq(1, window_size - 1, by = step_size)
  chr_clause <- .chrom_filter_sql(chrs)
  
  # Process each sample separately
  for (samp in all_samps) {
    samp_esc <- gsub("'", "''", samp)
    
    pos_view <- glue::glue("
      CREATE OR REPLACE TEMP VIEW temp_positions AS
      SELECT
          sample_name,
          chrom,
          start,
          COUNT(*) AS num_calls,
          {countsS}
      FROM {in_id}
      WHERE sample_name = '{samp_esc}' AND start > 0{chr_clause}
      GROUP BY sample_name, chrom, start
    ")
    DBI::dbExecute(ch3_db$con, pos_view)
    
    # Dynamic SUMs for counts & fractions at window level
    sum_counts <- paste(sprintf("SUM(%s_counts) AS %s_counts", labels, labels),
                        collapse = ",\n            ")
    frac_cols <- paste(sprintf(
      "CASE WHEN SUM(num_calls) = 0 THEN NULL ELSE SUM(%s_counts) * 1.0 / SUM(num_calls) END AS %s_frac",
      labels, labels), collapse = ",\n            ")
    
    for (offset in offsets) {
      win_sql <- glue::glue("
        CREATE TEMP TABLE temp_table AS
        WITH window_map AS (
          SELECT 
            sample_name,
            chrom,
            start,
            (start - ((start - {offset}) % {window_size})) AS temp_start,
            num_calls,
            {paste(sprintf('%s_counts', labels), collapse = ', ')}
          FROM temp_positions
        )
        SELECT 
          sample_name,
          chrom,
          temp_start AS start,
          temp_start + {window_size} - 1 AS \"end\",
          COUNT(*)              AS num_CpGs,   -- one row per position
          SUM(num_calls)        AS num_calls,
          {sum_counts},
          {frac_cols}
        FROM window_map
        GROUP BY sample_name, chrom, temp_start
        HAVING SUM(num_calls) >= {min_num_calls};
      ")
      
      DBI::dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_table;")
      DBI::dbExecute(ch3_db$con, win_sql)
      
      ins_sql <- glue::glue("INSERT INTO {out_id} SELECT * FROM temp_table;")
      DBI::dbExecute(ch3_db$con, ins_sql)
      
      DBI::dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_table;")
    }
    
    DBI::dbExecute(ch3_db$con, "DROP VIEW IF EXISTS temp_positions;")
  }
  
  end_time <- Sys.time()
  message("Windows table created as ", output_table,
          " (", round(as.numeric(end_time - start_time, "mins"), 2), " min).")
  
  ch3_db$current_table <- output_table
  ch3_db <- MethylSeqR:::.ch3helper_cleanup(ch3_db)
  invisible(ch3_db)
}