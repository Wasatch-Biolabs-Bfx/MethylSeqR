#' Summarize per-position methylation calls into a DuckDB table
#'
#' Aggregates row-level CH3 calls (from \code{input_table}, typically \code{"calls"})
#' into a per-sample, per-genomic-position table (\code{output_table}, default \code{"positions"}).
#' Counts and fractions are created for the unmodified class and for any set of
#' user-defined modification codes and/or combinations (e.g., \code{"m"}, \code{"h"},
#' \code{"m + h"}, \code{"a"}).
#'
#' @section How modification codes work:
#' Pass \code{mod_code} as character values (single codes) or \code{"code1 + code2"} to
#' build combinations. Column labels are derived by removing spaces and \code{"+"}. For example:
#' \itemize{
#'   \item \code{mod_code = c("m", "h", "m + h")} produces columns \code{m_counts}, \code{h_counts}, \code{mh_counts}
#'   \item \code{mod_code = "a"} produces \code{a_counts}
#' }
#' The unmodified class is controlled by \code{unmod_code} (e.g., \code{"-"}), with a column
#' name prefix set by \code{unmod_label} (default \code{"c"} for \code{c_counts}, \code{c_frac}).
#'
#' @param ch3_db A path to a \code{.ch3.db} DuckDB file or a \code{"ch3_db"} object; a connection
#'   is opened via internal helpers and closed on return.
#' @param input_table Name of the source table containing call-level records (default \code{"calls"}).
#'   Must include at least: \code{sample_name}, \code{chrom}, \code{start}, \code{call_code}.
#' @param output_table Name of the destination positions table to create/extend (default \code{"positions"}).
#' @param chrs Character vector of chromosome/name filters. Rows whose \code{chrom} match any value
#'   are retained. Defaults to common human aliases (1–22, \code{chrX}, \code{chrY}, \code{chrM}, etc.).
#' @param samples Optional character vector of \code{sample_name}s to include. If \code{NULL},
#'   all samples present in \code{input_table} are processed.
#' @param mod_code Character vector of modification specifications to count. Each entry is either
#'   a single code (e.g., \code{"m"}, \code{"h"}, \code{"a"}) or a \code{"code1 + code2"} combination.
#'   Default: \code{c("m","h","m + h")}.
#' @param unmod_code The call code representing unmodified calls (default \code{"-"}).
#' @param unmod_label Label used to name unmodified columns (default \code{"c"} yielding
#'   \code{c_counts} and \code{c_frac}).
#' @param min_num_calls Minimum total calls required at a position to be written (default \code{1}).
#' @param temp_dir Directory for DuckDB temporary files (default \code{tempdir()}).
#' @param threads Integer DuckDB thread count. If \code{NULL}, uses an internal heuristic
#'   (typically all-but-one core).
#' @param memory_limit DuckDB memory limit string (e.g., \code{"16384MB"}). If \code{NULL},
#'   an internal heuristic (~80\% of RAM) is used.
#' @param overwrite If \code{TRUE} and \code{output_table} exists, it is dropped before writing.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Opens a DB connection and configures DuckDB pragmas (\code{temp_directory}, \code{threads}, \code{memory_limit}).
#'   \item Determines the list of samples to process (all by default, or the intersection with \code{samples}).
#'   \item Builds dynamic SQL to compute, per sample and position (\code{chrom}, \code{start}):
#'         \code{num_calls}, \code{<label>_counts} for each requested code/combination, and
#'         \code{<label>_frac} = \code{<label>_counts / num_calls}.
#'   \item Pre-creates the \code{output_table} schema with the appropriate dynamic columns, then
#'         inserts rows per sample. \code{end} is set equal to \code{start}.
#' }
#' Positions with \code{num_calls < min_num_calls} are skipped. Chromosome filtering is done in SQL.
#'
#' @return (Invisibly) a \code{"ch3_db"} object pointing to the same DB file with
#'   \code{current_table} set to \code{output_table}. The created table has columns:
#'   \itemize{
#'     \item \code{sample_name}, \code{chrom}, \code{start}, \code{end},
#'           \code{num_calls},
#'     \item for each label in \code{c(unmod_label, parsed(mod_code))}: \code{<label>_counts}, \code{<label>_frac}.
#'   }
#'
#' @examples
#' \dontrun{
#' # From a calls table, build per-position summaries for default m/h classes:
#' summarize_ch3_positions(
#'   ch3_db       = "my_db.ch3.db",
#'   input_table  = "calls",
#'   output_table = "positions"
#' )
#'
#' # Restrict to specific samples and chromosomes, require at least 5 calls:
#' summarize_ch3_positions(
#'   ch3_db       = "my_db.ch3.db",
#'   samples      = c("Cortical_Neurons","Astrocytes"),
#'   chrs         = c("chr1","chrX"),
#'   min_num_calls = 5
#' )
#'
#' # Count a novel code "a" and an m+h combination; rename unmodified to 'u':
#' summarize_ch3_positions(
#'   ch3_db       = "my_db.ch3.db",
#'   mod_code     = c("a", "m + h"),
#'   unmod_code   = "-",
#'   unmod_label  = "u"
#' )
#'
#' # Recreate positions table (drop if exists):
#' summarize_ch3_positions(
#'   ch3_db       = "my_db.ch3.db",
#'   overwrite    = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{make_ch3_db}}, \code{\link{summarize_ch3_regions}},
#' \code{\link{summarize_ch3_windows}}, \code{\link{calc_ch3_diff}}
#'
#' @importFrom DBI dbExecute dbExistsTable dbGetQuery dbQuoteIdentifier dbRemoveTable
#' @importFrom glue glue
#'
#' @export

summarize_ch3_positions <- function(ch3_db,
                                    input_table  = "calls",
                                    output_table = "positions",
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
  
  # Determine sample list
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
  
  # Pre-create schema once with dynamic columns
  count_nulls <- paste(sprintf("CAST(NULL AS BIGINT) AS %s_counts", labels), collapse = ",\n      ")
  frac_nulls  <- paste(sprintf("CAST(NULL AS DOUBLE) AS %s_frac", labels),  collapse = ",\n      ")
  schema_sql <- glue::glue("
    CREATE TABLE IF NOT EXISTS {out_id} AS
    SELECT
      CAST(NULL AS VARCHAR) AS sample_name,
      CAST(NULL AS VARCHAR) AS chrom,
      CAST(NULL AS BIGINT)  AS start,
      CAST(NULL AS BIGINT)  AS \"end\",
      CAST(NULL AS BIGINT)  AS num_calls,
      {count_nulls},
      {frac_nulls}
    WHERE 1=0;
  ")
  DBI::dbExecute(ch3_db$con, schema_sql)
  
  chr_clause <- .chrom_filter_sql(chrs)
  
  # Process each sample separately (chunking by sample)
  for (samp in all_samps) {
    samp_esc <- gsub("'", "''", samp)
    
    view_sql <- glue::glue("
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
    DBI::dbExecute(ch3_db$con, view_sql)
    
    frac_cols <- paste(sprintf(
      "CASE WHEN num_calls = 0 THEN NULL ELSE %s_counts * 1.0 / num_calls END AS %s_frac",
      labels, labels), collapse = ",\n          ")
    
    insert_sql <- glue::glue("
      INSERT INTO {out_id}
      SELECT
        sample_name,
        chrom,
        start,
        start AS \"end\",
        num_calls,
        {paste(sprintf('%s_counts', labels), collapse = ', ')},
        {frac_cols}
      FROM temp_positions
      WHERE num_calls >= {min_num_calls};
    ")
    DBI::dbExecute(ch3_db$con, insert_sql)
    
    DBI::dbExecute(ch3_db$con, "DROP VIEW IF EXISTS temp_positions;")
  }
  
  end_time <- Sys.time()
  message("Positions table created as ", output_table,
          " (", round(as.numeric(end_time - start_time, "secs"), 1), "s).")
  
  ch3_db$current_table <- output_table
  ch3_db <- MethylSeqR:::.ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}