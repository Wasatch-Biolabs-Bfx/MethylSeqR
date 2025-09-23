#' Summarize methylation by annotated regions
#'
#' Aggregates call-level rows from \code{input_table} (typically \code{"calls"})
#' into per-sample region summaries written to \code{output_table} (default
#' \code{"regions"}). Regions are provided via a BED/TSV/CSV file with columns
#' \code{chrom}, \code{start}, \code{end}, and optional \code{region_name}.
#' For each region the function computes:
#' number of CpG positions (rows), total calls, per-class counts, and per-class
#' fractions.
#'
#' @section How modification codes work:
#' Pass \code{mod_code} as single codes (e.g. \code{"m"}, \code{"h"}, \code{"a"})
#' or combinations with \code{"+"} (e.g. \code{"m + h"}). Labels are created by
#' removing spaces and \code{"+"} (e.g. \code{"m + h"} → \code{"mh"}). For each
#' label the table includes \code{<label>_counts} and \code{<label>_frac}. The
#' unmodified class is defined by \code{unmod_code} (default \code{"-"}), named
#' using \code{unmod_label} (default \code{"c"} → \code{c_counts}, \code{c_frac}).
#'
#' @param ch3_db Path to a \code{.ch3.db} DuckDB file or a \code{"ch3_db"} object.
#'   A connection is opened via internal helpers and cleaned on return.
#' @param input_table Source table containing call-level records (default \code{"calls"}).
#'   Must contain at least: \code{sample_name}, \code{chrom}, \code{start}, \code{end}, \code{call_code}.
#' @param output_table Destination table name (default \code{"regions"}).
#' @param region_file BED/TSV/CSV path with columns \code{chrom}, \code{start}, \code{end}
#'   and optional \code{region_name}. If missing, \code{region_name} is synthesized as
#'   \code{"chrom_start_end"}.
#' @param join Join type between positions and regions: one of \code{"inner"}, \code{"left"},
#'   or \code{"right"} (default \code{"inner"}).
#' @param chrs Character vector of chromosome filters; rows whose \code{chrom} match any
#'   value are retained. Defaults to common human aliases (1–22, \code{chrX}, \code{chrY}, \code{chrM}, …).
#' @param samples Optional character vector of \code{sample_name}s to include. If \code{NULL},
#'   all samples present in \code{input_table} are processed.
#' @param mod_code Character vector of modification specs to count (single codes or
#'   \code{"code1 + code2"} combinations). Default \code{c("m","h","m + h")}.
#' @param unmod_code Call code representing unmodified (default \code{"-"}).
#' @param unmod_label Label used to name unmodified columns (default \code{"c"}).
#' @param min_num_calls Minimum total calls required at the \emph{region} level to be written
#'   (default \code{1}). Regions below this threshold are skipped.
#' @param temp_dir Directory for DuckDB temporary files (default \code{tempdir()}).
#' @param threads Integer DuckDB thread count. If \code{NULL}, an internal heuristic
#'   (typically all-but-one core) is used.
#' @param memory_limit DuckDB memory limit string (e.g. \code{"16384MB"}).
#'   If \code{NULL}, an internal heuristic (~80\% of RAM) is used.
#' @param overwrite If \code{TRUE} and \code{output_table} exists, it is dropped before writing.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Reads \code{region_file} (CSV/TSV/BED) and normalizes columns. If \code{region_name}
#'         is absent, it is synthesized. Basic chromosome-prefix harmonization is performed when
#'         DB positions and annotation disagree on presence of a \code{"chr"} prefix.
#'   \item Configures DuckDB pragmas (\code{temp_directory}, \code{threads}, \code{memory_limit}).
#'   \item Builds per-position counts (one row per \code{sample_name}, \code{chrom}, \code{start})
#'         for the requested classes (\code{<label>_counts}, plus fractions later).
#'   \item Joins positions to regions using the chosen \code{join} type and aggregates per region:
#'         \code{num_CpGs}, \code{num_calls}, \code{<label>_counts}, and
#'         \code{<label>_frac} = \code{<label>_counts / num_calls}.
#' }
#'
#' @return (Invisibly) a \code{"ch3_db"} object pointing to the same DB file with
#'   \code{current_table} set to \code{output_table}. The created table has columns:
#'   \itemize{
#'     \item \code{sample_name}, \code{region_name}, \code{chrom}, \code{start}, \code{end},
#'           \code{num_CpGs}, \code{num_calls},
#'     \item for each label in \code{c(unmod_label, parsed(mod_code))}:
#'           \code{<label>_counts}, \code{<label>_frac}.
#'   }
#'
#' @examples
#' \dontrun{
#' # Default m/h summary by regions
#' summarize_ch3_regions(
#'   ch3_db      = "my_db.ch3.db",
#'   region_file = "islands_hg38.bed"
#' )
#'
#' # Novel 'a' code and m+h combination, left join to keep empty regions
#' summarize_ch3_regions(
#'   ch3_db      = "my_db.ch3.db",
#'   region_file = "islands_hg38.csv",
#'   mod_code    = c("a","m + h"),
#'   join        = "left",
#'   min_num_calls = 10
#' )
#' }
#'
#' @seealso
#' \code{\link{make_ch3_db}},
#' \code{\link{summarize_ch3_positions}},
#' \code{\link{summarize_ch3_windows}},
#' \code{\link{calc_ch3_diff}}
#'
#' @importFrom DBI dbExecute dbExistsTable dbGetQuery dbWriteTable dbQuoteIdentifier dbRemoveTable
#' @importFrom glue glue
#' @importFrom readr read_csv read_tsv
#' @importFrom tools file_ext
#' @export

summarize_ch3_regions <- function(ch3_db,
                                  input_table  = "calls",
                                  output_table = "regions",
                                  region_file,
                                  join = c("inner","left","right"),
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
  join <- match.arg(join)
  
  # ---- Read annotation (CSV/TSV/BED) ------------------------------------------
  ext <- tools::file_ext(region_file)
  if (tolower(ext) == "csv") {
    annotation <- readr::read_csv(
      region_file,
      col_names = c("chrom","start","end","region_name"),
      show_col_types = FALSE
    )
  } else if (tolower(ext) %in% c("tsv","bed")) {
    annotation <- readr::read_tsv(
      region_file,
      col_names = c("chrom","start","end","region_name"),
      show_col_types = FALSE
    )
  } else {
    stop("Invalid file type: ", ext, ". Only CSV, TSV, or BED are supported.")
  }
  
  # If file has a header row that duplicates colnames, drop it
  if (nrow(annotation) > 0 && tolower(annotation$chrom[1]) %in% c("chr","chrom")) {
    annotation <- annotation[-1, , drop = FALSE]
  }
  
  # Synthesize region_name if missing
  if (!"region_name" %in% names(annotation) || all(is.na(annotation$region_name))) {
    annotation$region_name <- paste(annotation$chrom, annotation$start, annotation$end, sep = "_")
  }
  
  # ---- Open DB and set resource caps ------------------------------------------
  ch3_db <- MethylSeqR:::.ch3helper_connectDB(ch3_db)
  
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
  
  # ---- Determine sample list ---------------------------------------------------
  samp_query <- sprintf("SELECT DISTINCT sample_name FROM %s WHERE sample_name IS NOT NULL", in_id)
  all_samps <- DBI::dbGetQuery(ch3_db$con, samp_query)[,1]
  if (!is.null(samples)) {
    all_samps <- intersect(all_samps, samples)
  }
  if (length(all_samps) == 0) stop("No samples found.")
  
  # ---- Build dynamic count SQL like positions/windows --------------------------
  specs   <- .parse_mod_specs(mod_code)
  cntsql  <- .build_pos_count_sql(unmod_code, unmod_label, specs)
  countsS <- cntsql$select_counts_pos
  labels  <- cntsql$labels_all
  
  # ---- Prepare output schema with dynamic columns ------------------------------
  count_nulls <- paste(sprintf("CAST(NULL AS BIGINT) AS %s_counts", labels), collapse = ",\n      ")
  frac_nulls  <- paste(sprintf("CAST(NULL AS DOUBLE) AS %s_frac",  labels), collapse = ",\n      ")
  schema_sql <- glue::glue("
    CREATE TABLE IF NOT EXISTS {out_id} AS
    SELECT
      CAST(NULL AS VARCHAR) AS sample_name,
      CAST(NULL AS VARCHAR) AS region_name,
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
  
  # ---- Upload annotation to temp table -----------------------------------------
  DBI::dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_annotation;")
  DBI::dbWriteTable(ch3_db$con, "temp_annotation", annotation, temporary = TRUE)
  
  # ---- Chromosome filter clause ------------------------------------------------
  chr_clause <- .chrom_filter_sql(chrs)
  
  # ---- Harmonize 'chr' prefix if needed (simple heuristic) ---------------------
  # Check whether positions have 'chr' prefix
  has_chr_positions <- DBI::dbGetQuery(ch3_db$con, sprintf("
    SELECT CASE WHEN EXISTS (
      SELECT 1 FROM %s WHERE chrom LIKE 'chr%%' LIMIT 1
    ) THEN 1 ELSE 0 END AS v
  ", in_id))$v[1] == 1
  
  has_chr_annot <- any(grepl("^chr", annotation$chrom))
  
  if (has_chr_positions && !has_chr_annot) {
    # add 'chr' to annotation
    DBI::dbExecute(ch3_db$con, "UPDATE temp_annotation SET chrom = 'chr' || CAST(chrom AS VARCHAR);")
  } else if (!has_chr_positions && has_chr_annot) {
    # strip 'chr' from annotation
    DBI::dbExecute(ch3_db$con, "UPDATE temp_annotation SET chrom = REGEXP_REPLACE(chrom, '^chr', '');")
  }
  
  # ---- Process per sample: position counts -> region aggregation ---------------
  for (samp in all_samps) {
    samp_esc <- gsub("'", "''", samp)
    
    # Per-position counts for this sample
    pos_view <- glue::glue("
      CREATE OR REPLACE TEMP VIEW temp_positions AS
      SELECT
          sample_name,
          chrom,
          start,
          \"end\",
          COUNT(*) AS num_calls,
          {countsS}
      FROM {in_id}
      WHERE sample_name = '{samp_esc}' AND start > 0{chr_clause}
      GROUP BY sample_name, chrom, start, \"end\"
    ")
    DBI::dbExecute(ch3_db$con, pos_view)
    
    # Build dynamic SUMs and fractions at region level
    sum_counts <- paste(sprintf("COALESCE(SUM(p.%s_counts), 0) AS %s_counts", labels, labels),
                        collapse = ",\n          ")
    frac_cols <- paste(sprintf(
      "CASE WHEN COALESCE(SUM(p.num_calls),0) = 0 THEN NULL ELSE COALESCE(SUM(p.%s_counts),0) * 1.0 / COALESCE(SUM(p.num_calls),0) END AS %s_frac",
      labels, labels), collapse = ",\n          ")
    
    # Choose JOIN keyword
    join_kw <- switch(join,
                      inner = "JOIN",
                      left  = "LEFT JOIN",
                      right = "RIGHT JOIN")
    
    reg_sql <- glue::glue("
      CREATE TEMP TABLE temp_regions AS
      SELECT
        p.sample_name,
        a.region_name,
        a.chrom,
        a.start,
        a.\"end\",
        COUNT(*)                       AS num_CpGs,
        COALESCE(SUM(p.num_calls), 0)  AS num_calls,
        {sum_counts},
        {frac_cols}
      FROM temp_positions p
      {join_kw} temp_annotation a
        ON p.chrom = a.chrom
       AND CAST(p.start AS DOUBLE) BETWEEN CAST(a.start AS DOUBLE) AND CAST(a.\"end\" AS DOUBLE)
      GROUP BY p.sample_name, a.region_name, a.chrom, a.start, a.\"end\"
      HAVING COALESCE(SUM(p.num_calls), 0) >= {min_num_calls};
    ")
    
    DBI::dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_regions;")
    DBI::dbExecute(ch3_db$con, reg_sql)
    DBI::dbExecute(ch3_db$con, glue::glue("INSERT INTO {out_id} SELECT * FROM temp_regions;"))
    DBI::dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_regions;")
    DBI::dbExecute(ch3_db$con, "DROP VIEW IF EXISTS temp_positions;")
  }
  
  # ---- Finish ------------------------------------------------------------------
  message("Regions table created as ", output_table,
          " (", round(as.numeric(Sys.time() - start_time, "secs"), 1), "s).")
  
  ch3_db$current_table <- output_table
  ch3_db <- MethylSeqR:::.ch3helper_cleanup(ch3_db)
  invisible(ch3_db)
}
