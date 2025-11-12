#' Create a CH3 DuckDB from Parquet CH3 files (with optional sample naming)
#'
#' Build a DuckDB database containing filtered methylation call data from one or
#' more \code{.ch3} parquet files. Inputs may be individual files, directories
#' (expanded to all \code{*.ch3} files), or a **named character vector** where
#' names are used as \code{sample_name}s in the output.
#'
#' @section Input forms:
#' \itemize{
#'   \item \strong{Files}: \code{c("a.ch3", "b.ch3")}
#'   \item \strong{Directories}: \code{c("dir_of_ch3s/")} (expands to all \code{*.ch3})
#'   \item \strong{Named files/dirs}: \code{c(SampleA = "a.ch3", SampleB = "dir/")}
#'         — names become \code{sample_name}. If a name is not provided, the
#'         filename stem (without \code{.ch3}) is used.
#' }
#'
#' @param ch3_files Character vector of CH3 parquet file paths and/or directories.
#'   May be a \emph{named} vector to assign \code{sample_name}s explicitly; any
#'   entry of the form \code{NAME=PATH} is also accepted. Directories are scanned
#'   (non-recursively) for \code{*.ch3} files. Must not be empty.
#' @param db_name Path (without or with \code{.ch3.db} extension) for the DuckDB
#'   database to be created; \code{.ch3.db} is appended if missing.
#' @param chrom Optional chromosome filter. Either a single string (e.g.,
#'   \code{"chr1"}) or a character vector (e.g., \code{c("chr1","chr2","chrX")}).
#'   If \code{NULL}, all chromosomes are included.
#' @param min_read_length Minimum read length to keep (default \code{50}).
#' @param min_call_prob Minimum call probability to keep (default \code{0.9}).
#' @param min_base_qual Minimum base quality to keep (default \code{10}).
#' @param flag Optional numeric flag value to require; if \code{NULL}, no flag filter.
#'
#' @details
#' \strong{What it does}
#' \itemize{
#'   \item Expands \code{ch3_files} (handling directories and named entries) into a mapping
#'         of source files and optional \code{sample_name}s.
#'   \item Configures DuckDB pragmas for temp directory, thread count (all-but-one core),
#'         and a memory limit (~50% of detected RAM).
#'   \item Drops any existing tables in the target DB.
#'   \item Reads all input \code{.ch3} parquet files in a single pass and creates a
#'         table \code{calls} with columns:
#'         \code{sample_name}, \code{chrom}, \code{start}, \code{end}, \code{read_position},
#'         \code{call_code}, \code{read_length}, \code{call_prob}, \code{base_qual}, \code{flag}.
#'         When names are not given for inputs, \code{sample_name} defaults to the file stem.
#'   \item Applies pushdown filters based on \code{chrom}, \code{min_read_length},
#'         \code{min_call_prob}, \code{min_base_qual}, and \code{flag}.
#' }
#'
#' \strong{Side effects and performance}
#' \itemize{
#'   \item Creates (or overwrites) a DuckDB file at \code{db_name}.
#'   \item Uses a temp directory for DuckDB spills under \code{tempdir()}.
#'   \item A temporary in-memory table \code{file_map} may be created for input mapping.
#' }
#'
#' @return (Invisibly) a list of class \code{"ch3_db"} with elements:
#'   \itemize{
#'     \item \code{db_file}: path to the created DuckDB file,
#'     \item \code{current_table}: \code{NULL} (set by downstream functions),
#'     \item \code{con}: connection is closed by cleanup and set to \code{"none"}.
#'   }
#'   The database contains at least the \code{calls} table.
#'
#' @examples
#' \dontrun{
#' # 1) Directory of CH3 files (non-recursive scan for *.ch3)
#' make_ch3_db(ch3_files = "path/to/ch3_dir",
#'             db_name   = "my_db")
#'
#' # 2) Explicit files (auto-sample names from stems)
#' make_ch3_db(ch3_files = c("A.ch3", "B.ch3"),
#'             db_name   = "two_samples.ch3.db",
#'             min_read_length = 100, min_base_qual = 10)
#'
#' # 3) Named inputs (sample_name set from names)
#' make_ch3_db(
#'   ch3_files = c(
#'     Sample1      = "../CH3/Sample1.ch3",
#'     Sample2  = "../CH3/Sample2.ch3",
#'     Sample3   = "../CH3/Sample3.combined.ch3"
#'   ),
#'   db_name = "My_DB",
#'   min_base_qual = 10,
#'   min_read_length = 100
#' )
#'
#' # 4) Filter to specific chromosomes
#' make_ch3_db(
#'   ch3_files = c(S1 = "A.ch3", S2 = "B.ch3"),
#'   db_name   = "chr1_chrX_only",
#'   chrom     = c("chr1","chrX")
#' )
#' }
#'
#' @seealso
#' \code{\link{summarize_ch3_positions}}, \code{\link{summarize_ch3_regions}},
#' \code{\link{summarize_ch3_windows}}, \code{\link{get_ch3_dbinfo}},
#' \code{\link{get_ch3_tableinfo}}, \code{\link{calc_ch3_diff}}
#'
#' @importFrom DBI dbConnect dbExecute dbListTables dbWriteTable dbQuoteIdentifier
#' @importFrom duckdb duckdb
#' @importFrom glue glue
#' @importFrom parallel detectCores
#'
#' @export


make_ch3_db <- function(ch3_files, 
                        db_name,
                        chrom = NULL,                 # can be a single value or a character vector
                        min_read_length = 50, 
                        min_call_prob  = 0.9, 
                        min_base_qual  = 10, 
                        flag = NULL)
{
  start_time <- Sys.time()
  if (length(ch3_files) == 0) stop("No files provided.")
  
  # --- Accept named vectors: convert names to "name=path" entries -----------------
  if (!is.null(names(ch3_files))) {
    nm <- names(ch3_files)
    named_idx <- nzchar(nm)
    ch3_files[named_idx] <- paste0(nm[named_idx], "=", ch3_files[named_idx])
    names(ch3_files) <- NULL
  }
  
  # --- resolve entries: expand dirs and parse `name=path` -------------------------
  parse_entry <- function(entry) {
    if (grepl("=", entry, fixed = TRUE)) {
      parts <- strsplit(entry, "=", fixed = TRUE)[[1]]
      nm <- trimws(parts[1])
      p  <- trimws(paste(parts[-1], collapse = "="))  # keep '=' if in path
    } else {
      nm <- NA_character_
      p  <- entry
    }
    list(name = nm, path = p)
  }
  entries <- lapply(ch3_files, parse_entry)
  
  expand_paths <- function(nm, p) {
    if (!file.exists(p)) stop("Path does not exist: ", p)
    if (dir.exists(p)) {
      fs <- list.files(p, pattern = "\\.ch3$", full.names = TRUE, recursive = FALSE)
      if (length(fs) == 0) stop("No .ch3 files in directory: ", p)
      data.frame(sample_name = if (!is.na(nm)) nm else NA_character_,
                 file = normalizePath(fs, winslash = "/"),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(sample_name = nm,
                 file = normalizePath(p, winslash = "/"),
                 stringsAsFactors = FALSE)
    }
  }
  
  df_files <- do.call(
    rbind,
    Map(expand_paths, lapply(entries, `[[`, "name"), lapply(entries, `[[`, "path"))
  )
  if (nrow(df_files) == 0) stop("No .ch3 files found.")
  
  # --- derive robust fallback names in R (not SQL) -------------------------------
  derive_sample_from_path <- function(p) {
    b <- basename(p)
    # strip case-insensitive ".ch3"
    stem <- sub("(?i)\\.ch3$", "", b, perl = TRUE)
    # then strip a trailing "-<digits>" if present
    stem <- sub("-[0-9]+$", "", stem, perl = TRUE)
    stem
  }
  
  need_name <- is.na(df_files$sample_name) | !nzchar(df_files$sample_name)
  df_files$sample_name[need_name] <- vapply(df_files$file[need_name], derive_sample_from_path, character(1))
  
  # also store basename to make the SQL JOIN robust across path formats
  df_files$base <- basename(df_files$file)
  
  # --- DB setup -------------------------------------------------------------------
  if (!grepl("\\.ch3\\.db$", db_name)) db_name <- paste0(db_name, ".ch3.db")
  cat("Building Database...\n")
  
  ch3_db <- list(db_file = db_name, current_table = NULL, con = NULL)
  class(ch3_db) <- "ch3_db"
  ch3_db$con <- DBI::dbConnect(duckdb::duckdb(ch3_db$db_file), read_only = FALSE)
  
  # Resource caps: fast spill dir, all-but-one cores, ~50% RAM (absolute units)
  detect_caps <- function(frac = 0.80) {
    cores <- tryCatch(parallel::detectCores(), error = function(e) 2)
    threads <- max(1, cores - 1)
    total_bytes <- NA_real_
    if (file.exists("/proc/meminfo")) {
      ln <- tryCatch(readLines("/proc/meminfo"), error = function(e) "")
      mt <- sub(".*MemTotal:\\s*([0-9]+) kB.*", "\\1", grep("^MemTotal:", ln, value = TRUE))
      if (nzchar(mt)) total_bytes <- as.numeric(mt) * 1024
    } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
      out <- suppressWarnings(system("sysctl -n hw.memsize", intern = TRUE))
      if (length(out)) total_bytes <- as.numeric(out)
    } else if (.Platform$OS.type == "windows") {
      out <- suppressWarnings(system("wmic computersystem get TotalPhysicalMemory /value", intern = TRUE))
      if (length(out)) {
        val <- sub(".*TotalPhysicalMemory=([0-9]+).*", "\\1", paste(out, collapse = ""))
        if (grepl("^[0-9]+$", val)) total_bytes <- as.numeric(val)
      }
    }
    if (!is.finite(total_bytes)) total_bytes <- 8 * 1024^3
    limit_mb <- max(1024, floor(frac * total_bytes / 1024^2))
    list(threads = threads, memory_limit = paste0(limit_mb, "MB"))
  }
  
  # Use 50% of RAM (returns e.g. "16384MB")
  caps <- detect_caps(0.50)
  
  tmp_dir <- file.path(tempdir(), "duckdb_tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  DBI::dbExecute(ch3_db$con, sprintf("PRAGMA temp_directory='%s';", tmp_dir))
  DBI::dbExecute(ch3_db$con, sprintf("PRAGMA threads=%d;", caps$threads))
  DBI::dbExecute(ch3_db$con, sprintf("PRAGMA memory_limit='%s';", caps$memory_limit))
  
  # Drop existing tables
  for (tbl in DBI::dbListTables(ch3_db$con)) {
    DBI::dbExecute(ch3_db$con, paste0("DROP TABLE ", DBI::dbQuoteIdentifier(ch3_db$con, tbl)))
  }
  
  # --- small mapping table (path -> optional sample_name) -------------------------
  DBI::dbWriteTable(ch3_db$con, "file_map", df_files, temporary = TRUE, overwrite = TRUE)
  
  # --- filters & columns (pushdown) ----------------------------------------------
  esc <- function(x) gsub("'", "''", x)
  filters <- c()
  if (!is.null(chrom)) {
    if (length(chrom) == 1) {
      filters <- c(filters, paste0("chrom = '", esc(chrom), "'"))
    } else {
      filters <- c(filters, paste0("chrom IN (", paste(sprintf("'%s'", esc(chrom)), collapse = ", "), ")"))
    }
  }
  filters <- c(filters, paste0("read_length >= ", as.numeric(min_read_length)))
  filters <- c(filters, paste0("call_prob  >= ", as.numeric(min_call_prob)))
  filters <- c(filters, paste0("base_qual  >= ", as.numeric(min_base_qual)))
  if (!is.null(flag)) filters <- c(filters, paste0("flag = ", as.numeric(flag)))
  where_clause <- if (length(filters)) paste("WHERE", paste(filters, collapse = " AND ")) else ""
  
  # ⬅️ Include read_position here
  wanted_sql <- paste(c("read_id", "chrom","start","\"end\"", "read_position",
                        "call_code","read_length","call_prob","base_qual","flag"),
                      collapse = ", ")
  
  # Paths for read_parquet
  file_list_sql <- paste0("['", paste(esc(df_files$file), collapse = "','"), "']")
  
  # --- one parallel scan; join to mapping; fallback name from filename -----------
  sql <- glue::glue("
    CREATE TABLE calls AS
    WITH src AS (
      SELECT *
      FROM read_parquet({file_list_sql}, filename = TRUE)
    ),
    tagged AS (
      SELECT
        -- use the sample_name we computed in R
        m.sample_name AS sample_name,
        {wanted_sql}
      FROM src s
      LEFT JOIN file_map m
        ON m.file = s.filename
        OR m.base = REGEXP_REPLACE(s.filename, '^.*[/\\\\\\\\]', '')
    )
    SELECT *
    FROM tagged
    {where_clause}
  ")
  DBI::dbExecute(ch3_db$con, sql)
  
  # --- finish --------------------------------------------------------------------
  total_seconds <- as.numeric(Sys.time() - start_time, units = "secs")
  if (total_seconds > 60) {
    message("Database created at ", ch3_db$db_file,
            "\nTime elapsed: ", round(total_seconds/60, 2), " minutes\n")
  } else {
    message("Database created at ", ch3_db$db_file,
            "\nTime elapsed: ", round(total_seconds, 2), " seconds\n")
  }
  
  ch3_db <- MethylSeqR:::.ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}