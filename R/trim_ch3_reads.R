#' Trim methylation calls at the read ends in a ch3 database
#'
#' Creates a new table in a ch3 database that keeps only methylation calls
#' within a central fraction of each read, optionally trimming from the
#' start, the end, or both. Trimming is expressed as a fraction of the
#' read length and applied via SQL filters on `read_position`.
#'
#' For a read of length `L` and `by_frac = f`:
#' \itemize{
#'   \item If \code{trim_start = TRUE}, positions with
#'         \code{read_position < L * f} are dropped.
#'   \item If \code{trim_end   = TRUE}, positions with
#'         \code{read_position >= L * (1 - f)} are dropped.
#' }
#' For example, with \code{read_length = 100} and \code{by_frac = 0.1},
#' positions \code{0–9} and/or \code{90–99} are trimmed depending on
#' which end(s) are enabled.
#'
#' @param ch3_db An object identifying the ch3 database, as accepted by
#'   \code{MethylSeqR:::.ch3helper_connectDB()} (e.g. a file path or an
#'   existing ch3 database object).
#' @param by_frac Numeric scalar between 0 and 1 (exclusive) giving the
#'   fraction of the read length to trim from each enabled end. For
#'   example, \code{by_frac = 0.1} trims 10\% of the read length from
#'   the start and/or end.
#' @param trim_start Logical; whether to trim from the start (low
#'   \code{read_position}) of each read.
#' @param trim_end Logical; whether to trim from the end (high
#'   \code{read_position}) of each read.
#' @param input_table Character scalar; name of the input table
#'   containing methylation calls (default \code{"calls"}). Must contain
#'   at least \code{read_position} and \code{read_length} columns.
#' @param output_table Character scalar; name of the output table to
#'   create/overwrite with the trimmed calls
#'   (default \code{"calls_trimmed"}).
#'
#' @details
#' This function issues a \code{CREATE OR REPLACE TABLE} SQL statement
#' against the ch3 database. The resulting table contains a subset of
#' rows from \code{input_table} where the read positions fall within the
#' untrimmed interior region defined by \code{by_frac}, \code{trim_start},
#' and \code{trim_end}.
#'
#' If \code{by_frac} is not strictly between 0 and 1, an error is thrown.
#' If both \code{trim_start} and \code{trim_end} are \code{FALSE}, the
#' function also errors, as there would be nothing to trim.
#'
#' On success, \code{ch3_db$current_table} is updated to \code{output_table}
#' and a short timing message is printed.
#'
#' @return
#' Invisibly returns the updated \code{ch3_db} object, with
#' \code{current_table} set to \code{output_table}. The function is
#' called for its side effects of creating/replacing the trimmed table
#' in the database.
#'
#' @importFrom DBI dbExecute
#' @importFrom glue glue
#'
#' @examples
#' #' \dontrun{
#' # Trim the last 10% of each read and write to a new table
#' ch3_db <- "my_sample.ch3.sqlite"
#'
#' ch3_db <- trim_ch3_reads(
#'   ch3_db       = ch3_db,
#'   by_frac      = 0.1,
#'   trim_start   = FALSE,
#'   trim_end     = TRUE,
#'   input_table  = "calls",
#'   output_table = "calls_trimmed_10pct_end"
#' )
#'
#' # Trim 5% from both the start and end of each read
#' ch3_db <- trim_ch3_reads(
#'   ch3_db       = ch3_db,
#'   by_frac      = 0.05,
#'   trim_start   = TRUE,
#'   trim_end     = TRUE,
#'   input_table  = "calls",
#'   output_table = "calls_trimmed_5pct_both"
#' )
#' }
#'
#' @export
trim_ch3_reads <- function(ch3_db, 
                           by_frac = .1,
                           trim_start = FALSE,
                           trim_end   = TRUE,
                           input_table  = "calls",
                           output_table = "calls_trimmed") {
  start_time <- Sys.time()
  
  # sanity check
  if (by_frac <= 0 || by_frac >= 1) {
    stop("by_frac must be between 0 and 1 (e.g. 0.1 for trimming last 10%).")
  }
  
  ch3_db <- MethylSeqR:::.ch3helper_connectDB(ch3_db)
  
  # Build filter conditions
  conds <- c()
  
  if (trim_start) {
    # e.g. read_length=100, by_frac=.1 → drop read_position < 10
    conds <- c(conds, glue::glue("read_position >= read_length * {by_frac}"))
  }
  
  if (trim_end) { #  SQL: keep only calls that start *before* the trimmed tail
    # e.g. read_length=100, by_frac=.1 → drop read_position >= 90
    conds <- c(conds, glue::glue("read_position < read_length * (1 - {by_frac})"))
  }
  
  # If neither is TRUE, it's a no-op
  if (length(conds) == 0) {
    stop("Both trim_start and trim_end are set to FALSE — nothing to trim.")
  }
  
  where_clause <- paste(conds, collapse = " AND ")
  
  query <- glue::glue("
    CREATE OR REPLACE TABLE {output_table} AS
    SELECT *
    FROM {input_table}
    WHERE {where_clause}
  ")
  
  DBI::dbExecute(ch3_db$con, query)
  
  elapsed <- Sys.time() - start_time
  message(sprintf(
    "trim_ch3_reads: wrote '%s' (by_frac = %.3f, start=%s, end=%s) in %.2f sec",
    output_table, by_frac, trim_start, trim_end,
    as.numeric(elapsed, units = 'secs')
  ))
  
  ch3_db$current_table <- output_table
  
  ch3_db <- MethylSeqR:::.ch3helper_closeDB(ch3_db)
  
  invisible(ch3_db)
  
  
}