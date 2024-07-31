#' Open Methylation Data from a Directory
#'
#' This function opens methylation data stored in a specified directory and applies filters
#' based on sample names, chromosomes, read length, base quality, and call probability.
#'
#' @param path Path to the database directory containing methylation data.
#'
#' @param samples A vector of sample names to extract methylation data from.
#'                Default is "all", which retrieves data from all available samples.
#'
#' @param chrs A vector of chromosome names to extract. Default is "all",
#'              which retrieves data from all available chromosomes.
#'
#' @param min_read_length Minimum read length to allow. Filters out any reads
#'                        shorter than this length. Default is 100.
#'
#' @param min_call_prob Minimum probability for base call. Default is 0.9.
#'
#' @param min_base_qual Minimum Phred quality score for base calls. Default is 10.
#'
#' @param max_memory Maximum memory allocation for the database connection.
#'                   Default is "8GB".
#'
#' @param max_threads Maximum number of CPU threads to use for processing. Default is 4.
#'
#' @return A DuckDB object containing all methylation data that meets the specified
#'         filtering criteria.
#'
#' @examples
#' data = open_dat("path/to/database")
#' data = open_dat("path/to/database", samples = c("sample1", "sample2"), chrs = c("chr1", "chr2"))
#'
#' @importFrom DBI dbConnect dbExecute
#' @import duckdb
#'
#' @export
open_dat <- function(path,
                     samples = "all",
                     chrs = "all",
                     min_read_length = 100,
                     min_call_prob = .9,
                     min_base_qual = 10,
                     max_memory = "8GB",
                     max_threads = 4)
{
  # Perform checks and calcs for arguments
  stopifnot(
    "A single folder must be provided" =
      length(path) == 1)

  if (all(samples == "all"))
    samples <-
      list.dirs(path = path,
                full.names = FALSE, recursive = FALSE)

  if (all(chrs == "all"))
    chrs <-
      list.dirs(path = paste0(path, "/", samples),
                full.names = FALSE, recursive = FALSE) |>
      unique()

  stopifnot(
    "Could not find one or more samples in your path" =
      dir.exists(paste0(path, "/", samples)))

  # Set up db connection
  db_con <- dbConnect(duckdb(tempfile()))

  # Set memory and processor limits
  dbExecute(
    db_con,
    sprintf("SET memory_limit = '%s';
             SET threads TO %i;",
            max_memory, max_threads))

  # open (and return) data
  open_dataset(
    path,
    partitioning = c("sample_name", "chrom")) |>
    filter(
      sample_name %in% samples,
      chrom %in% chrs,
      call_prob >= min_call_prob,
      read_length >= min_read_length,
      base_qual >= min_base_qual) |>
    to_duckdb(con = db_con)
}
