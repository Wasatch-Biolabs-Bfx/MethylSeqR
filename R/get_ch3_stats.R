#' Compute Statistics for a Ch3 File
#'
#' This function computes various statistics for a given Ch3 file stored in 
#' Parquet format, including total calls, total reads, mean read length, 
#' CpG coverage at different read thresholds, flag distributions, and high-quality 
#' call counts based on probability and base quality thresholds.
#'
#' @param ch3_file Character. Path to the Ch3 file in Parquet format.
#' @param log_file Character (optional). Path to a log file where the output 
#' will be written. If NULL, results are printed to the console.
#' @param min_reads Numeric vector. A set of thresholds for reporting CpG coverage 
#' at different minimum read counts. Default: \code{c(1, 5, 10, 15)}.
#' @param call_prob_threshold Numeric. The minimum modification probability 
#' to consider a high-confidence call. Default: \code{0.9}.
#' @param base_qual_threshold Numeric. The minimum base quality required 
#' to count as a high-quality call. Default: \code{10}.
#' @param silently Logical. If \code{TRUE}, suppresses console output. 
#' Default: \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{num_calls}}{Total number of modification calls in the file.}
#'   \item{\code{num_reads}}{Total number of unique reads in the file.}
#'   \item{\code{cpg_coverage}}{A matrix with CpG coverage counts at different 
#'   \code{min_reads} thresholds.}
#'   \item{\code{flag_counts}}{A data frame with the count and percentage of 
#'   calls per flag value.}
#'   \item{\code{high_conf_calls}}{A data frame with the count and percentage 
#'   of calls with modification probability above \code{call_prob_threshold}.}
#'   \item{\code{high_qual_calls}}{A data frame with the count and percentage 
#'   of calls with base quality above \code{base_qual_threshold}.}
#'   \item{\code{avg_read_length}}{The mean read length.}
#' }
#'
#' @details 
#' The function reads the Ch3 file as a dataset and computes basic statistics 
#' about the calls and reads. It also evaluates CpG coverage based on different 
#' read count thresholds, distributions of flag values, and counts of high-quality 
#' calls based on user-defined probability and quality thresholds.
#'
#' If a \code{log_file} is provided, the results are written to it. Otherwise, 
#' they are printed to the console unless \code{silently = TRUE}.
#'
#' 
#' @examples
#' \dontrun{
#' get_ch3_stats("example.ch3.parquet")
#' get_ch3_stats("example.ch3.parquet", log_file = "stats.log", silently = TRUE)
#' }
#'
#' @importFrom dplyr select summarise pull n_distinct count mutate collect
#' @importFrom arrow open_dataset
#' 
#' @export

get_ch3_stats <- function(ch3_file, 
                          log_file = NULL,
                          min_reads = c(1, 5, 10, 15), 
                          call_prob_threshold = 0.9, 
                          base_qual_threshold = 10,
                          silently = FALSE) 
{
  # Write to log_file if provided
  if (!is.null(log_file)) {
    sink(log_file, split = TRUE)
    on.exit(sink(), add = TRUE)
  }
  
  dataset <- open_dataset(ch3_file, format = "parquet")
  
  cat("=================================================\n",
      "               Ch3 File Statistics               \n",
      "=================================================\n", 
      "File: ", ch3_file, "\n",
      sep = "")
  
  # Open the Parquet file as a table
  num_rows <- 
    dataset |>
    summarize(
      n = n()) |>
    pull(as_vector = TRUE)
  
  # Unique read_id count
  num_reads <- 
    dataset |> 
    select(read_id) |>
    summarise(
      n = n_distinct(read_id)) |> 
    pull(as_vector = TRUE)
  
  if (!silently) {
    cat("\nGeneral Statistics:\n",
        "  Total Calls:\t\t", prettyNum(num_rows, big.mark = ","), "\n",
        "  Total Reads:\t\t", prettyNum(num_reads, big.mark = ","), "\n", 
        sep = "")
  }
  
  # Average read_length
  avg_read_length <- 
    dataset |> 
    summarise(
      avg_length = mean(read_length, na.rm = TRUE)) |> 
    collect()
  
  if (!silently) {
    cat("  Mean Calls/Read:\t",
        prettyNum(round(num_rows / num_reads)), "\n",
        "  Mean Read Length:\t", 
        prettyNum(round(avg_read_length$avg_length), big.mark = ","), "\n",
        sep = "")
  }
  
  # Unique CpG coverage
  
  unique_cpgs <- dataset |> 
    count(chrom, start) |> 
    collect()
  
  total_cpgs <- nrow(unique_cpgs)
  
  cpg_coverage <- sapply(min_reads, 
                         \(x) 
                         {
                           covered <- sum(unique_cpgs$n >= x)
                           percent <- covered / total_cpgs * 100
                           return(c(covered, percent))
                         })
  
  if (!silently) {
    cat("\nModified Base Coverage:\n")
    for (i in seq_along(min_reads)) {
      cat("  Cov >= ", min_reads[i], ":\t\t", 
          prettyNum(cpg_coverage[1, i], big.mark = ","), 
          "\t(", round(cpg_coverage[2, i], 2), "%)\n", sep = "")
    }
  }
  
  # Calls by flag value
  flag_counts <- 
    dataset |> 
    count(flag) |> 
    mutate(percent = n / sum(n) * 100) |> 
    collect()
  
  if(!silently) {
    cat("\nCalls by Flag Value:\n")
    for (i in seq_len(nrow(flag_counts))) {
      cat("  Flag ", flag_counts$flag[i], 
          ":\t\t", prettyNum(flag_counts$n[i], big.mark = ","), 
          "\t(", round(flag_counts$percent[i], 2), "%)\n", sep = "")
    }
  }
  
  # Calls with call_prob >= threshold
  high_conf_calls <- 
    dataset |> 
    summarise(
      count = sum(call_prob >= call_prob_threshold, na.rm = TRUE),
      percent = count / n() * 100) |> 
    collect()
  
  if (!silently) {
    cat("\nHigh Quality Calls:\n  Mod Prob >= ", call_prob_threshold, ":\t",
        prettyNum(high_conf_calls$count, big.mark = ","), 
        "\t(", round(high_conf_calls$percent, 2), "%)\n", 
        sep = "")
  }
  
  # Calls with base_qual >= threshold
  high_qual_calls <- 
    dataset |> 
    summarise(
      count = sum(base_qual >= base_qual_threshold, na.rm = TRUE),
      percent = count / n() * 100) |> 
    collect()
  
  if (!silently) {
    cat("  Base Qual >= ", base_qual_threshold, ":\t", 
        prettyNum(high_qual_calls$count, big.mark = ","), 
        "\t(", round(high_qual_calls$percent, 2), "%)\n", 
        sep = "")
  }
  
  # Return results invisibly
  invisible(
    list(num_calls = num_rows,
         num_reads = num_reads,
         cpg_coverage = cpg_coverage,
         flag_counts = flag_counts,
         high_conf_calls = high_conf_calls,
         high_qual_calls = high_qual_calls,
         avg_read_length = avg_read_length))
}