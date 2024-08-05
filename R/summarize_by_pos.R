#' Summarize Methylation Data by Reference Positions
#'
#' This function summarizes methylation data by reference positions, calculating
#' coverage and counts of different methylation states for each sample. It can handle
#' multiple methylation score types as specified by the user.
#'
#' @param modseq_dat A duckdb object containing methylation data, typically opened
#'                   using the `open_dat()` function.
#'
#' @param score A character vector specifying the methylation score types to summarize.
#'              Default is "mh", but can also include "m", "h", or any combination of
#'              these scores.
#'
#' @return A duckdb object containing methylation scores summarized by each reference position,
#'         including counts and fractions of methylated and unmethylated bases for each sample.
#'
#' @examples
#' summarized_data <- summarize_by_pos(data)
#' summarized_data_m <- summarize_by_pos(data, score = "m")
#' summarized_data_multiple <- summarize_by_pos(data, score = c("m", "mh", "h"))
#'
#' @export
summarize_by_pos <- function(modseq_dat,
                             score = "mh")
{
  score_cols <- expand.grid(score, c("counts", "frac")) |>
    unite("x") |>
    unlist() |>
    unname()

  modseq_dat |>
    summarize(
      .by = c(sample_name, chrom, ref_position),
      cov = n(),
      c_counts = sum(as.integer(call_code == "-"),
                     na.rm = TRUE),
      m_counts = sum(as.integer(call_code == "m"),
                     na.rm = TRUE),
      h_counts = sum(as.integer(call_code == "h"),
                     na.rm = TRUE),
      mh_counts = sum(as.integer(call_code %in% c("m", "h")),
                      na.rm = TRUE)) |>
    mutate(
      m_frac = m_counts / cov,
      h_frac = h_counts / cov,
      mh_frac = mh_counts / cov) |>
    select(
      sample_name, chrom, ref_position, cov, c_counts, !!!score_cols) %>%
    collect()
}
