#' Summarize methylation data by reference positions
#'
#' @param modseq_dat Methylation data in a duckdb object as opened by open_dat().
#' @param score A methylation score type. Default is mh, can also include m, h, or a vector of score types wanted.
#' @return A duckdb object of methylation scores summarized by each position.
#' @examples
#' summarize_by_pos(data)
#' summarize_by_pos(data, score = "m")
#' summarize_by_pos(data, score = c("m", "mh", "h")
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
      sample_name, chrom, ref_position, cov, c_counts, !!!score_cols)
}
