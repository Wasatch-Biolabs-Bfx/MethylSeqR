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
