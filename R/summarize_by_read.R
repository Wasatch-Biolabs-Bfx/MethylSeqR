
summarize_by_read <- function(modseq_dat,
                              score = "mh",
                              include_ID = FALSE)
{ 
  score_cols <- 
    expand.grid(score, c("counts", "frac")) |>
    unite("x") |>
    unlist() |>
    unname()

# determine which columns to be collected later once we have summarized by read
if (include_ID == TRUE) {
  selected_columns = c("sample_name", "chrom", "first_CG_pos", 
                       "last_CG_pos", "read_length", score_cols, "read_id")
} else if (include_ID == FALSE) {
  selected_columns = c("sample_name", "chrom", "first_CG_pos", 
                       "last_CG_pos", "read_length", score_cols)
}

modseq_dat |>
  summarize(
    .by = c(sample_name, chrom, read_id, read_length),
    cov = n(),
    c_counts  = sum(as.integer(call_code == "-"),
                    na.rm = TRUE),
    m_counts  = sum(as.integer(call_code == "m"),
                    na.rm = TRUE),
    h_counts  = sum(as.integer(call_code == "h"),
                    na.rm = TRUE),
    mh_counts = sum(as.integer(call_code %in% c("m", "h")),
                    na.rm = TRUE),
    first_CG_pos = min(ref_position),
    last_CG_pos = max(ref_position)) |>
  mutate(
    m_frac = m_counts / cov,
    h_frac = h_counts / cov,
    mh_frac = mh_counts / cov,
    na.rm = TRUE) |>
  select(any_of(selected_columns))
}