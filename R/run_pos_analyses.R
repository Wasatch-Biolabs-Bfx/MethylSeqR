run_pos_analyses <- function(modseq_dat,
                             cases,
                             controls) {
  # wrapper
  modseq_dat |>
    summarize_by_pos() |>
    calc_mod_diff(
      cases = cases,
      controls = controls,
      calc_type = "fast_fisher") |>
    compute()
}
