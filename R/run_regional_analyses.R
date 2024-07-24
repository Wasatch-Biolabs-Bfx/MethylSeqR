run_regional_analyses <- function(modseq_dat,
                             annot_file,
                             cases,
                             controls) {
  # wrapper function
  modseq_dat |>
    summarize_by_pos() |>
    aggregate_regions(
      annot_file = annot_file) |>
    calc_mod_diff(
      cases = cases,
      controls = controls,
      calc_type = "fast_fisher") |>
    compute()
}
