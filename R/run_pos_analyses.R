#' Run Positional Analyses on Methylation Data
#'
#' This function serves as a wrapper to perform positional analyses on methylation
#' data by summarizing it by position and calculating the difference in methylation
#' between specified case and control groups using the fast Fisher's exact test.
#'
#' @param modseq_dat A data frame containing methylation data. It should include relevant
#'                    columns such as sample names, methylation counts, and coverage.
#'
#' @param cases A character vector of sample names designated as cases for analysis.
#'
#' @param controls A character vector of sample names designated as controls for analysis.
#'
#' @return A data frame containing the results of the positional analyses, including
#'         differences in methylation between cases and controls, as well as p-values.
#'
#' @examples
#' results <- run_pos_analyses(modseq_data, cases = c("case1", "case2"), controls = c("control1", "control2"))
#'
#' @importFrom dplyr summarize
#'
#' @importFrom magrittr %>%
#'
#' @export
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
