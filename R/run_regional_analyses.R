#' Run Regional Analyses Wrapper on Methylation Data
#'
#' This function serves as a wrapper to perform regional analyses on methylation
#' data by summarizing it by position, aggregating the data into specified regions
#' using an annotation file, and calculating the difference in methylation
#' between designated case and control groups using the fast Fisher's exact test.
#'
#' @param modseq_dat A data frame containing methylation data. It should include relevant
#'                    columns such as sample names, methylation counts, and coverage.
#'
#' @param annot_file A character string representing the path to an annotation file
#'                   that specifies the regions to aggregate the methylation data into.
#'                   The file should contain columns like 'chrom', 'start', 'end',
#'                   and optionally 'region_name'.
#'
#' @param cases A character vector of sample names designated as cases for analysis.
#'
#' @param controls A character vector of sample names designated as controls for analysis.
#'
#' @return A data frame containing the results of the regional analyses, including
#'         differences in methylation between cases and controls, as well as p-values.
#'
#' @examples
#' results <- run_regional_analyses(modseq_data, "regions.tsv", cases = c("case1", "case2"), controls = c("control1", "control2"))
#'
#' @importFrom dplyr summarize
#'
#' @importFrom magrittr %>%
#'
#' @export
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
