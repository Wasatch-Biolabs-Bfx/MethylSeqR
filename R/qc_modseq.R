#' Quality Control Wrapper for Methylation Data
#'
#' This function performs a series of quality control analyses on methylation data,
#' including coverage statistics, methylation statistics, correlation analysis,
#' clustering, PCA, and methylation density estimation.
#'
#' @param modseq_dat A data frame containing methylation data. It should include
#'                    relevant columns such as sample names, coverage, and methylation fractions.
#'
#' @param plot A logical value indicating whether to plot the results. Defaults to TRUE.
#'
#' @return The function performs various quality control analyses and generates plots
#'         as specified by the `plot` parameter. It does not return any values.
#'
#' @examples
#' qc_modseq(data)
#'
#' @importFrom dplyr print
#'
#' @importFrom magrittr %>%
#'
#' @export
qc_modseq <- function(modseq_dat, plot = TRUE) {
  print("calculating coverage stats...")
  get_cov_stats(modseq_dat = modseq_dat, plot = plot)
  print("calculating mod stats...")
  get_mod_stats(modseq_dat = modseq_dat, plot = plot)
  print("calculating correlations...")
  cor_modseq(modseq_dat = modseq_dat, plot = plot)
  print("calculating cluster...")
  cluster_modseq(modseq_dat = modseq_dat)
  print("running pca...")
  pca_modseq(modseq_dat = modseq_dat)
  print("calculating methylation density...")
  density_modseq(modseq_dat = modseq_dat)
}
