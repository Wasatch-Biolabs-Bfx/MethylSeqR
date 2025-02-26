#' Quality Control Wrapper for Methylation Data
#'
#' This function serves as a wrapper for various quality control analyses on methylation data.
#' It sequentially calculates coverage statistics, modification statistics, correlation analysis,
#' and performs principal component analysis (PCA).
#'
#' @param ch3_db A database connection or object containing methylation data.
#' @param call_type A character string indicating the type of call to retrieve data (e.g., "positions", "regions").
#' @param plot A logical value indicating whether to generate plots for the statistical analyses.
#'               Defaults to TRUE. Set to FALSE to suppress plots.
#'
#' @return This function does not return a value. It performs calculations and may produce plots based on the
#'         results of the analysis.
#'
#' @details The function first calculates coverage statistics using `get_cov_stats()`, then computes modification
#'          statistics with `get_mod_stats()`, follows up with correlation analysis via `cor_modseq()`, and finally
#'          runs principal component analysis with `pca_modseq()`.
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Run quality control wrapper
#'  run_qc(ch3_db, call_type = "positions")
#'   
#' @export
run_qc <- function(ch3_db, call_type = "positions", plot = TRUE) {
  message("calculating coverage stats...")
  get_cov_stats(ch3_db, call_type, plot = plot)
  message("calculating mod stats...")
  get_mod_stats(ch3_db, call_type, plot = plot)
  message("calculating correlations...")
  cor_modseq(ch3_db, call_type, plot = plot)
  message("running pca...")
  pca_modseq(ch3_db, call_type)
}
