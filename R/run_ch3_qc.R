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
#' @param max_rows The maximum amount of rows wanted for calculation. This argument can help analysis run faster when there is a lot of data.
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
run_ch3_qc <- function(ch3_db, 
                       call_type = "positions", 
                       plot = TRUE, 
                       max_rows = NULL) {
  
  accepted_tables = c("calls", "positions", "regions", "windows")
  
  if (!(call_type %in% accepted_tables)) {
    stop("Quality control plots can currently only be made for calls, positions, regions, or a windows table in your database. Please select one of these.\n")
  }
  
  cat(paste0("Running quality control on ", call_type, " table."))
  cat("\n")
  
  start_time = Sys.time()
  
  if (call_type == "calls") {
    suppressMessages(suppressWarnings(
      summarize_ch3_positions(ch3_db)
    ))
    call_type = "positions" # switch call_type to positions now
  }
  
  message("calculating coverage stats...")
  plot_ch3_cov(ch3_db, call_type, plot = plot, max_rows = max_rows)
  message("calculating mod stats...")
  plot_ch3_modfrac(ch3_db, call_type, plot = plot, max_rows = max_rows)
  message("calculating correlations...")
  calc_ch3_samplecor(ch3_db, call_type, plot = plot, max_rows = max_rows)
  message("running pca...")
  plot_ch3_pca(ch3_db, call_type, max_rows = max_rows)
  
  end_time = Sys.time()
  
  message("QC complete! Time elapsed: ", end_time - start_time, "\n")
}
