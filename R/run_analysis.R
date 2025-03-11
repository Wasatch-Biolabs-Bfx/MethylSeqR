#' Run Differential Analysis on Methylation Data
#'
#' This function performs a differential analysis on methylation data based on the specified call type
#' and applies the appropriate summarization method. It supports analysis of positions, regions, or windows.
#' The function handles the summarization of methylation data and performs differential modification analysis
#' based on case-control comparisons.
#'
#' @param ch3_db A `ch3_db` object representing the DuckDB database containing methylation data.
#'   The database should include necessary tables for the analysis, such as positions, regions, or windows.
#' @param call_type A character string specifying the type of data to analyze. Must be one of:
#'   \code{"positions"}, \code{"regions"}, or \code{"windows"}. This determines the summarization approach to use.
#' @param region_file A character string specifying the path to the region annotation file (required if `call_type` is `"regions"`).
#'   This file should be in a supported format (e.g., BED, CSV, TSV).
#' @param window_size An integer specifying the window size for summarizing methylation data if `call_type` is `"windows"`.
#'   The default value is 1000.
#' @param step_size An integer specifying the step size for sliding windows if `call_type` is `"windows"`.
#'   The default value is 10.
#' @param cases A character vector of sample names to be used as cases in the differential analysis.
#'   This argument is required and cannot be NULL.
#' @param controls A character vector of sample names to be used as controls in the differential analysis.
#'   This argument is required and cannot be NULL.
#' @param mod_type A character string specifying the modification type to analyze. The default is `"mh"`,
#'   which includes both methylation and hydroxymethylation.
#'   Other options are `"c"` for unmodified cytosine, `"m"` for methylation, and `"h"` for hydroxymethylation.
#' @param calc_type A character string specifying the type of statistical test to use for the differential analysis.
#'   The default is `"fast_fisher"`, but other calculation methods can be implemented.
#'
#' @details
#' This function first summarizes the methylation data by the specified call type (positions, regions, or windows).
#' It then proceeds with a differential modification analysis between the provided case and control samples.
#' The analysis is tailored based on the selected modification type (`mod_type`) and calculation method (`calc_type`).
#'
#' @return The result of the differential analysis, typically in the form of a table or data frame with
#'   calculated statistics and p-values for each position, region, or window, depending on the `call_type`.
#'   The result is printed to the console.
#'
#' @import dplyr
#' @importFrom DBI dbExecute
#'
#' @export
run_analysis <- function(ch3_db,
                         call_type,
                         region_file = NULL,
                         window_size = 1000,
                         step_size = 10,
                         cases,
                         controls,
                         mod_type = "mh",
                         calc_type = "fast_fisher") {
  
  if (missing(call_type) || is.null(call_type)) {
    stop("Error: The 'call_type' argument is required and cannot be NULL. Please specify 'positions', 'regions', or 'windows'.")
  }
  
  if (missing(cases) || is.null(cases)) {
    stop("Error: The 'cases' argument is required and cannot be NULL. Please specify which samples are to be cases for differential analysis.")
  }
  
  if (missing(controls) || is.null(controls)) {
    stop("Error: The 'controls' argument is required and cannot be NULL. Please specify which samples are to be controls for differential analysis.")
  }
  
  get_db_stats(ch3_db)
  
  cat("\nRunning Analyses...\n")
  
  # First, summarize by call type requested
  if (call_type == "positions") { 
    summarize_positions(ch3_db)
  } else if (call_type == "regions") {
    if (is.null(region_file)) { # check to make sure user included an annotation file
      stop(paste("Error: region annotation file missing.\n
                 Please include path to annotation file in the region_file argument of this function.\n"))
    }
    summarize_regions(ch3_db, region_file)
  } else if (call_type == "windows") {
    summarize_windows(ch3_db, window_size, step_size)
  }
  
  # Second, run differential modification analysis...
  cat("\nBeginning Differential Analysis...\n")
  calc_mod_diff(ch3_db,
                call_type,
                cases,
                controls,
                mod_type,
                calc_type)
    
  
}