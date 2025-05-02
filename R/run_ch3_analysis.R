#' Run Differential Analysis on Methylation Data
#'
#' This function performs a differential analysis on methylation data based on the specified call type
#' and applies the appropriate summarization method. It supports analysis of positions, regions, or windows.
#' The function handles the summarization of methylation data and performs differential modification analysis
#' based on case-control comparisons.
#'
#' @param ch3_db A `ch3_db` object representing the DuckDB database containing methylation data.
#'   The database should include necessary tables for the analysis, such as positions, regions, or windows.
#' @param out_path The directory in which the "Mod_Diff_Analysis_Results" directory containing result data will be written out too. 
#' If the user does not provide a directory, the working directory will be used.
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
#' @param p_val_max The p value threshold in which significant differentially modified positions/regions/windows 
#' will be written out in the final directory.
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
#' @importFrom dplyr select filter matches pivot_longer pivot_wider unite semi_join
#' @importFrom utils write.csv
#' @importFrom stats setNames
#' @importFrom fs dir_create
#'
#' @export

run_ch3_analysis <- function(ch3_db,
                         out_path,
                         call_type,
                         region_file = NULL,
                         window_size = 1000,
                         step_size = 10,
                         cases,
                         controls,
                         mod_type = "mh",
                         calc_type = "fast_fisher",
                         p_val_max = 0.05) {
  
  if (missing(call_type) || is.null(call_type)) {
    stop("Error: The 'call_type' argument is required and cannot be NULL. Please specify 'positions', 'regions', or 'windows'.")
  }
  
  if (missing(cases) || is.null(cases)) {
    stop("Error: The 'cases' argument is required and cannot be NULL. Please specify which samples are to be cases for differential analysis.")
  }
  
  if (missing(controls) || is.null(controls)) {
    stop("Error: The 'controls' argument is required and cannot be NULL. Please specify which samples are to be controls for differential analysis.")
  }
  
  # get_db_stats(ch3_db)
  
  cat("\nRunning Analyses...\n")
  
  # Record the start time
  start_time <- Sys.time()
  
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
  
  # Third, create directory for user with meth_diff table, all methylation data, and the significant regions/windows/posiitons meth data...
  if (is.null(out_path)) {
    out_path = getwd() 
  }
  
  cat(paste0("\nWriting out differential analysis results to ", out_path, "\n"))
  
  new_dir <- file.path(out_path, "Mod_Diff_Analysis_Results")
  
  if (!dir.exists(new_dir)) {
    dir.create(new_dir, recursive = TRUE)
  }
  
  # Write out the mod_diff table from analysis
  mod_diff_path <- file.path(new_dir, "mod_diff.csv")
  # Write the data frame to a CSV file
  mod_diff = get_table(ch3_db, "mod_diff")
  write.csv(mod_diff, mod_diff_path, row.names = FALSE)
  
  # Write out all methylation data form all samples
  all_CGs_path <- file.path(new_dir, "All_CpGs.csv")
  data = get_table(ch3_db, call_type)
  
  cat("\nWriting out all CpG data...\n")
  df_wide <- data |>
    select(sample_name, chrom, start, end, matches("_frac")) |>
    pivot_longer(cols = matches("_frac"),  # pivot both mh_frac and h_frac
                 names_to = "frac_type",     # create a new column 'frac_type'
                 values_to = "value") |>    # pivot values into 'value' column
    unite("sample_frac", sample_name, frac_type, sep = ".") |>  # combine sample and frac_type
    pivot_wider(names_from = sample_frac, values_from = value) |> # spread the new 'sample_frac' into columns
    select(chrom, start, end, everything()) |>  # Ensure chrom, start, end come first
    select(chrom, start, end, sort(names(.)[4:length(.)]))  # Sort the sample columns in alphabetical order
  
  write.csv(df_wide, all_CGs_path, row.names = FALSE)
  
  # Write out only the significant regions of interests methylation data...
  cat("\nWriting out significant CpG data based on differential analysis...\n")
  cat(paste0("p-value threshold: ", p_val_max, "\n"))
  Sig_CGs_path <- file.path(new_dir, "Sig_CpGs.csv")
  
  sig = mod_diff |>
    filter(p_val <= p_val_max) |>
    select("chrom", "start", "end")
  
  # Now, filter the data to only keep rows that match chrom, start, end in sig
  df_wide_significant <- df_wide |>
    semi_join(sig, by = c("chrom", "start", "end"))
  
  write.csv(df_wide_significant, Sig_CGs_path, row.names = FALSE)
  
  cat("\nRun Analysis Complete!\n")
  
  # Record the end time
  end_time <- Sys.time()
  
  # Calculate the elapsed time
  elapsed_time <- end_time - start_time
  
  # Print out the elapsed time
  cat("Time elapsed:", elapsed_time, "\n")
  
}