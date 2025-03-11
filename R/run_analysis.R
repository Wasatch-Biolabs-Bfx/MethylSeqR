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