#' Calculate Differential Methylation
#'
#' This function calculates differential methylation between specified case and control groups using various statistical methods. 
#' The results are stored in a DuckDB database for further analysis.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param call_type A string representing the name of the table in the database from which to pull the data. 
#' Default is "positions".
#' @param cases A character vector containing the sample names for the case group.
#' @param controls A character vector containing the sample names for the control group.
#' @param mod_type A string indicating the type of modification to analyze. 
#' Default is "mh" for methylation/hydroxymethylation.
#' @param calc_type A string specifying the statistical method to use for calculating p-values. 
#' Options include "fast_fisher", "r_fisher", and "log_reg". Default is "fast_fisher".
#' @param collapse_windows A boolean specifiying whether the user would like to collapse windows and create a more defined windows of signficance.
#'
#' @details
#' The function connects to the specified DuckDB database and retrieves methylation data from the specified call type table. 
#' It summarizes the data for cases and controls, calculates p-values based on the specified method, and stores the results in the 
#' "meth_diff" table. 
#'
#' @return A list containing the updated "ch3_db" object with the latest tables in the database, including "meth_diff".
#'
#' @import dbplyr
#' @import duckdb
#' 
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Get methylation statistics for the 'positions' call type without plotting
#'  calc_mod_diff(ch3_db = ch3_db, 
#'                call_type = "positions",
#'                cases = c("Blood1_chr21", "Blood2_chr21", "Blood3_chr21"),
#'                controls = c("Sperm1_chr21", "Sperm2_chr21", "Sperm3_chr21")))
#'
#' @export
calc_mod_diff <- function(ch3_db,
                          call_type = "positions",
                          cases,
                          controls,
                          mod_type = "mh",
                          calc_type = "fast_fisher",
                          collapse_windows = FALSE)
{
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # Finish up: purge extra tables & update table list and close the connection
  keep_tables = c("calls","positions", "regions", "windows", "mod_diff")
  on.exit(.helper_purgeTables(db_con, keep_tables), add = TRUE)
  on.exit(dbExecute(db_con, "VACUUM;"), add = TRUE)  # <-- Ensure space is reclaimed
  on.exit(.helper_closeDB(database), add = TRUE)

  # check for windows function
  if (!dbExistsTable(db_con, call_type)) { # add db_con into object and put in every function...
    stop(paste0(call_type, " table does not exist. Build it with summarize_positions, summarize_regions, or summarize_windows."))
  }
  
  if (dbExistsTable(db_con, "mod_diff"))
    dbRemoveTable(db_con, "mod_diff")
  dbExecute(db_con, "VACUUM;")  # <-- Add this to free space immediately
  
  # Set stat to use
  mod_counts_col <- paste0(mod_type[1], "_counts")
  
  # Label cases and controls
  in_dat <-
    tbl(db_con, call_type) |>
    select(
      sample_name, any_of(c("region_name", "chrom", "start", "end")),
      c_counts, mod_counts = !!mod_counts_col) |>
    mutate(
      exp_group = case_when(
        sample_name %in% cases ~ "case",
        sample_name %in% controls ~ "control",
        TRUE ~ NA)) |>
    filter(
      !is.na(exp_group))
  
  # Calculate p-vals and diffs
  result <- switch(calc_type,
                   fast_fisher = .calc_diff_fisher(in_dat,
                                                   calc_type = "fast_fisher"),
                   r_fisher    = .calc_diff_fisher(in_dat,
                                                   calc_type = "r_fisher"),
                   log_reg     = .calc_diff_logreg(in_dat)) |>
    rename_with(
      ~ gsub("mod", mod_type[1], .x))
  
  # working on...
  # if (call_type == "windows" && collapse_windows == TRUE) {
  #   cat("Collapsing windows...\n\n")
  #   result <- .collapse_windows(result, calc_type) 
  #   
  #   result |>
  #     mutate(p_adjust = p.adjust(p_val, method = "BH")) |>
  #     arrange(p_adjust) |>
  #     dbWriteTable(
  #       conn = db_con,
  #       name = "mod_diff",
  #       append = TRUE
  #     )
  # } else {
    # Collect the result and write to the database
  result |>
    collect() |>
    mutate(p_adjust = p.adjust(p_val, method = "BH")) |>
    arrange(p_adjust) |>
    dbWriteTable(
      conn = db_con, 
      name = "mod_diff", 
      append = TRUE
    )

  message("Mod diff analysis complete - mod_diff table successfully created!")
  print(head(tbl(db_con, "mod_diff")))

  invisible(database)
}

.collapse_windows <- function(result, calc_type)
{
  # filter first for significance, unchanged p-value,  whats a significant P-value?
  # filter for all + (hypermethylated), collapse windows, then filter for - (hypomethylated) thn collapse
  # a user can choose the difference tthreshold  
  # throw it back inot a regional analysis
  
  # Collect the result
  result <- result |> 
    collect() |> 
    arrange(chrom, start)
  
  print(head(result))
  
  
  # Initialize lists to store results
  collapsed <- list()
  current_window <- NULL
  
  result$type <- ifelse(result$meth_diff > 0, "hypermethylated", 
                     ifelse(result$meth_diff < 0, "hypomethylated", "unchanged"))
  
  result <- result |> drop_na(meth_diff)
  
  # Initialize the progress bar
  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent ETA: :eta",
    total = nrow(result), clear = FALSE, width = 60
  )
  
  for (i in 1:nrow(result)) {
    row <- result[i,]
    
    if (is.null(current_window)) {
      current_window <- row
      next
    }
    
    # Check if current row overlaps with current window and has same toi_meth
    if (current_window$chrom == row$chrom &&
        current_window$type == row$type &&
        row$start <= current_window$end) {
      
      # Merge windows by updating the end position and aggregating metrics
      current_window$end <- max(current_window$end, row$end)
      current_window$c_counts_control <- current_window$c_counts_control + row$c_counts_control
      current_window$c_counts_case <- current_window$c_counts_case + row$c_counts_case
      current_window$mh_counts_control <- current_window$mh_counts_control + row$mh_counts_control
      current_window$mh_counts_case <- current_window$mh_counts_case + row$mh_counts_case
      
      # Recalculate metrics
      # current_window$cov_control <- current_window$c_counts_control + current_window$mh_counts_control
      # current_window$cov_case <- current_window$c_counts_case + current_window$mh_counts_case
      current_window$mh_frac_case <- current_window$mh_counts_case / (current_window$c_counts_case + current_window$mh_counts_case) # divided by coverage of cases
      current_window$mh_frac_control <- current_window$mh_counts_control / (current_window$c_counts_control + current_window$mh_counts_control) # divided by coverage of controls
      current_window$meth_diff <- current_window$mh_frac_case - current_window$mh_frac_control
      #current_window$abs_meth_diff <- abs(current_window$meth_diff)
      
      # YOU'LL NEED TO CHANGE THIS PART TO GET THE CORRECT NEW PVALUE!!! This just takes the minimum p-value and adjusted p-value 
      current_window$p_val <- min(current_window$p_val, row$p_val)
      
    } else {
      # Store current window and start a new one
      collapsed[[length(collapsed) + 1]] <- current_window
      current_window <- row
    }
    
    # Update the progress bar
    pb$tick()
  }
  
  cat("\n")
  
  # Don't forget to add the last window
  if (!is.null(current_window)) {
    collapsed[[length(collapsed) + 1]] <- current_window
  }
  
  # Convert list to data frame
  result <- bind_rows(collapsed) |>
    mutate(dmr_length = end - start)
  
  return(result)
}

## Calculate p-values using fisher exact tests. If there are multiple samples,
## they will be combined.
.calc_diff_fisher <- function(in_dat,
                              calc_type)
{
  # Combine replicates and pivot wider
  dat <-
    in_dat |>
    summarize(
      .by = c(exp_group, any_of(c("region_name", "chrom", "start", "end"))),
      c_counts = sum(c_counts, na.rm = TRUE),
      mod_counts = sum(mod_counts, na.rm = TRUE)) |>
    pivot_wider(
      names_from = exp_group,
      values_from = c(c_counts, mod_counts),
      values_fill = 0)
  
  # Extract matrix and calculate p-vals
  pvals <-
    dat |>
    select(!any_of(c("region_name", "chrom", "start", "end"))) |>
    distinct() |>
    mutate(
      mod_frac_case = mod_counts_case /
        (mod_counts_case + c_counts_case),
      mod_frac_control = mod_counts_control /
        (mod_counts_control + c_counts_control),
      meth_diff = mod_counts_case /
        (mod_counts_case + c_counts_case) -
        mod_counts_control /
        (mod_counts_control + c_counts_control)) |>
    collect() |>
    mutate(
      p_val = switch(
        calc_type,
        fast_fisher = .fast_fisher(
          q = mod_counts_case,
          m = mod_counts_case + mod_counts_control,
          n = c_counts_case + c_counts_control,
          k = mod_counts_case + c_counts_case),
        r_fisher = .r_fisher(
          a = mod_counts_control,
          b = mod_counts_case,
          c = c_counts_control,
          d = c_counts_case)))
  
  dat |>
    inner_join(
      pvals,
      by = join_by(c_counts_case, c_counts_control,
                   mod_counts_case, mod_counts_control),
      copy = TRUE)
  
  
}


.fast_fisher <- function(q, m, n, k)
{
  # derived from https://github.com/al2na/methylKit/issues/96
  
  mat <- cbind(q, m, n, k)
  
  apply(mat, 1,
        \(qmnk)
        {
          dhyper_val <- 0.5 * dhyper(x = qmnk[1], m = qmnk[2],
                                     n = qmnk[3], k = qmnk[4])
          
          pval_right <- phyper(q = qmnk[1], m = qmnk[2],
                               n = qmnk[3], k = qmnk[4],
                               lower.tail = FALSE) + dhyper_val
          
          pval_left  <- phyper(q = qmnk[1] - 1, m = qmnk[2],
                               n = qmnk[3], k = qmnk[4],
                               lower.tail = TRUE) + dhyper_val
          
          return(ifelse(test = pval_right > pval_left,
                        yes  = pval_left * 2,
                        no   = pval_right * 2))
        })
}


.r_fisher <- function(a, b, c, d)
{
  mat <- cbind(a, b, c, d)
  
  apply(mat, 1,
        \(x)
        {
          fisher.test(matrix(x, 2))$p.val
        })
}


.calc_diff_logreg <- function(in_dat)
{
  pvals <-
    in_dat |>
    mutate(
      cov = c_counts + mod_counts,
      mod_frac = mod_counts / cov) |>
    collect() |>
    summarize(
      .by = c(chrom, start),
      mean_cov = mean(cov, na.rm = TRUE),
      mean_frac_case = mean(mod_frac[exp_group == "case"]),
      mean_frac_ctrl = mean(mod_frac[exp_group == "control"]),
      mean_diff = mean_frac_case - mean_frac_ctrl,
      p_val = .logreg(mod_frac, cov, exp_group))
  
  # Pivot wider, add pvals, and return
  in_dat |>
    pivot_wider(
      id_cols = c(chrom, ref_position),
      names_from = sample_name,
      values_from = c(c_counts, mod_counts),
      values_fill = 0) |>
    inner_join(
      pvals, by = join_by(chrom, ref_position),
      copy = TRUE)
}


.logreg <- function(mod_frac,
                    cov,
                    exp_group)
{
  exp_group <- as.numeric(factor(exp_group))
  fit <- glm.fit(exp_group, mod_frac,
                 weights = cov / sum(cov), family = binomial())
  deviance <- fit$null.deviance - fit$deviance
  
  pchisq(deviance, 1, lower.tail = FALSE)
}