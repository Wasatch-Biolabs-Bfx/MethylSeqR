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
#'
#' @details
#' The function connects to the specified DuckDB database and retrieves methylation data from the specified call type table. 
#' It summarizes the data for cases and controls, calculates p-values based on the specified method, and stores the results in the 
#' "meth_diff" table. 
#'
#' @return A list containing the updated "ch3_db" object with the latest tables in the database, including "meth_diff".
#' 
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Get methylation statistics for the 'positions' call type without plotting
#'  calc_ch3_diff(ch3_db = ch3_db, 
#'                call_type = "positions",
#'                cases = c("Blood1_chr21", "Blood2_chr21", "Blood3_chr21"),
#'                controls = c("Sperm1_chr21", "Sperm2_chr21", "Sperm3_chr21")))
#'                
#' @importFrom DBI dbConnect dbDisconnect dbExistsTable dbRemoveTable dbExecute dbWriteTable
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl select any_of mutate case_when filter pull unique summarize pivot_wider inner_join rename_with collect arrange
#' @importFrom tidyr pivot_wider
#' @importFrom stats fisher.test p.adjust dhyper phyper glm.fit pchisq
#'
#' @export

calc_ch3_diff <- function(ch3_db,
                          call_type = "positions",
                          cases,
                          controls,
                          mod_type = "mh",
                          calc_type = "fast_fisher")
{
  start_time <- Sys.time()
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)

  # check for windows function
  if (!dbExistsTable(ch3_db$con, call_type)) { # add db_con into object and put in every function...
    stop(paste0(call_type, " table does not exist. Build it with summarize_positions, summarize_regions, or summarize_windows."))
  }
  
  mod_diff_table <- paste0("mod_diff_", call_type)
  
  if (dbExistsTable(ch3_db$con, mod_diff_table)) {
    dbRemoveTable(ch3_db$con, mod_diff_table)
  }
  
  cat("Running differential analysis...")
  cat("\n")
  # Set stat to use
  mod_counts_col <- paste0(mod_type[1], "_counts")
  
  # Label cases and controls
  in_dat <-
    tbl(ch3_db$con, call_type) |>
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
  
  # Get unique sample names in the data
  all_samples <- unique(pull(in_dat, sample_name))
  if (any(!cases %in% all_samples)) {
    stop(paste(
      "Check case names - some case samples are missing from the data:",
      paste(cases[!cases %in% all_samples], collapse = ", ")
    ))
  }
  
  if (any(!controls %in% all_samples)) {
    stop(paste(
      "Check control names- some control samples are not in the data:",
      paste(controls[!controls %in% all_samples], collapse = ", ")
    ))
  }
  
  # Calculate p-vals and diffs
  result <- switch(calc_type,
                   fast_fisher = .calc_diff_fisher(in_dat,
                                                   calc_type = "fast_fisher"),
                   r_fisher    = .calc_diff_fisher(in_dat,
                                                   calc_type = "r_fisher"),
                   log_reg     = .calc_diff_logreg(in_dat)) |>
    rename_with(
      ~ gsub("mod", mod_type[1], .x))
  
  
  result |>
    collect() |>
    mutate(p_adjust = p.adjust(p_val, method = "BH")) |>
    arrange(p_adjust) |>
    dbWriteTable(
      conn = ch3_db$con, 
      name = mod_diff_table, 
      append = TRUE
    )

  end_time <- Sys.time()
  cat("\n")
  if (call_type == "windows") {
    message(paste0("Mod diff analysis complete! ", 
                   mod_diff_table, 
                   " table successfully created!\nTime elapsed: ", 
                   end_time - start_time, 
                   "\n"))
    message("Call collapse_ch3_windows() to collapse significant windows.")
  } else {
    message(paste0("Mod diff analysis complete! ", 
                   mod_diff_table, 
                   " table successfully created!\nTime elapsed: ", 
                   end_time - start_time, 
                   "\n"))
    
    # Print a preview of what table looks like
    print(head(tbl(ch3_db$con, mod_diff_table)))
    
    ch3_db$current_table = mod_diff_table
    ch3_db <- .ch3helper_cleanup(ch3_db)
    invisible(ch3_db)
  }
  
  # if (call_type == "windows" && collapse_windows == TRUE) {
  #   cat("\nCollapsing Windows...\n")
  #   .collapse_windows(db_con)
  #   message("collapsed_windows table successfully created!\n")
  # }
  
  print(head(tbl(ch3_db$con, mod_diff_table)))

  invisible(ch3_db)
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

.fast_fisher <- function(q, m, n, k) {
  # Calculate values once
  dhyper_val <- 0.5 * dhyper(q, m, n, k)
  
  pval_right <- phyper(q, m, n, k, lower.tail = FALSE) + dhyper_val
  pval_left  <- phyper(q - 1, m, n, k, lower.tail = TRUE) + dhyper_val
  
  # Return min tail * 2
  pmin(pval_right, pval_left) * 2
}

# old fast fisher
# .fast_fisher <- function(q, m, n, k)
# {
#   # derived from https://github.com/al2na/methylKit/issues/96
#   
#   mat <- cbind(q, m, n, k)
#   
#   apply(mat, 1,
#         \(qmnk)
#         {
#           dhyper_val <- 0.5 * dhyper(x = qmnk[1], m = qmnk[2],
#                                      n = qmnk[3], k = qmnk[4])
#           
#           pval_right <- phyper(q = qmnk[1], m = qmnk[2],
#                                n = qmnk[3], k = qmnk[4],
#                                lower.tail = FALSE) + dhyper_val
#           
#           pval_left  <- phyper(q = qmnk[1] - 1, m = qmnk[2],
#                                n = qmnk[3], k = qmnk[4],
#                                lower.tail = TRUE) + dhyper_val
#           
#           return(ifelse(test = pval_right > pval_left,
#                         yes  = pval_left * 2,
#                         no   = pval_right * 2))
#         })
# }


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