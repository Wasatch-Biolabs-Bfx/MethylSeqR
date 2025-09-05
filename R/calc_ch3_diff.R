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
#' @importFrom dplyr tbl select any_of mutate case_when filter pull summarize inner_join join_by rename_with collect arrange
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
    stop(call_type, " table does not exist. Build it with summarize_ch3_positions(), ",
         "summarize_ch3_regions(), or summarize_ch3_windows().")
  }
  
  # Discover available *_counts columns and validate mod_type
  cols <- colnames(dplyr::tbl(ch3_db$con, call_type))
  counts_cols <- grep("_counts$", cols, value = TRUE)
  available_labels <- sub("_counts$", "", counts_cols)
  
  mod_type <- mod_type[1]
  mod_counts_col <- paste0(mod_type, "_counts")
  if (!mod_counts_col %in% cols) {
    stop("Requested mod_type '", mod_type, "' not found in ", call_type, ". ",
         "Available mod types are: ",
         paste(sort(available_labels), collapse = ", "), ".\n",
         "Tip: rebuild ", call_type, " with summarize_* and include mod_code = '", mod_type, "'.")
  }
  if (!"num_calls" %in% cols) {
    stop(call_type, " is missing required column 'num_calls'. ",
         "Recreate it with the latest summarize_* function.")
  }
  
  mod_diff_table <- paste0("mod_diff_", call_type)
  if (DBI::dbExistsTable(ch3_db$con, mod_diff_table)) {
    DBI::dbRemoveTable(ch3_db$con, mod_diff_table)
  }
  
  cat("Running differential analysis...\n")
  
  in_dat <-
    dplyr::tbl(ch3_db$con, call_type) |>
    dplyr::select(
      sample_name,
      dplyr::any_of(c("region_name", "chrom", "start", "end")),
      num_calls,
      mod_counts = dplyr::all_of(mod_counts_col)
    ) |>
    dplyr::mutate(
      exp_group = dplyr::case_when(
        sample_name %in% cases ~ "case",
        sample_name %in% controls ~ "control",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(exp_group))
  
  # Check sample names present
  all_samples <- unique(dplyr::pull(in_dat, sample_name))
  if (any(!cases %in% all_samples)) {
    missing <- paste(cases[!cases %in% all_samples], collapse = ", ")
    stop("Check case names - some case samples are missing from the data: ", missing)
  }
  if (any(!controls %in% all_samples)) {
    missing <- paste(controls[!controls %in% all_samples], collapse = ", ")
    stop("Check control names - some control samples are not in the data: ", missing)
  }
  
  # Compute p-values / diffs
  result <- switch(calc_type,
                   fast_fisher = .calc_diff_fisher(in_dat, calc_type = "fast_fisher"),
                   r_fisher    = .calc_diff_fisher(in_dat, calc_type = "r_fisher"),
                   log_reg     = .calc_diff_logreg(in_dat),
                   stop("Unknown calc_type: ", calc_type, ". Use 'fast_fisher', 'r_fisher', or 'log_reg'.")
  ) |>
    # rename mod_* columns to e.g. a_* or mh_* to reflect the chosen mod_type
    dplyr::rename_with(~ gsub("^mod", mod_type, .x))
  
  
  result |>
    dplyr::collect() |>
    dplyr::mutate(p_adjust = stats::p.adjust(p_val, method = "BH")) |>
    dplyr::arrange(p_adjust) |>
    DBI::dbWriteTable(conn = ch3_db$con, name = mod_diff_table, append = TRUE)

  end_time <- Sys.time()
  total_time_difftime <- end_time - start_time
  
  # Convert the total_time_difftime object to numeric seconds for a reliable comparison
  total_seconds <- as.numeric(total_time_difftime, units = "secs")
  
  if (total_seconds > 60) {
    # If greater than 60 seconds, convert to numeric minutes for display
    total_minutes <- as.numeric(total_time_difftime, units = "mins")
    message("Mod diff analysis complete! ", mod_diff_table, " table successfully created!",
            "\nTime elapsed: ", round(total_minutes, 2), " minutes\n")
  } else {
    # Otherwise, display in numeric seconds
    message("Mod diff analysis complete! ", mod_diff_table, " table successfully created!", 
            "\nTime elapsed: ", round(total_seconds, 2), " seconds\n")
  }
  
  if (call_type == "windows") {
    message("Call collapse_ch3_windows() to collapse significant windows.\n")
  } 
  
  # Print a preview of what table looks like
  print(head(dplyr::tbl(ch3_db$con, mod_diff_table)))
  
  ch3_db$current_table = mod_diff_table
  ch3_db <- .ch3helper_cleanup(ch3_db)
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
    dplyr::summarize(
      .by = c(exp_group, any_of(c("region_name", "chrom", "start", "end"))),
      num_calls = sum(num_calls, na.rm = TRUE),
      mod_counts = sum(mod_counts, na.rm = TRUE)) |>
    dplyr::mutate(
      c_counts = num_calls - mod_counts) |>
    pivot_wider(
      names_from = exp_group,
      values_from = c(num_calls, mod_counts, c_counts),
      values_fill = 0)
  
  # Extract matrix and calculate p-vals
  pvals <-
    dat |>
    dplyr::select(!any_of(c("region_name", "chrom", "start", "end"))) |>
    distinct() |>
    dplyr::mutate(
      mod_frac_case = mod_counts_case /
        (num_calls_case),
      mod_frac_control = mod_counts_control /
        (num_calls_control),
      meth_diff = mod_counts_case /
        (num_calls_case) -
        mod_counts_control /
        (num_calls_control)) |>
    collect() |>
    dplyr::mutate(
      p_val = switch(
        calc_type,
        fast_fisher = .fast_fisher(
          q = mod_counts_case,
          m = mod_counts_case + mod_counts_control,
          n = c_counts_case + c_counts_control,
          k = num_calls_case),
        r_fisher = .r_fisher(
          a = mod_counts_control,
          b = mod_counts_case,
          c = c_counts_control,
          d = c_counts_case)))
  
  dat |>
    inner_join(
      pvals,
      by = join_by(num_calls_case, num_calls_control,
                   mod_counts_case, mod_counts_control,
                   c_counts_case,   c_counts_control),
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
      cov = num_calls + mod_counts,
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
      values_from = c(num_calls, mod_counts),
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