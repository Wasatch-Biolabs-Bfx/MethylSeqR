calc_mod_diff <- function(modseq_dat, 
                          cases, 
                          controls, 
                          mod_type = "mh",
                          calc_type = "fast_fisher") 
{
  # Set stat to use
  mod_counts_col <- paste0(mod_type[1], "_counts")
  
  # Label cases and controls
  in_dat <- 
    modseq_dat |>
    select(
      sample_name, any_of(c("chrom", "ref_position", "region_name")), 
      cov, c_counts, mod_counts = !!mod_counts_col) |>
    mutate(
      exp_group = case_when(
        sample_name %in% cases ~ "case",
        sample_name %in% controls ~ "control",
        TRUE ~ NA)) |>
    filter(
      !is.na(exp_group))
  
  # Calculate p-vals and diffs
  switch(calc_type,
         fast_fisher = .calc_diff_fisher(in_dat,
                                         calc_type = "fast_fisher"),
         r_fisher    = .calc_diff_fisher(in_dat,
                                         calc_type = "r_fisher"),
         log_reg     = .calc_diff_logreg(in_dat)) |>
    rename_with(
      ~ gsub("mod", mod_type[1], .x))
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
      .by = c(exp_group, any_of(c("chrom", "ref_position", "region_name"))),
      c_counts = sum(c_counts),
      mod_counts = sum(mod_counts)) |>
    pivot_wider(
      names_from = exp_group, 
      values_from = c(c_counts, mod_counts),
      values_fill = 0) 
  
  # Extract matrix and calculate p-vals
  pvals <- 
    dat |>
    select(!any_of(c("chrom", "ref_position", "region_name"))) |>
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
      .by = c(chrom, ref_position),
      mean_cov = mean(cov),
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
  
  