summarize_by_region <- function(modseq_dat, annotation)
{
  # Make sure annotation file is formatted correctly
  if (ncol(annotation) < 3 || ncol(annotation) > 4) {
    stop("Annotation must be in correct format. ex. chrom, start, end, region_name OR chrom, start, end")
  }

  # Create regional dataframe
  print("stuck 1")
  modseq_dat <- right_join(modseq_dat, annotation, copy = TRUE,
                           by = join_by(chrom, between(ref_position, start, end))) |>
    summarize(mean_cov = mean(cov, na.rm = TRUE))
  
  print("stuck 2")
  
  if ("mh_frac" %in% colnames(modseq_dat)) {
    modseq_dat = modseq_dat |>
      summarize(mean_mh_frac = sum(cov * mh_frac) / sum(cov))
  }
  
  print("stuck 3")
  
  if ("m_frac" %in% colnames(modseq_dat)) {
    modseq_dat = modseq_dat |>
      summarize(mean_m_frac = sum(cov * m_frac) / sum(cov))
  }
  
  print("stuck 4")
  if ("h_frac" %in% colnames(modseq_dat)) {
    modseq_dat = modseq_dat |>
      summarize(mean_h_frac = sum(cov * h_frac) / sum(cov))
  }
}
