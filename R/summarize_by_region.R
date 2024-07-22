aggregate_regions <- function(modseq_dat, 
                              annot_file)
{
  # Read annotation
  annotation <- 
    read_tsv(annot_file, 
             col_names = c("chrom", "start", "end", "region_name"), 
             show_col_types = FALSE) 
  
  stopifnot("Invalid annotation format. File must have chr, start, end." = 
              ncol(annotation) >= 3)
  
  if (ncol(annotation == 3)) {
    annotation <- 
      annotation |>
      mutate(
        region_name = paste(chrom, start, end, sep = "_"))
  }
  
  annotation <- 
    annotation |>
    reframe(
      .by = c(region_name, chrom),
      ref_position = start:end)
    
  # Create regional dataframe
  modseq_dat |>
    right_join(
      annotation, 
      by = join_by(chrom, ref_position),
      copy = TRUE) |>
    summarize( 
      .by = c(sample_name, region_name),
      cov = sum(cov),
      across(ends_with("_counts"), sum),
      across(ends_with("_frac"), ~ sum(.x * cov) / sum(cov)))
}
