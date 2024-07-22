cluster_modseq <- function(modseq_dat)
{
  # Decide if regional or positional data
  regional_dat <- "region_name" %in% colnames(modseq_dat)
  
  if (regional_dat) {
    data_matrix <- 
      na.omit(modseq_dat) |>
      select(
        c(region_name, sample_name, mean_mh_frac)) |>
      pivot_wider(
        names_from = sample_name, 
        values_from = mean_mh_frac) |>
      select(
        -region_name) |>
      as.matrix()
  } else {
    data_matrix <- 
      na.omit(modseq_dat) |>
      slice_sample(
        n = 5000000) |>
      mutate(
        chr_pos = paste(chrom, ref_position, sep = "_")) |>
      pivot_wider(
        names_from = sample_name, values_from = mh_frac) |>
      select(
        -cov, -chrom, -ref_position) |>
      select(-chr_pos) |>
      as.matrix()
  }
  
  # Compute the distance matrix
  distance_matrix <- dist(t(data_matrix))
  
  # Perform hierarchical clustering
  hc <- hclust(na.omit(distance_matrix), method = "complete")
  
  # Convert the clustering object to a dendrogram
  dend <- as.dendrogram(hc)
  
  # Plot the dendrogram using base R plot
  plot(dend, main = "Dendrogram of Methylation Samples")
}
