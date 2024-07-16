cluster_modseq <- function(modseq_dat){
  #library(dplyr)
  #library(dendextend)
  
  # decide if regional or positional data
  regional_dat = "region_name" %in% colnames(modseq_dat)
  
  if (regional_dat) {
    test_wide <- na.omit(modseq_dat) |>
      select(c(region_name, sample_name, mean_mh_frac)) |>
      pivot_wider(names_from = sample_name, 
                  values_from = mean_mh_frac)
    # Remove the 'region' column to focus on methylation values
    data_matrix <- as.matrix(test_wide |> 
                               select(-region_name))
  } else {
     test_wide <- modseq_dat %>%
      slice_sample(n = 5000000) %>%
      mutate(chr_pos = paste(chrom, ref_position, sep = "_")) %>%
      pivot_wider(names_from = sample_name, values_from = mh_frac) %>%
      dplyr::select(-cov, -chrom, -ref_position) %>%
      na.omit()
     
      # Remove the 'region' column to focus on methylation values
      data_matrix <- as.matrix(test_wide %>% dplyr::select(-chr_pos))
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
