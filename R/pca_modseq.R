pca_modseq <- function(modseq_dat) 
{
  # clean
  modseq_dat = na.omit(modseq_dat)
  
  # decide if per base or per region 
  regional_dat = "region_name" %in% colnames(modseq_dat)
  
  if (regional_dat) {
    # Aggregate mean_mh_frac by sample and region_name
    test_wide <- 
      modseq_dat |>
      select(c(region_name, sample_name, mean_mh_frac)) |>
      pivot_wider(names_from = sample_name, values_from = mean_mh_frac) |> 
      na.omit()
  } else {
    # Aggregate mean_mh_frac by sample and region_name
    test_wide <- 
      modseq_dat |>
      slice_sample(n = 500000) |>
      mutate(chr_pos = paste(chrom, ref_position, sep = "_")) |>
      pivot_wider(id_cols = chr_pos, names_from = sample, 
                  values_from = mh_frac) |>
      na.omit()
  }
  
  # clean dataframe & remove the first region name column or position column...
  scaled_data <- scale(test_wide[ , -1])
  
  # Now you can perform PCA
  pca_result <- prcomp(t(scaled_data))
  print(summary(pca_result))

  # Plot PCA
  # Prepare PCA data for plotting
  pca_data <- data.frame(pca_result$x)
  pca_data$sample <- rownames(pca_data)
  
  print(pca_data)

  # Plot PCA
  ggplot(pca_data, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot of Methylation Data", 
       x = "Principal Component 1", y = "Principal Component 2")
}
