pca_modseq <- function(ch3_db, call_type = "positions") {
  # Connect to the database
  db_con <- helper_connectDB(ch3_db)
  
  # Retrieve the data from the database
  modseq_dat <- tbl(db_con, call_type) %>% collect()  # Collect data to bring it into memory
  
  # Omit any missing values
  modseq_dat <- na.omit(modseq_dat)
  
  if (call_type == "regions") {
    # Aggregate mean_mh_frac by sample and region_name
    test_wide <- modseq_dat %>%
      select(c(region_name, sample_name, mh_frac)) %>%
      pivot_wider(names_from = sample_name, values_from = mh_frac) %>%
      na.omit() %>%
      as.data.frame()  # Convert to dataframe
  } else {
    # Aggregate mean_mh_frac by chr_pos and sample_name
    test_wide <- modseq_dat %>%
      mutate(chr_pos = paste(chrom, ref_position, sep = "_")) %>%
      pivot_wider(id_cols = chr_pos, names_from = sample_name, values_from = mh_frac) %>%
      na.omit() %>%
      as.data.frame()  # Convert to dataframe
  }
  
  # Ensure the collected data has the correct structure
  if (ncol(test_wide) <= 1) {
    stop("The data doesn't have enough columns for PCA after processing.")
  }
  
  # Scale the data: remove the first column (e.g., chr_pos or region_name)
  scaled_data <- scale(test_wide[, -1])
  
  # Perform PCA
  pca_result <- prcomp(t(scaled_data))  # Transpose data because samples should be rows
  pca_summary <- summary(pca_result)
  print(pca_summary)
  
  # Extract variance explained by PC1 and PC2
  pc1_var <- round(pca_summary$importance[2, 1] * 100, 2)
  pc2_var <- round(pca_summary$importance[2, 2] * 100, 2)
  
  # Prepare PCA data for plotting
  pca_data <- data.frame(pca_result$x)
  pca_data$sample_name <- rownames(pca_data)  # Add sample names
  
  print(pca_data)
  
  # Plot PCA
  print(ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_name)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = "PCA Plot of Methylation Data",
      x = paste0("PC1 (", pc1_var, "% variance)"),
      y = paste0("PC2 (", pc2_var, "% variance)")
    ))
  
  # Finish up: close the connection
  dbDisconnect(db_con, shutdown = TRUE)
}
