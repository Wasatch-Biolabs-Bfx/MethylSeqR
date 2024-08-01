#' Perform Principal Component Analysis (PCA) on Methylation Data
#'
#' This function performs PCA on methylation data to reduce dimensionality
#' and visualize the variance explained by the principal components.
#'
#' @param modseq_dat A data frame containing methylation data. It should include
#'                    sample names, methylation fractions, and may include region names.
#'
#' @return A PCA plot showing the first two principal components and the
#'         variance explained by each component. The PCA summary is also printed
#'         to the console.
#'
#' @examples
#' pca_modseq(data)
#'
#' @importFrom dplyr na.omit select mutate
#'
#' @importFrom tidyr pivot_wider
#'
#' @importFrom stats prcomp
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal
#'
#' @export
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
      dplyr::select(c(region_name, sample_name, mh_frac)) |>
      pivot_wider(names_from = sample_name,
                  values_from = mh_frac) |>
      na.omit()
  } else {
    # Aggregate mean_mh_frac by sample and region_name
    test_wide <-
      modseq_dat |>
      # slice_sample(n = 500000) |>
      mutate(chr_pos = paste(chrom, ref_position, sep = "_")) |>
      pivot_wider(id_cols = chr_pos,
                  names_from = sample_name,
                  values_from = mh_frac) |>
      na.omit()
  }

  # clean dataframe & remove the first region name column or position column...
  scaled_data <- scale(test_wide[ , -1])

  # Now you can perform PCA
  pca_result <- prcomp(t(scaled_data))
  pca_summary <- summary(pca_result)
  print(pca_summary)

  # Extract variance explained by PC1 and PC2
  pc1_var <- round(pca_summary$importance[2, 1] * 100, 2)
  pc2_var <- round(pca_summary$importance[2, 2] * 100, 2)

  # Plot PCA
  # Prepare PCA data for plotting
  pca_data <- data.frame(pca_result$x)
  pca_data$sample_name <- rownames(pca_data)

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
}
