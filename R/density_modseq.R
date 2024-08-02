#' Generate and Plot Density of Methylation Data
#'
#' This function generates density plots for pairs of methylation samples
#' and visualizes them using raster plots.
#'
#' @param modseq_dat A data frame containing methylation data. The data frame can
#' either contain positional or regional data. If it contains regional data, it must
#' have a column named \code{region_name}.
#'
#' @return None. The function prints density raster plots for pairs of methylation samples.
#'
#' @examples
#' \dontrun{
#' density_modseq(modseq_dat)
#' }
#'
#' @import dplyr tidyr ggplot2 MASS
#'
#' @importFrom MASS kde2d
#'
#' @export
density_modseq <- function(modseq_dat) {
  # decide if per base or per region
  regional_dat = "region_name" %in% colnames(modseq_dat)

  if (regional_dat) {

    plotting_data <- na.omit(modseq_dat) |>
      dplyr::select(c(region_name,
                      sample_name,
                      mh_frac)) |>
      pivot_wider(names_from = sample_name,
                  values_from = mh_frac) |>
      na.omit()
  } else {
    plotting_data <- modseq_dat |>
      slice_sample(n = 5000000) |>
      mutate(chr_pos = paste(chrom,
                             ref_position,
                             sep = "_")) |>
      pivot_wider(id_cols = chr_pos,
                  names_from = sample_name,
                  values_from = mh_frac) |>
      na.omit()
  }
  # Get all pairs of sample columns
  sample_cols <- colnames(plotting_data)[-1]
  # Exclude the first column (region_name or chr_pos)
  pairs <- combn(sample_cols, 2, simplify = FALSE)
  # Prepare an empty list to store density data for each pair
  all_density_data <- list()

  for (pair in pairs) {
    # go through each set of samples, create a density plot for each
    sample1 <- pair[1]
    sample2 <- pair[2]

    kde <- kde2d(plotting_data[[sample1]],
                 plotting_data[[sample2]],
                 n = 100,
                 h = 0.67)
    # Adjust h parameter if needed

    # Prepare data for plotting
    density_data <- expand.grid(x = kde$x, y = kde$y)
    density_data$z <- as.vector(kde$z)
    density_data$sample1 <- sample1
    density_data$sample2 <- sample2

    all_density_data[[paste(sample1, sample2, sep = "_")]] <- density_data
  }

  # Combine all density data into one data frame
  combined_density_data <- bind_rows(all_density_data)

  # Plot using ggplot2 and geom_raster
  print(ggplot(combined_density_data, aes(x = x,
                                          y = y,
                                          fill = z)) +
    geom_raster() +
    scale_fill_viridis_c(option = "turbo") +  # Using a color scale from the viridis package
    labs(
      x = paste(sample1, "Methylation Signal"),
      y = paste(sample2, "Methylation Signal"),
      fill = "Density",
      title = "Raster Plot of Methylation Signals"
    ) +
    theme_minimal() +
    facet_wrap(~ sample1 + sample2,
               scales = "free"))
}
