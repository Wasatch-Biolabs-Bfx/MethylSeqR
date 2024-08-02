#' Generate and Plot Density of Methylation Data
#'
#' This function generates density plots for pairs of methylation samples
#' and visualizes them using raster plots.
#'
#' @param modseq_dat A data frame containing methylation data. The data frame can
#' either contain positional or regional data. If it contains regional data, it must
#' have a column named \code{region_name}.
#'
#' @param samples A character vector specifying which samples to compare. If set to \code{"all"},
#' all samples in the data will be compared. If fewer than two samples are specified, the function will
#' stop execution with an error message.
#'
#' @param save A logical value indicating whether to save the plots to PNG files. If \code{TRUE},
#' the plots will be saved with filenames based on the sample names.
#'
#' @return A list of ggplot objects representing the density plots for each pair of samples, if \code{save}
#' is \code{FALSE}. If \code{save} is \code{TRUE}, the plots are saved as PNG files and the function
#' returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#' density_modseq(modseq_dat)
#' density_modseq(your_data, samples = c("barcode81", "barcode84", "barcode87"))
#' }
#'
#' @import dplyr tidyr ggplot2 MASS
#'
#' @importFrom MASS kde2d
#'
#' @export
density_modseq <- function(modseq_dat, samples = c("all"), save = FALSE) {
  # Decide if per base or per region
  regional_dat <- "region_name" %in% colnames(modseq_dat)

  if (regional_dat) {
    plotting_data <- na.omit(modseq_dat) %>%
      dplyr::select(c(region_name, sample_name, mh_frac)) %>%
      pivot_wider(names_from = sample_name, values_from = mh_frac) %>%
      na.omit()
  } else {
    plotting_data <- modseq_dat %>%
      slice_sample(n = 5000000) %>%
      mutate(chr_pos = paste(chrom, ref_position, sep = "_")) %>%
      pivot_wider(id_cols = chr_pos, names_from = sample_name, values_from = mh_frac) %>%
      na.omit()
  }

  # Filter samples if specified
  sample_cols <- colnames(plotting_data)[-1]
  if (!("all" %in% samples)) {
    sample_cols <- intersect(sample_cols, samples)
  }

  # Check if only one sample is provided
  if (length(sample_cols) < 2) {
    stop("At least two samples must be provided for comparison.")
  }

  # Get all pairs of sample columns
  pairs <- combn(sample_cols, 2, simplify = FALSE)

  for (pair in pairs) {
    sample1 <- pair[1]
    sample2 <- pair[2]

    kde <- kde2d(plotting_data[[sample1]], plotting_data[[sample2]], n = 100, h = 0.67)

    density_data <- expand.grid(x = kde$x, y = kde$y)
    density_data$z <- as.vector(kde$z)

    plot <- ggplot(density_data, aes(x = x, y = y, fill = z)) +
      geom_raster() +
      scale_fill_viridis_c(option = "turbo") +
      labs(
        x = paste(sample1, "Methylation Signal"),
        y = paste(sample2, "Methylation Signal"),
        fill = "Density",
        title = paste("Density Plot:", sample1, "vs", sample2)
      ) +
      theme_minimal() +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 12, hjust = 0.5))

    print(plot)

    if (save) {
      ggsave(paste0("density_plot_", sample1, "_vs_", sample2, ".png"), plot = plot, height = 6, width = 8)
    }
  }
}

