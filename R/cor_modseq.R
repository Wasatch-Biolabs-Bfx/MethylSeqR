#' Compute and Optionally Plot Correlation Matrix of Methylation Data
#'
#' This function computes the Pearson correlation matrix of methylation data
#' and optionally plots the correlation matrix.
#'
#' @param modseq_dat A data frame containing methylation data. The data frame can
#' either contain positional or regional data. If it contains regional data, it must
#' have a column named \code{region_name}.
#'
#' @param plot A logical value indicating whether to plot the correlation matrix.
#' Default is \code{FALSE}.
#'
#' @return The correlation matrix of the methylation data.
#'
#' @examples
#' \dontrun{
#' correlation_matrix <- cor_modseq(modseq_dat, plot = TRUE)
#' }
#'
#' @import dplyr tidyr ggplot2 reshape2
#'
#' @importFrom stats cor
#'
#' @export
cor_modseq <- function(modseq_dat,
                       plot = FALSE)
{
  # decide if regional or positional data
  regional_dat = "region_name" %in% colnames(modseq_dat)

  if (regional_dat) {
    # Aggregate mean_mh_frac by sample and region_name
    dat_wide <- modseq_dat |>
      pivot_wider(names_from = sample_name,
                  values_from = mean_mh_frac)

    # Compute Correlation
    sample_names <- unique(modseq_dat$sample_name)
    sample_names <- sample_names[!is.na(sample_names)]  # Remove NA values
    numeric_columns <- dat_wide[, sample_names, drop = FALSE]

    # Calculate correlation matrix
    correlation_matrix <- cor(numeric_columns,
                              use = "pairwise.complete.obs",
                              method = "pearson")

    # View the correlation matrix
    print(correlation_matrix)
  } else {
    dat_wide <- modseq_dat |>
      mutate(chr_pos = paste(chrom,
                             ref_position,
                             sep = "_")) |>
      pivot_wider(id_cols = chr_pos,
                  names_from = sample_name,
                   values_from = mh_frac)


    #slice_sample(n = 5000000) |>
    # just grab a single chromosme

    print(head(dat_wide))

    # Compute Correlation
    numeric_columns <- dat_wide[, unique(modseq_dat$sample_name)]
    print(head(numeric_columns))

    # Calculate correlation matrix
    correlation_matrix <- cor(numeric_columns,
                              use = "pairwise.complete.obs",
                              method = "pearson")

    # # View the correlation matrix
    print(correlation_matrix)
  }


  # Plot Correlation Matrix with Correlation Values
  if (plot) {
    melted_cor <- melt(correlation_matrix)

    print(ggplot(data = melted_cor, aes(x = Var1,
                                  y = Var2,
                                  fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)),
              color = "black") + # Add correlation values
    scale_fill_gradient2(low = "blue",
                         high = "red",
                         mid = "white",
                         midpoint = 0,
                         limits = c(-1, 1),
                         space = "Lab",
                         name = "Pearson\nCorrelation") +
    scale_y_discrete(limits = rev(levels(factor(melted_cor$Var2)))) +
    theme_minimal() +
    labs(title = "Sample Correlation Matrix",
         x = "Sample",
         y = "Sample") +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1)))
  }
}
