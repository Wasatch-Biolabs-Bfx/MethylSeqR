#' Calculate and Plot Coverage Statistics
#'
#' This function calculates and optionally plots statistics for coverage data from
#' methylation sequencing experiments. It can handle both positional and regional
#' methylation data.
#'
#' @param modseq_dat A data frame containing methylation data. It should have a
#' column named \code{cov} for positional data or \code{mean_cov} for regional data.
#'
#' @param plot Logical, if \code{TRUE}, the function will generate a histogram of
#' the coverage data. Default is \code{FALSE}.
#'
#' @return If \code{plot} is \code{FALSE}, the function prints summary statistics
#' and percentiles of the coverage data. If \code{plot} is \code{TRUE}, it prints a
#' histogram of the log-transformed coverage data.
#'
#' @examples
#' \dontrun{
#' get_cov_stats(modseq_dat)
#' get_cov_stats(modseq_dat, plot = TRUE)
#' }
#'
#' @import dplyr ggplot2
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal
#'
#' @export
get_cov_stats <- function(modseq_dat,
                          plot = FALSE)
{
  # Checks
  stopifnot("Invalid dataframe format. A 'cov' or 'mean_cov' column must be present." =
              any(c("cov", "mean_cov") %in% colnames(modseq_dat)))

  # Clean dataframe
  modseq_dat <- na.omit(modseq_dat)

  # Decide if per base or per region
  regional_dat = "region_name" %in% colnames(modseq_dat)

  # if (!regional_dat) {
    cov = pull(modseq_dat, cov)
  # } else {
  #   cov = pull(modseq_dat, mean_cov)
  # }

  qts <- c(seq(0, 0.9, 0.1), 0.95, 0.99, 0.995, 0.999, 1)

  if (!plot) {
    title <- "read coverage statistics per base\n"

    if (regional_dat) {
      title <- "read coverage statistics per region\n"
    }

    cat(title)
    cat("summary:\n")
    print( summary( cov ) )
    cat("percentiles:\n")
    print(quantile(cov, p=qts ))
    cat("\n")
  } else {
    x_title <- "log10 of read coverage per base"
    if (regional_dat) {
      x_title <- "log10 of read coverage per region"
    }

    # Create a data frame from your list
    plot <- data.frame(coverage = log10(cov))

    # Create the histogram
    print(ggplot(plot, aes(x = coverage)) +
      geom_histogram(binwidth = 0.25, fill = "chartreuse4",
                     color = "black", linewidth = 0.25) +
      labs(title = "Histogram of CpG Coverage",
           x = x_title, y = "Frequency") +
      theme_minimal())
  }
}
