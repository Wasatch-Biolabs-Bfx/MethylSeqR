get_cov_stats <- function(modseq_dat, plot = FALSE) {
  # clean dataframe
  modseq_dat <- na.omit(modseq_dat)
  # decide if per base or per region
  regional_dat = "region_name" %in% colnames(modseq_dat)

  if (!regional_dat) {
    cov = modseq_dat$cov
  } else {
    cov = modseq_dat$mean_cov
  }

  qts=seq(0,0.9,0.1) # get quantiles
  qts=c(qts,0.95,0.99,0.995,0.999,1)

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

  } else if (plot) { # lets make a lil plot
    x_title <- "log10 of read coverage per base"
    if (regional_dat) {
      x_title <- "log10 of read coverage per region"
    }

    # Create a data frame from your list
    plot <- data.frame(coverage = log10(cov))
    # Create the histogram
    ggplot(plot, aes(x = coverage)) +
      geom_histogram(binwidth = 0.25, fill = "chartreuse4", color = "black", linewidth = 0.25) +
      labs(title = "Histogram of CpG Coverage", x = x_title, y = "Frequency") +
      theme_minimal()

  }

}
