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

  if (!regional_dat) {
    cov = pull(modseq_dat, cov)
  } else {
    cov = pull(modseq_dat, mean_cov)
  }

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
    ggplot(plot, aes(x = coverage)) +
      geom_histogram(binwidth = 0.25, fill = "chartreuse4", 
                     color = "black", linewidth = 0.25) +
      labs(title = "Histogram of CpG Coverage", 
           x = x_title, y = "Frequency") +
      theme_minimal()
  }
}
