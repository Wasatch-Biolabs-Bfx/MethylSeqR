
get_mod_stats <- function(modseq_dat, 
                          plot = FALSE) 
{
  # decide if per base or per region 
  regional_dat = "region_name" %in% colnames(modseq_dat)
  
  # grab mh frac info- prioritize mh_frac over m_frac
  if ("mh_frac" %in% colnames(modseq_dat)) {
    goodMeth = 100 * pull(modseq_dat, mh_frac)
  } else if ("mean_mh_frac" %in% colnames(modseq_dat)) {
    goodMeth = 100 * pull(modseq_dat, mean_mh_frac)
  } else if ("m_frac" %in% colnames(modseq_dat)) {
    goodMeth = 100 * pull(modseq_dat, m_frac)
  } else if ("mean_m_frac" %in% colnames(modseq_dat)) {
    goodMeth = 100 * pull(modseq_dat, mean_m_frac)
  }
  
  qts <- c(seq(0, 0.9, 0.1), 0.95, 0.99, 0.995, 0.999, 1) 
  
  if (!plot) { # if only stats are wanted
    title <- "Methylation statistics per base\n"
    if (regional_dat) {
      title <- "Methylation statistics per region\n"
    }
    
    cat(title)
    cat("Summary:\n")
    print(summary(goodMeth, p=qts))
    cat("percentiles:\n")
    print(quantile(goodMeth, p=qts))
    cat("\n")
    
  } else if (plot) { # if they want a plot
    x_title <- "% methylation per base"
    if (regional_dat) {
      x_title <- "% methylation per region"
    }
    
    # Create a data frame from your list
    plot <- data.frame(methylation_value = goodMeth)
    
    # Create the histogram
    ggplot(plot, aes(x = methylation_value)) +
    geom_histogram(binwidth = 10, fill = "cornflowerblue", 
                   color = "black", linewidth = 0.25) +
    labs(title = "Histogram of % CpG Methylation", 
         x = x_title, y = "Frequency") +
    theme_minimal()
  }
}