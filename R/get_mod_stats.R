#' Calculate and Plot Methylation Statistics
#'
#' This function calculates and optionally plots statistics for methylation data from
#' sequencing experiments. It can handle both positional and regional methylation data.
#'
#' @param modseq_dat A data frame containing methylation data. It should have one of
#' the following columns: \code{mh_frac}, \code{mean_mh_frac}, \code{m_frac}, or \code{mean_m_frac}.
#'
#' @param plot Logical, if \code{TRUE}, the function will generate a histogram of
#' the methylation data. Default is \code{FALSE}.
#'
#' @return If \code{plot} is \code{FALSE}, the function prints summary statistics
#' and percentiles of the methylation data. If \code{plot} is \code{TRUE}, it prints a
#' histogram of the methylation data.
#'
#' @examples
#' \dontrun{
#' get_mod_stats(modseq_dat)
#' get_mod_stats(modseq_dat, plot = TRUE)
#' }
#'
#' @import dplyr ggplot2
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal
#'
#' @export
get_mod_stats <- function(ch3_db,
                          call_type = c("positions", "regions"),
                          plot = FALSE)
{
  
  # If a character file name is provided, then make ch3 class obj
  db_con = helper_connectDB(ch3_db)
  
  if (length(call_type) > 1) {
    call_type = c("positions")
  }
  
  # Check for specific table and connect to it in the database
  if (!dbExistsTable(db_con, call_type)) {
    stop(paste0(call_type, " Table does not exist. You can create it by..."))
  }
  
  modseq_dat = tbl(db_con, call_type) 
  
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
    print(ggplot(plot, aes(x = methylation_value)) +
    geom_histogram(binwidth = 10, fill = "cornflowerblue",
                   color = "black", linewidth = 0.25) +
    labs(title = "Histogram of % CpG Methylation",
         x = x_title, y = "Frequency") +
    theme_minimal())
  }
  
  # Finish up: close the connection
  dbDisconnect(db_con, shutdown = TRUE)
}
