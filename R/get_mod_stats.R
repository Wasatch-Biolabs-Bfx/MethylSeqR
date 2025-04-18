#' Get Methylation Statistics from the ch3 Database
#'
#' This function retrieves and calculates methylation statistics (mean methylation fractions) 
#' from a specified table in the ch3 database. It can either return summary statistics or plot 
#' a histogram of the methylation values, depending on the user's preference. 
#'
#' @param ch3_db A string. The path to the database containing ch3 files from nanopore data.
#' @param call_type A character vector specifying the type of data to retrieve from the database. 
#'                  Default is "positions". Can also be "regions".
#' @param plot A logical value. If TRUE, a histogram of methylation values is plotted. Default is FALSE.
#' @param save_path Pathway to save the plot to. Usually .pdf or .png.
#'
#' @details
#' The function connects to the specified database, checks for the existence of the relevant table, 
#' and retrieves methylation fraction data. If the table contains methylation data, it prioritizes 
#' the `mh_frac` column over others. Depending on the `call_type`, it can compute statistics 
#' for either per base or per region. If `plot` is set to TRUE, it generates a histogram of the 
#' methylation values.
#'
#' @note The function assumes that the database has tables named according to the `call_type` parameter 
#' (e.g., "positions", "regions"). It also expects specific columns for methylation data to exist.
#'
#' @import DBI dplyr ggplot2
#'
#' @return If `plot` is FALSE, the function prints a summary of the methylation statistics and quantiles. 
#' If `plot` is TRUE, it displays a histogram of methylation values.
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Get methylation statistics for the 'positions' call type without plotting
#'  get_mod_stats(ch3_db = ch3_db, call_type = "positions")
#'
#' @export
get_mod_stats <- function(ch3_db,
                          call_type = c("positions", "regions"),
                          plot = TRUE,
                          save_path = NULL)
{
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # Finish up: update table list and close the connection
  on.exit(.helper_closeDB(database), add = TRUE)
  
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
    } else if (call_type == "windows") {
      title <- "Methylation statistics per window\n"
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
    } else if (call_type == "windows") {
      x_title <- "% methylation per window"
    }
    
    # Create a data frame from your list
    plot <- data.frame(methylation_value = goodMeth)
    
    # Create the histogram
    p <- ggplot(plot, aes(x = methylation_value)) +
            geom_histogram(binwidth = 10, fill = "cornflowerblue",
                           color = "black", linewidth = 0.25) +
            labs(title = "Histogram of % CpG Methylation",
                 x = x_title, y = "Frequency") +
            theme_minimal()
    print(p)
    
    # Save the plot if save_path is specified
    if (!is.null(save_path)) {
      ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
      cat("Statistics plot saved to ", save_path, "\n")
      }
  }
}
