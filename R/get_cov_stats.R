#' Calculate and Plot Coverage Statistics
#'
#' This function calculates and optionally plots statistics for coverage data from
#' methylation sequencing experiments. It can handle both positional and regional
#' methylation data.
#'
#' @param ch3_db A data base either linking to the file name or of class ch3_db.
#' @param call_type Either positions or regions data to analyze coverage on.
#' @param plot Logical, if \code{TRUE}, the function will generate a histogram of
#' the coverage data. Default is \code{FALSE}.
#' @param save_path Pathway to save the plot to. Usually .pdf or .png.
#'
#' @return If \code{plot} is \code{FALSE}, the function prints summary statistics
#' and percentiles of the coverage data. If \code{plot} is \code{TRUE}, it prints a
#' histogram of the log-transformed coverage data.
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylseqR")
#'  
#'  # Get coverage statistics for the 'positions' call type without plotting
#'  get_cov_stats(ch3_db = ch3_db, call_type = "positions")
#'  
#'  # Get coverage statistics for the 'positions' call type with plotting
#'  get_cov_stats(ch3_db = ch3_db, call_type = "positions", plot = TRUE)
#'
#' @import dplyr ggplot2
#'
#' @export

get_cov_stats <- function(ch3_db,
                          call_type = "positions",
                          plot = FALSE,
                          save_path = NULL)
{
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  tryCatch(
    {
      if (length(call_type) > 1) {
        call_type = c("positions")
      }
      
      # Check for specific table and connect to it in the database
      if (!dbExistsTable(db_con, call_type)) {
        stop(paste0(call_type, " Table does not exist. You can create it by..."))
      }
      
      modseq_dat = tbl(db_con, call_type) 
      
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
        p <- ggplot(plot, aes(x = coverage)) +
                geom_histogram(binwidth = 0.25, fill = "chartreuse4",
                               color = "black", linewidth = 0.25) +
                labs(title = "Histogram of CpG Coverage",
                     x = x_title, y = "Frequency") +
                theme_minimal()
        print(p)
        
        # Save the plot if save_path is specified
        if (!is.null(save_path)) {
          ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
          cat("Coverage plot saved to ", save_path, "\n")
        }
      }
    }, 
    error = function(e) 
    {
      # Print custom error message
      message("An error occurred: ", e$message)
    }, 
    finally = 
      {
        # Finish up: close the connection
        .helper_closeDB(database)
      })
}
