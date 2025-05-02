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
#' @param max_rows The maximum amount of rows wanted for calculation. This argument can help analysis run faster when there is a lot of data.
#'
#' @return If \code{plot} is \code{FALSE}, the function prints summary statistics
#' and percentiles of the coverage data. If \code{plot} is \code{TRUE}, it prints a
#' histogram of the log-transformed coverage data.
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Get coverage statistics for the 'positions' call type without plotting
#'  plot_ch3_cov(ch3_db = ch3_db, call_type = "positions")
#'
#' @importFrom DBI dbConnect dbDisconnect dbExistsTable
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl summarise pull sql
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal ggsave
#' @importFrom stats quantile na.omit
#'
#' @export

plot_ch3_cov <- function(ch3_db,
                          call_type = "positions",
                          plot = TRUE,
                          save_path = NULL,
                         max_rows = NULL)
{
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  if (length(call_type) > 1) {
    call_type = c("positions")
  }
  
  # Check for specific table and connect to it in the database
  if (!dbExistsTable(ch3_db$con, call_type)) {
    stop(paste0(call_type, " Table does not exist in the database. Check spelling or make sure you create it first.\n"))
  }
  
  # Determine total number of rows first
  total_rows <- tbl(ch3_db$con, call_type) |> summarise(n = n()) |> pull(n)
  
  # Sample in SQL if max_rows is given and valid
  if (!is.null(max_rows)) {
    if (max_rows > total_rows) {
      stop(paste0("Requested max_rows (", max_rows, 
                  ") exceeds available rows in the table (", total_rows, ")."))
    }
    
    modseq_dat <- tbl(ch3_db$con, sql(paste0(
      "SELECT * FROM ", call_type, 
      " USING SAMPLE ", max_rows, " ROWS"
    )))
  } else {
    modseq_dat <- tbl(ch3_db$con, call_type)
  }
  
  # Checks
  stopifnot("Invalid dataframe format. A 'num_calls' or 'mean_num_calls' column must be present." =
              any(c("num_calls", "mean_num_calls") %in% colnames(modseq_dat)))
  
  # Clean dataframe
  modseq_dat <- na.omit(modseq_dat)
  
  # Decide if per base or per region
  regional_dat = "region_name" %in% colnames(modseq_dat)
  
  # if (!regional_dat) {
  num_calls = pull(modseq_dat, num_calls)
  # } else {
  #   num_calls = pull(modseq_dat, mean_num_calls)
  # }
  
  qts <- c(seq(0, 0.9, 0.1), 0.95, 0.99, 0.995, 0.999, 1)
  
  if (!plot) {
    title <- "read coverage statistics per base\n"
    
    if (regional_dat) {
      title <- "read coverage statistics per region\n"
    } else if (call_type == "windows") {
      title <- "read coverage statistics per window\n"
    }
    
    cat(title)
    cat("summary:\n")
    print( summary( num_calls ) )
    cat("percentiles:\n")
    print(quantile(num_calls, p=qts ))
    cat("\n")
  } else {
    x_title <- "log10 of read coverage per base"
    if (regional_dat) {
      x_title <- "log10 of read coverage per region"
    } else if (call_type == "windows") {
      x_title <- "log10 of read coverage per window"
    }
    
    # Create a data frame from your list
    plot <- data.frame(coverage = log10(num_calls))
    
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
  
  ch3_db <- .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}