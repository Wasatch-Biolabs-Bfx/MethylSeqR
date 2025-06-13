#' Plot Differential Methylation Volcano Raster
#'
#' Creates a volcano-style raster plot of differential methylation results using database-backed plotting via `dbplot`. 
#' This function connects to a DuckDB database, retrieves the specified differential methylation table, filters and 
#' transforms the data, and generates a raster plot with `meth_diff` on the x-axis and `-log10(p-value)` on the y-axis.
#'
#' @param ch3_db A `ch3_db` object or a character string representing the file path to a DuckDB database.
#'        The database must contain a table with differential methylation results.
#' @param table A character string specifying the name of the table in the database containing the differential methylation data.
#'        The table must contain at least the following columns: `meth_diff` and `p_val`.
#'
#' @return Invisibly returns the `ch3_db` object after closing the database connection. The function prints the ggplot2 raster plot to the active R graphics device.
#' 
#' @details
#' The plot uses `dbplot::dbplot_raster()` to efficiently create a raster visualization of large-scale methylation difference data. 
#' It applies a log10 transformation to the `p_val` column and uses a color gradient to show the density of observations in each bin.
#' A message is printed to indicate the time taken to generate the plot.
#'
#' @examples
#' \dontrun{
#' plot_ch3_diff("my_methylation.ch3.db", "mod_diff_windows")
#' }
#'
#' @importFrom dplyr filter mutate
#' @importFrom dbplot dbplot_raster
#' @importFrom ggplot2 scale_fill_viridis_c labs theme_minimal
#' @importFrom DBI dbExistsTable
#' @export
plot_ch3_diff <- function(ch3_db,
                          table) {
  start_time <- Sys.time()
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  # check for differential methylation table
  if (!dbExistsTable(ch3_db$con, table)) { # add db_con into object and put in every function...
    stop(paste0(table, " table does not exist. Build it with calc_ch3_diff()!"))
  }
  
  # Connect to the table
  tbl_diff <- tbl(ch3_db$con, table)
  
  # Plot using dbplot_raster
  plot <- tbl(ch3_db$con, table) |>
    filter(
      !is.na(meth_diff),
      !is.nan(meth_diff),
      !is.na(p_val),
      p_val > 0
    ) |>
    mutate(log_p = -log10(p_val)) |>
    dbplot::dbplot_raster(meth_diff, log_p)
  
  plot <- plot +
    scale_fill_viridis_c(option = "rocket", direction = -1, begin = 0.2, end = 0.9) +  # Fancy color gradient- add option = "mako"
    labs(
      title = "Differential Methylation Volcano Plot",
      x = "Methylation Difference (case-control)",
      y = "-log10(p-value)",
      fill = "Count"
    ) +
    theme_minimal()  #+
    #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    #annotate("text", x = 0.9, y = -log10(0.05) + 5, label = "p = 0.05", color = "red", size = 4)
  
  
  print(plot)

  end_time <- Sys.time()
  
  total_time_difftime <- end_time - start_time
  
  # Convert the total_time_difftime object to numeric seconds for a reliable comparison
  total_seconds <- as.numeric(total_time_difftime, units = "secs")
  
  if (total_seconds > 60) {
    # If greater than 60 seconds, convert to numeric minutes for display
    total_minutes <- as.numeric(total_time_difftime, units = "mins")
    message("Differential Methylation Plotted!",
            "\nTime elapsed: ", round(total_minutes, 2), " minutes\n")
  } else {
    # Otherwise, display in numeric seconds
    message("Differential Methylation Plotted!",
            "\nTime elapsed: ", round(total_seconds, 2), " seconds\n")
  }
  ch3_db <- .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}