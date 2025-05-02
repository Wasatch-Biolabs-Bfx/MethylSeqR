#' Perform PCA on Methylation Data
#'
#' This function performs Principal Component Analysis (PCA) on methylation data retrieved from a DuckDB database.
#' It aggregates the methylation fraction data based on the specified call type and prepares it for PCA analysis.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param call_type A string representing the name of the table in the database from which to pull the data. 
#' Default is "positions".
#' @param save_path Pathway to save the plot to. Usually .pdf or .png.
#' @param max_rows The maximum amount of rows wanted for calculation. This argument can help analysis run faster when there is a lot of data.
#'
#' @details
#' The function connects to the specified DuckDB database, retrieves the methylation data from the specified call type table, 
#' and aggregates the data either by regions or chromosomal positions. PCA is then performed on the aggregated data, and 
#' a scatter plot of the first two principal components is generated.
#'
#' @return A PCA plot is produced, showing the first two principal components of the methylation data.
#' The function also prints a summary of the PCA results and the PCA data used for plotting.
#' 
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  # Calculate PCA
#'  plot_ch3_pca(ch3_db)
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl select mutate collect
#' @importFrom tidyr pivot_wider
#' @importFrom stats prcomp na.omit
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs ggsave
#' 
#' @export

plot_ch3_pca <- function(ch3_db, 
                       call_type = "positions",
                       save_path = NULL,
                       max_rows = NULL) {
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)

  # If max_rows is specified, check table size and sample rows randomly in SQL
  if (!is.null(max_rows)) {
    row_count <- dbGetQuery(ch3_db$con, paste0("SELECT COUNT(*) as n FROM ", call_type))$n
    
    if (row_count < max_rows) {
      stop(paste0("Table '", call_type, "' only has ", row_count, " rows, which is fewer than max_rows = ", max_rows, ". Pick less rows."))
    }
    
    # Use SQL random sampling with ORDER BY RANDOM()
    modseq_dat <- dbGetQuery(ch3_db$con, paste0("SELECT * FROM ", call_type, " ORDER BY RANDOM() LIMIT ", max_rows))
  } else {
    # Retrieve full table if max_rows is not specified
    modseq_dat <- tbl(ch3_db$con, call_type) |> collect()
  }
  
  # Omit any missing values
  modseq_dat <- na.omit(modseq_dat)
  
  if (call_type == "regions") {
    # Aggregate mean_mh_frac by sample and region_name
    test_wide <- modseq_dat |>
      select(c(region_name, sample_name, mh_frac)) |>
      pivot_wider(names_from = sample_name, values_from = mh_frac) |>
      na.omit() |>
      as.data.frame()  # Convert to dataframe
  } else if (call_type == "windows") {
    # Aggregate mean_mh_frac by chr_pos and sample_name
    test_wide <- modseq_dat |>
      mutate(window = paste(chrom, start, end, sep = "_")) |>
      pivot_wider(id_cols = window, names_from = sample_name, values_from = mh_frac) |>
      na.omit() |>
      as.data.frame()  # Convert to dataframe
  } else {
    # Aggregate mean_mh_frac by chr_pos and sample_name
    test_wide <- modseq_dat |>
      mutate(chr_pos = paste(chrom, start, end, sep = "_")) |>
      pivot_wider(id_cols = chr_pos, names_from = sample_name, values_from = mh_frac) |>
      na.omit() |>
      as.data.frame()  # Convert to dataframe
  }
  
  # Ensure the collected data has the correct structure
  if (ncol(test_wide) <= 1) {
    stop("The data doesn't have enough columns for PCA after processing.")
  }
  
  # Scale the data: remove the first column (e.g., chr_pos or region_name)
  scaled_data <- scale(test_wide[, -1])
  
  # Perform PCA
  pca_result <- prcomp(t(scaled_data))  # Transpose data because samples should be rows
  pca_summary <- summary(pca_result)
  print(pca_summary)
  
  # Extract variance explained by PC1 and PC2
  pc1_var <- round(pca_summary$importance[2, 1] * 100, 2)
  pc2_var <- round(pca_summary$importance[2, 2] * 100, 2)
  
  # Prepare PCA data for plotting
  pca_data <- data.frame(pca_result$x)
  pca_data$sample_name <- rownames(pca_data)  # Add sample names
  
  print(pca_data)
  
  # Plot PCA
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_name)) +
          geom_point(size = 3) +
          theme_minimal() +
          labs(
            title = "PCA Plot of Methylation Data",
            x = paste0("PC1 (", pc1_var, "% variance)"),
            y = paste0("PC2 (", pc2_var, "% variance)")
          )
  print(p)
  # Save the plot if save_path is specified
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
    cat("PCA plot saved to ", save_path, "\n")
  }
  
  ch3_db <- .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}
