#' Correlation Matrix of Modified Sequence Data
#'
#' This function calculates and optionally plots a correlation matrix for methylation or other modification fraction data 
#' from genomic positions. It supports position-based, region-based, and window-based calls and provides visualization 
#' using ggplot2.
#'
#' @param ch3_db A string. The path to the database containing ch3 files from nanopore data.
#' @param call_type A character vector specifying the type of data to retrieve from the database. Options are "positions", 
#'   "regions", or "windows". Default is "positions".
#' @param plot A logical value. If TRUE, a correlation heatmap with correlation values is plotted. Default is FALSE.
#' @param save_path A string. The file path to save the plot (e.g., .pdf or .png). If NULL, the plot is not saved. Default is NULL.
#' @param plot_sample_order A character vector specifying the desired order of samples in the plot. Default is NULL, 
#'   which uses the default order.
#' @param plot_title A string. The title of the correlation heatmap. Default is "Sample Correlation Matrix".
#' @param max_rows The maximum amount of rows wanted for calculation. This argument can help analysis run faster when there is a lot of data.
#'
#' @details
#' This function connects to the ch3 files database, retrieves data based on the `call_type` parameter, and computes the 
#' Pearson correlation matrix for each sample. For "regions", data is aggregated by sample and region name. For "windows", 
#' data is aligned to common genomic windows before computing the correlation matrix. If `plot = TRUE`, the function 
#' generates a ggplot2 heatmap with correlation values displayed. Optionally, the plot can be saved to a specified file path.
#'
#' @note The function assumes that the database contains tables named according to the `call_type` parameter (e.g., "positions", 
#' "regions", "windows"). It calculates the correlation matrix using the `pairwise.complete.obs` method, which handles missing data.
#' 
#' @return A correlation matrix of modification fractions across samples. If `plot = TRUE`, a ggplot object of the correlation 
#' matrix heatmap is also printed. If `save_path` is specified, the plot is saved to the given file path.
#'
#' @examples
#' # Specify the path to the database
#' ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#' 
#' # Run the correlation matrix function using the 'positions' call type and plot the results
#' calc_ch3_samplecor(ch3_db = ch3_db, call_type = "positions")
#'
#' @importFrom DBI dbConnect dbDisconnect dbExistsTable dbGetQuery
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl select distinct arrange left_join mutate pivot_wider pull collect
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 scale_y_discrete theme_minimal labs theme element_text ggsave
#' @importFrom stats cor
#' 
#' @export

calc_ch3_samplecor <- function(ch3_db,
                       call_type = c("positions"),
                       plot = TRUE,
                       save_path = NULL,
                       plot_sample_order = NULL,
                       plot_title = "Sample Correlation Matrix",
                       max_rows = NULL)
{
  # Open the database connection
  database <- .ch3helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  # Specify on exit what to do...
  # Finish up: update table list and close the connection
  on.exit(.ch3helper_closeDB(database), add = TRUE)
  
  if (length(call_type) > 1) {
    call_type = c("positions")
  }
  
  # Check if the call_type table exists in the database
  if (!dbExistsTable(db_con, call_type)) {
    stop(paste0(call_type, " Table does not exist. You can create it by..."))
  }
  
  # If max_rows is specified, check table size and sample rows randomly in SQL
  if (!is.null(max_rows)) {
    row_count <- dbGetQuery(db_con, paste0("SELECT COUNT(*) as n FROM ", call_type))$n
    
    if (row_count < max_rows) {
      stop(paste0("Table '", call_type, "' only has ", row_count, " rows, which is fewer than max_rows = ", max_rows, ". Pick less rows."))
    }
    
    # Use SQL random sampling with ORDER BY RANDOM()
    modseq_dat <- dbGetQuery(db_con, paste0("SELECT * FROM ", call_type, " ORDER BY RANDOM() LIMIT ", max_rows))
  } else {
    # Retrieve full table if max_rows is not specified
    modseq_dat <- tbl(db_con, call_type) |> collect()
  }
  
  # Retrieve the 'positions' table
  # modseq_dat <- tbl(db_con, call_type) |> collect()
  
  if (call_type == "regions") {
    print("regional data")
    # Aggregate mean_mh_frac by sample and region_name
    dat_wide <- modseq_dat |>
      pivot_wider(id_cols = region_name,
                  names_from = sample_name,
                  values_from = mh_frac,
                  values_fn = mean)
    
    # Compute Correlation
    numeric_columns <- dat_wide[, unique(modseq_dat$sample_name)]
    
    print(dim(numeric_columns))
    # Calculate correlation matrix
    correlation_matrix <- cor(numeric_columns,
                              use = "pairwise.complete.obs",
                              method = "pearson")
    # View the correlation matrix
    print(correlation_matrix)
  } else if (call_type == "windows") {
    
    # Step 1: Define common windows
    common_windows <- modseq_dat |>
      select(chrom, start, end) |>
      distinct() |>
      arrange(chrom, start, end)
    
    # Step 2: Align data to common windows
    aligned_data <- common_windows |>
      left_join(modseq_dat, by = c("chrom", "start", "end")) |>
      pivot_wider(
        id_cols = c(chrom, start, end),
        names_from = sample_name,
        values_from = mh_frac
      )
    
    # Step 3: Replace NA values with 0
    aligned_data[is.na(aligned_data)] <- 0
    
    # Step 4: Compute correlation matrix
    numeric_columns <- aligned_data |> select(-chrom, -start, -end)
    correlation_matrix <- cor(numeric_columns, 
                              use = "pairwise.complete.obs", 
                              method = "pearson")
    
    # Print the correlation matrix
    print(correlation_matrix)
  }
  else {
    # Create a 'chr_pos' identifier for each genomic position
    dat_wide <- modseq_dat |>
      mutate(chr_pos = paste(chrom, start, end, sep = "_")) |>
      pivot_wider(id_cols = chr_pos,
                  names_from = sample_name,  # Each sample (barcode) becomes a column
                  values_from = mh_frac)     # Values are the 'mh_frac' for each sample
    
    # Convert the sample columns to numeric (if needed)
    numeric_columns <- dat_wide |>
      select(-chr_pos) |>
      mutate(across(everything(), as.numeric))  # Convert all columns to numeric
    
    # Check if conversion was successful
    if (!all(sapply(numeric_columns, is.numeric))) {
      stop("Some columns could not be converted to numeric.")
    }
    
    # Calculate the correlation matrix
    correlation_matrix <- cor(numeric_columns, 
                              use = "pairwise.complete.obs", 
                              method = "pearson")
    # Print the correlation matrix
    print(correlation_matrix)
  }
  
  
  # Plot Correlation Matrix with Correlation Values
  if (plot) {
    melted_cor <- as.data.frame(correlation_matrix) |>
      mutate(Var1 = rownames(correlation_matrix)) |>  # Add row names as a new column
      pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value")  
    # Pivot to long format
    
    # Define the desired order for samples
    if (!is.null(plot_sample_order)) {
      sample_order <- plot_sample_order
      
      # Reorder the factors for Var1 and Var2
      melted_cor$Var1 <- factor(melted_cor$Var1, levels = sample_order)
      melted_cor$Var2 <- factor(melted_cor$Var2, levels = sample_order)
    }
    
    
    p <- ggplot(data = melted_cor, aes(x = Var1,
                                       y = Var2,
                                       fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 2)),
                color = "black") + # Add correlation values
      scale_fill_gradient2(low = "blue",
                           high = "red",
                           mid = "white",
                           midpoint = 0,
                           limits = c(-1, 1),
                           space = "Lab",
                           name = "Pearson\nCorrelation") +
      scale_y_discrete(limits = rev(levels(factor(melted_cor$Var2)))) +
      theme_minimal() +
      labs(title = plot_title,
           x = "Sample",
           y = "Sample") +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       hjust = 1))
    print(p)
    
    # Save the plot if save_path is specified
    if (!is.null(save_path)) {
      ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
      cat("Correlation plot saved to ", save_path, "\n")
    }
  }
      
}