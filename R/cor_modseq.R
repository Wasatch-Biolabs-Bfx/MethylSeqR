#' Correlation Matrix of Modified Sequence Data
#'
#' This function calculates and optionally plots a correlation matrix for methylation or other modification fraction data 
#' from genomic positions. It can handle both position-based and region-based calls and supports visualization using ggplot2.
#'
#' @param ch3_db A string. The path to the database containing ch3 files from nanopore data.
#' @param call_type A character vector specifying the type of data to retrieve from the database. Default is "positions". Can alo be "regions".
#' @param plot A logical value. If TRUE, a correlation heatmap with correlation values is plotted. Default is FALSE.
#' @param save_path Pathway to save the plot to. Usually .pdf or .png.
#'
#' @details
#' This function connects to the ch3 files database, retrieves either position or region-based modification fraction data, 
#' and computes the Pearson correlation matrix for each sample. If `call_type` is "regions", it aggregates data by sample 
#' and region name before computing the correlation matrix. If `plot = TRUE`, the function generates a heatmap using 
#' ggplot2 with the correlation values displayed.
#'
#' @note The function assumes that the database has tables named according to the `call_type` parameter (e.g., "positions", 
#' "regions"). The correlation matrix is calculated using the `pairwise.complete.obs` method, which handles missing data.
#'
#' @importFrom stats cor
#' @import DBI tidyr ggplot2
#' 
#' @return A correlation matrix of modification fractions across samples. If `plot = TRUE`, a ggplot object of the correlation 
#' matrix heatmap is also printed.
#'
#' @examples
#' cor_modseq(ch3_db = "path/to/database.db", call_type = "positions", plot = TRUE)
#'
#' @export

cor_modseq <- function(ch3_db,
                       call_type = c("positions"),
                       plot = FALSE)
{
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
  
  if (length(call_type) > 1) {
    call_type = c("positions")
  }
  
  tryCatch(
    {
      # Check if the call_type table exists in the database
      if (!dbExistsTable(db_con, call_type)) {
        stop(paste0(call_type, " Table does not exist. You can create it by..."))
      }
      
      # Retrieve the 'positions' table
      modseq_dat <- tbl(db_con, call_type) %>% collect()
      
      if (call_type == "regions") {
        print("regional data")
        # Aggregate mean_mh_frac by sample and region_name
        dat_wide <- modseq_dat |>
          pivot_wider(id_cols = region_name,
                      names_from = sample_name,
                      values_from = mh_frac)
        # Compute Correlation
        numeric_columns <- dat_wide[, unique(modseq_dat$sample_name)]
        # Calculate correlation matrix
        correlation_matrix <- cor(numeric_columns,
                                  use = "pairwise.complete.obs",
                                  method = "pearson")
        # View the correlation matrix
        print(correlation_matrix)
      } else {
        # Create a 'chr_pos' identifier for each genomic position
        dat_wide <- modseq_dat %>%
          mutate(chr_pos = paste(chrom, ref_position, sep = "_")) %>%
          pivot_wider(id_cols = chr_pos,
                      names_from = sample_name,  # Each sample (barcode) becomes a column
                      values_from = mh_frac)     # Values are the 'mh_frac' for each sample
        
        # Convert the sample columns to numeric (if needed)
        numeric_columns <- dat_wide %>%
          select(-chr_pos) %>%
          mutate(across(everything(), as.numeric))  # Convert all columns to numeric
        
        # Check if conversion was successful
        if (!all(sapply(numeric_columns, is.numeric))) {
          stop("Some columns could not be converted to numeric.")
        }
        
        # Calculate the correlation matrix
        correlation_matrix <- cor(numeric_columns, use = "pairwise.complete.obs", method = "pearson")
        # Print the correlation matrix
        print(correlation_matrix)
      }
      
      
      # Plot Correlation Matrix with Correlation Values
      if (plot) {
        melted_cor <- as.data.frame(correlation_matrix) %>%
          mutate(Var1 = rownames(correlation_matrix)) %>%  # Add row names as a new column
          pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value")  # Pivot to long format
        
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
                labs(title = "Sample Correlation Matrix",
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
      
    }, 
    error = function(e)
    {
      # Print custom error message
      message("An error occurred: ", e$message)
      # Optionally, re-throw the error if needed
      # stop(e)
    },
    finally = 
      {
        # Finish up: close the connection
        .helper_closeDB(database)
      }
  )
}
