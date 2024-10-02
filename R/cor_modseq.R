library(tidyr)
library(reshape2)

cor_modseq <- function(ch3_db,
                       call_type = c("positions"),
                       plot = FALSE)
{
  # Connect to the database
  db_con <- helper_connectDB(ch3_db)
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
          pivot_wider(names_from = sample_name,
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
        # melted_cor <- melt(correlation_matrix)
        # Replace `melt(correlation_matrix)` with `pivot_longer()`
        melted_cor <- as.data.frame(correlation_matrix) %>%
          mutate(Var1 = rownames(correlation_matrix)) %>%  # Add row names as a new column
          pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value")  # Pivot to long format
        
        print(ggplot(data = melted_cor, aes(x = Var1,
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
                                                 hjust = 1)))
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
        dbDisconnect(db_con, shutdown = TRUE)
      }
  )
}
