#' Perform PCA on Methylation Data
#'
#' This function performs Principal Component Analysis (PCA) on methylation (or other) data
#' retrieved from a DuckDB database. It aggregates the chosen value column based on the
#' specified call type and prepares it for PCA analysis.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param call_type A string representing the name of the table in the database from which to pull the data.
#'   Default is "positions".
#' @param value Column to use as the measurement for PCA (e.g., `mh_frac`, `m_frac`). Accepts a bare column
#'   name or a single string. Default: `m_frac`.
#' @param save_path Path to save the plot (e.g., .pdf or .png). If NULL, the plot is not saved.
#' @param max_rows Optional maximum number of rows to sample from the table (for speed on large datasets).
#'
#' @details
#' The function connects to the specified DuckDB database, retrieves data from `call_type`, reshapes
#' to a features × samples matrix (features are regions/windows/positions; columns are samples), scales
#' features, then runs PCA on samples (transpose before `prcomp`).
#'
#' @return Produces a PCA plot (PC1 vs PC2) and prints a PCA summary and the PCA scores.
#'
#' @examples
#' \dontrun{
#' # Default (m_frac)
#'  plot_ch3_pca(ch3_db)
#'  # Use mh_frac instead
#'  plot_ch3_pca(ch3_db, call_type = "regions", value = mh_frac)
#'  # Or as a string
#'  plot_ch3_pca(ch3_db, call_type = "windows", value = "mh_frac")
#' }
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl select mutate collect
#' @importFrom tidyr pivot_wider
#' @importFrom stats prcomp na.omit
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs ggsave
#' @importFrom rlang ensym as_name
#' @export
plot_ch3_pca <- function(ch3_db,
                         call_type = "positions",
                         value = m_frac,
                         save_path = NULL,
                         max_rows = NULL)
{
  start_time <- Sys.time()
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  # Resolve `value` (supports bare name or string)
  value_sym <- rlang::ensym(value)
  value_col <- rlang::as_name(value_sym)
  
  # If max_rows is specified, check table size and sample rows randomly in SQL
  if (!is.null(max_rows)) {
    row_count <- DBI::dbGetQuery(ch3_db$con, paste0("SELECT COUNT(*) as n FROM ", call_type))$n
    if (row_count < max_rows) {
      .ch3helper_closeDB(ch3_db)
      stop(paste0("Table '", call_type, "' only has ", row_count,
                  " rows, which is fewer than max_rows = ", max_rows, ". Pick fewer rows."))
    }
    # Use SQL random sampling with ORDER BY RANDOM()
    modseq_dat <- DBI::dbGetQuery(ch3_db$con, paste0("SELECT * FROM ", call_type, " ORDER BY RANDOM() LIMIT ", max_rows))
  } else {
    # Retrieve full table if max_rows is not specified
    modseq_dat <- dplyr::tbl(ch3_db$con, call_type) |> dplyr::collect()
  }
  
  # --- Check that the requested `value` column exists ---
  if (!(value_col %in% colnames(modseq_dat))) {
    .ch3helper_closeDB(ch3_db)
    stop(paste0(
      "The specified `value` column '", value_col, "' does not exist in the '", call_type, "' table.\n",
      "Available columns: ", paste(colnames(modseq_dat), collapse = ", "), "."
    ))
  }
  
  # Omit any rows with missing values (keeps PCA stable)
  modseq_dat <- stats::na.omit(modseq_dat)
  
  # Build a features × samples matrix depending on call_type
  if (call_type == "regions") {
    test_wide <- modseq_dat |>
      dplyr::select(region_name, sample_name, !!value_sym) |>
      tidyr::pivot_wider(names_from = sample_name, values_from = !!value_sym) |>
      stats::na.omit() |>
      as.data.frame()
  } else if (call_type == "windows") {
    test_wide <- modseq_dat |>
      dplyr::mutate(window = paste(chrom, start, end, sep = "_")) |>
      tidyr::pivot_wider(id_cols = window, names_from = sample_name, values_from = !!value_sym) |>
      stats::na.omit() |>
      as.data.frame()
  } else {
    test_wide <- modseq_dat |>
      dplyr::mutate(chr_pos = paste(chrom, start, end, sep = "_")) |>
      tidyr::pivot_wider(id_cols = chr_pos, names_from = sample_name, values_from = !!value_sym) |>
      stats::na.omit() |>
      as.data.frame()
  }
  
  # Ensure the collected data has the correct structure
  if (ncol(test_wide) <= 1) {
    .ch3helper_closeDB(ch3_db)
    stop("The data doesn't have enough columns for PCA after processing.")
  }
  
  # Scale features (drop the feature ID column)
  scaled_data <- scale(test_wide[, -1])
  
  # Perform PCA on samples (transpose so samples are rows)
  pca_result  <- stats::prcomp(t(scaled_data))
  pca_summary <- summary(pca_result)
  print(pca_summary)
  
  # Variance explained for PC1 and PC2
  pc1_var <- round(pca_summary$importance[2, 1] * 100, 2)
  pc2_var <- round(pca_summary$importance[2, 2] * 100, 2)
  
  # PCA scores for plotting
  pca_data <- data.frame(pca_result$x)
  pca_data$sample_name <- rownames(pca_data)
  print(pca_data)
  
  # Plot PCA
  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = sample_name)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PCA Plot of Methylation Data",
      x = paste0("PC1 (", pc1_var, "% variance)"),
      y = paste0("PC2 (", pc2_var, "% variance)")
    )
  print(p)
  
  # Save the plot if requested
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
    cat("PCA plot saved to ", save_path, "\n")
  }
  
  end_time <- Sys.time()
  message("Time elapsed: ", end_time - start_time, "\n")
  
  .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}
