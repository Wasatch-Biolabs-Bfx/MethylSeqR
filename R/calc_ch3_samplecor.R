#' Correlation Matrix of Modified Sequence Data
#'
#' This function calculates and optionally plots a correlation matrix for methylation or other
#' modification fraction data from genomic positions. It supports position-based, region-based,
#' and window-based calls and provides visualization using ggplot2.
#'
#' @param ch3_db A string. The path to the database containing ch3 files from nanopore data.
#' @param call_type A character vector specifying the type of data to retrieve from the database.
#'   Options are "positions", "regions", or "windows". Default is "positions".
#' @param value Column in the table to use as the measurement (e.g., mh_frac, m_frac).
#'   Accepts a bare column name or a string. Default: mh_frac.
#' @param agg_fun Function used to aggregate when multiple rows map to the same cell during pivot
#'   (e.g., duplicates in regions or windows). Default: mean.
#' @param plot Logical. If TRUE, plot a correlation heatmap. Default: TRUE.
#' @param save_path Optional file path to save the plot.
#' @param plot_sample_order Optional character vector setting sample order in the heatmap.
#' @param plot_title Title for the heatmap.
#' @param max_rows Optional integer to randomly sample rows from the table for faster runs.
#'
#' @details
#' The function connects to the ch3 database, pulls data for the selected `call_type`, reshapes
#' to a samples x samples matrix, and computes Pearson correlations using
#' `use = "pairwise.complete.obs"`. The `value` column controls which measurement is used.
#'
#' @return Invisibly returns the connection wrapper; prints the correlation matrix and (optionally)
#'   a ggplot heatmap.
#'
#' @examples
#' \dontrun{
#' # Use mh_frac (default)
#' calc_ch3_samplecor(ch3_db = "my_data.ch3.db", call_type = "positions")
#' # Use m_frac
#' calc_ch3_samplecor(ch3_db = "my_data.ch3.db", call_type = "regions", value = m_frac)
#' # Or as a string
#' calc_ch3_samplecor(ch3_db = "my_data.ch3.db", call_type = "windows", value = "m_frac")
#' }
#'
#' @importFrom DBI dbConnect dbDisconnect dbExistsTable dbGetQuery
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl select distinct arrange left_join mutate pull collect across where
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 scale_y_discrete theme_minimal labs theme element_text ggsave
#' @importFrom rlang enquo quo_is_symbol as_name sym is_string
#' @importFrom stats cor
#' @export
calc_ch3_samplecor <- function(ch3_db,
                               call_type = c("positions"),
                               value = m_frac,
                               agg_fun = mean,
                               plot = TRUE,
                               save_path = NULL,
                               plot_sample_order = NULL,
                               plot_title = "Sample Correlation Matrix",
                               max_rows = NULL)
{
  start_time <- Sys.time()
  
  # Tidy-eval: allow bare name or string for `value`
  value_quo <- rlang::enquo(value)
  if (rlang::quo_is_symbol(value_quo)) {
    value_col <- rlang::as_name(value_quo)
  } else {
    # If passed as string like value = "m_frac"
    value_str <- tryCatch(eval.parent(substitute(value)), error = function(e) NULL)
    if (is.character(value_str) && length(value_str) == 1) {
      value_col <- value_str
      value_quo <- rlang::sym(value_col)
    } else {
      stop("`value` must be a column name (bare) or a single string.")
    }
  }
  
  # Open the database connection
  ch3_db <- .ch3helper_connectDB(ch3_db)
  
  if (length(call_type) > 1) {
    call_type <- "positions"
  }
  
  if (!dbExistsTable(ch3_db$con, call_type)) {
    stop(paste0(call_type, " table does not exist in the database."))
  }
  
  # Retrieve data (optionally sampled)
  if (!is.null(max_rows)) {
    row_count <- dbGetQuery(ch3_db$con, paste0("SELECT COUNT(*) as n FROM ", call_type))$n
    if (row_count < max_rows) {
      stop(paste0("Table '", call_type, "' only has ", row_count,
                  " rows, which is fewer than max_rows = ", max_rows, ". Pick fewer rows."))
    }
    modseq_dat <- dbGetQuery(ch3_db$con, paste0("SELECT * FROM ", call_type, " ORDER BY RANDOM() LIMIT ", max_rows))
  } else {
    modseq_dat <- dplyr::tbl(ch3_db$con, call_type) |> dplyr::collect()
  }
  
  # --- New check for value column existence ---
  if (!(value_col %in% colnames(modseq_dat))) {
    .ch3helper_closeDB(ch3_db)
    stop(paste0("The specified `value` column '", value_col,
                "' does not exist in the '", call_type, "' table.\n",
                "Available columns: ", paste(colnames(modseq_dat), collapse = ", "), "."))
  }
  
  
  # Quick column presence check for the selected value column
  if (!(value_col %in% names(modseq_dat))) {
    on.exit(.ch3helper_closeDB(ch3_db), add = TRUE)
    stop("Selected `value` column '", value_col, "' is not present in the '", call_type, "' table.")
  }
  
  # Build a named values_fn for pivot_wider when needed
  values_fn_named <- setNames(list(agg_fun), value_col)
  
  if (call_type == "regions") {
    # Aggregate value by sample and region_name
    dat_wide <- modseq_dat |>
      tidyr::pivot_wider(
        id_cols   = region_name,
        names_from = sample_name,
        values_from = !!value_quo,
        values_fn  = values_fn_named
      )
    
    numeric_columns <- dat_wide[, unique(modseq_dat$sample_name), drop = FALSE]
    
    correlation_matrix <- stats::cor(numeric_columns,
                                     use = "pairwise.complete.obs",
                                     method = "pearson")
    
  } else if (call_type == "windows") {
    
    common_windows <- modseq_dat |>
      dplyr::select(chrom, start, end) |>
      dplyr::distinct() |>
      dplyr::arrange(chrom, start, end)
    
    aligned_data <- common_windows |>
      dplyr::left_join(modseq_dat, by = c("chrom", "start", "end")) |>
      tidyr::pivot_wider(
        id_cols    = c(chrom, start, end),
        names_from = sample_name,
        values_from = !!value_quo,
        values_fn   = values_fn_named
      )
    
    # Replace NA with 0 *only* for the measurement columns (keep coords intact)
    coord_cols <- c("chrom", "start", "end")
    measure_cols <- setdiff(names(aligned_data), coord_cols)
    aligned_data[measure_cols] <- lapply(aligned_data[measure_cols], function(x) { x[is.na(x)] <- 0; x })
    
    numeric_columns <- aligned_data |> dplyr::select(-chrom, -start, -end)
    
    correlation_matrix <- stats::cor(numeric_columns,
                                     use = "pairwise.complete.obs",
                                     method = "pearson")
    
  } else {
    # positions
    dat_wide <- modseq_dat |>
      dplyr::mutate(chr_pos = paste(chrom, start, end, sep = "_")) |>
      tidyr::pivot_wider(
        id_cols    = chr_pos,
        names_from = sample_name,
        values_from = !!value_quo,
        values_fn   = values_fn_named
      )
    
    # Try to ensure numerics (without breaking factors/characters that slipped in)
    numeric_columns <- dat_wide |>
      dplyr::select(-chr_pos) |>
      dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character)) |>
      dplyr::mutate(dplyr::across(dplyr::where(is.character), suppressWarnings(as.numeric)))
    
    if (!all(vapply(numeric_columns, is.numeric, logical(1)))) {
      stop("Some sample columns could not be converted to numeric. Check the input data.")
    }
    
    correlation_matrix <- stats::cor(numeric_columns,
                                     use = "pairwise.complete.obs",
                                     method = "pearson")
  }
  
  print(correlation_matrix)
  
  if (plot) {
    melted_cor <- as.data.frame(correlation_matrix) |>
      dplyr::mutate(Var1 = rownames(correlation_matrix)) |>
      tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value")
    
    if (!is.null(plot_sample_order)) {
      melted_cor$Var1 <- factor(melted_cor$Var1, levels = plot_sample_order)
      melted_cor$Var2 <- factor(melted_cor$Var2, levels = plot_sample_order)
    }
    
    p <- ggplot2::ggplot(melted_cor, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), color = "black") +
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                    midpoint = 0, limits = c(-1, 1),
                                    space = "Lab", name = "Pearson\nCorrelation") +
      ggplot2::scale_y_discrete(limits = rev(levels(factor(melted_cor$Var2)))) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = plot_title, x = "Sample", y = "Sample") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
    print(p)
    
    end_time <- Sys.time()
    message("Correlation analysis complete - Time elapsed: ", end_time - start_time, "\n")
    
    if (!is.null(save_path)) {
      ggplot2::ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
      cat("Correlation plot saved to ", save_path, "\n")
    }
  }
  
  .ch3helper_closeDB(ch3_db)
  invisible(ch3_db)
}
