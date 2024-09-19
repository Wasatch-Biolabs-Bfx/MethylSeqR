#' Cluster Methylation Data
#'
#' This function clusters methylation data using hierarchical clustering
#' and plots a dendrogram to visualize the clustering.
#'
#' @param modseq_dat A data frame containing methylation data. The data frame can
#' either contain positional or regional data. If it contains regional data, it must
#' have a column named \code{region_name}.
#'
#' @return None. The function plots a dendrogram.
#'
#' @examples
#' \dontrun{
#' cluster_modseq(modseq_dat)
#' }
#'
#' @import dplyr tidyr
#' @importFrom stats dist hclust as.dendrogram
#' @importFrom graphics plot
#'
#' @export
cluster_modseq <- function(ch3_db,
                           call_type = c("positions"))
{
  # If a character file name is provided, then make ch3 class obj
  db_con = helper_connectDB(ch3_db)
  
  # Check for specific table and connect to it in the database
  if (!dbExistsTable(db_con, call_type)) {
    stop(paste0(call_type, " Table does not exist. You can create it by..."))
  }
  
  modseq_dat = tbl(db_con, call_type) 

  if (call_type == "regions") {
    data_matrix <-
      na.omit(modseq_dat) |>
      dplyr::select(
        c(region_name, sample_name, mh_frac)) |>
      pivot_wider(
        names_from = sample_name,
        values_from = mh_frac) |>
      dplyr::select(
        -region_name) |>
      as.matrix()
  } else {
    data_matrix <-
      na.omit(modseq_dat) |>
      #slice_sample(
        #n = 50000000) |> # to make it faster with large positional data
      mutate(
        chr_pos = paste(chrom, ref_position, sep = "_")) |>
      pivot_wider(
        names_from = sample_name, values_from = mh_frac) |>
      dplyr::select(unique(modseq_dat$sample_name)) |>
      as.matrix()
  }

  # Compute the distance matrix
  distance_matrix <- dist(t(data_matrix))

  # Perform hierarchical clustering
  hc <- hclust(na.omit(distance_matrix), method = "complete")

  # Convert the clustering object to a dendrogram
  dend <- as.dendrogram(hc)

  # Plot the dendrogram using base R plot
  print(plot(dend, main = "Dendrogram of Methylation Samples"))
  
  # Finish up: close the connection
  dbDisconnect(db_con, shutdown = TRUE)
}
