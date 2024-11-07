#' Create CH3 Archive from Methylation Data
#'
#' This function reads methylation data from a specified file and creates an archive of the data 
#' in CH3 format, storing it in a designated output path. It processes the data by selecting relevant 
#' columns and mutating the reference position based on the modification strand.
#'
#' @param file_name A string representing the path to the input file containing methylation data.
#' @param sample_name A string representing the name of the sample. This name will be included in 
#' the output files.
#' @param out_path A string representing the path where the output CH3 dataset will be stored.
#'
#' @details
#' The function uses the Arrow package to read the input TSV dataset and select specific columns. 
#' It then adjusts the reference position based on the modification strand before writing the dataset 
#' to the specified output path. The output files are named using the provided `sample_name` and a 
#' unique index.
#'
#' @return A string representing the path of the created archive file, useful for piping or assignment.
#' The function is designed to be used as a part of a data processing pipeline.
#'
#' @import arrow
#' @import dplyr
#' @import dbplyr
#' @import tidyr
#'
#' @export
make_ch3_archive <- function(file_name, 
                             sample_name,
                             out_path) 
{
  start_time <- Sys.time()

  stopifnot("Invalid file_name" = 
            file.exists(file_name))

  # Read data as arrow table
  meth_data <- 
    open_tsv_dataset(
      file_name) |>
    select(
      read_id, chrom, ref_position, ref_mod_strand, read_length,
      call_prob, call_code, base_qual) |>
    mutate(
      sample_name = sample_name, 
      .before = read_id) |>
    mutate(
      ref_position = ifelse(ref_mod_strand == "-", 
                          ref_position - 1, ref_position)) |>
    write_dataset(
      path = out_path, 
      basename_template = paste0(sample_name, "-{i}.ch3"))

    elapsed_time <- round(as.numeric(Sys.time() - start_time, units = "secs"))

  # Report
  message("Wrote ", nrow(meth_data), " calls to the ", 
         sample_name, " table in ", out_path, 
          " (", elapsed_time, " secs)")

  # Return archive file name for piping/assignment only
  invisible(paste0(out_path, "/", sample_name, "-{i}.ch3"))
  }