## Make Archive

# Dependencies
library(arrow)
library(tidyr)
library(dplyr)
library(dbplyr)

# Function
make_ch3_archive <- function(file_name, 
                             sample_name,
                             out_path) 
{
  start_time <- Sys.time()

  stopifnot("Invalid file_name" = 
            !file.exists(file_name))

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