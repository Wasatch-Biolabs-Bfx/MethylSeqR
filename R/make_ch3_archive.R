#' Create a Compressed Methylation Call Archive (.ch3) from a calls file.
#'
#' Reads a delimited methylation call file and writes it to a compressed archive format.
#' Useful for standardizing and compressing methylation data across multiple samples.
#'
#' @param file_name Path to the input file (tab-delimited) containing methylation calls.
#' @param sample_name Name of the sample to embed in the output archive filenames.
#' @param out_path Output directory where the `.ch3` files will be written.
#' @param short_ids Logical; default is `TRUE`, shortens `read_id` by trimming prefixes to reduce file size.
#'
#' @return (Invisibly) the path template to the output archive files.
#' Also prints a message with timing information.
#' 
#' @details The input file should contain specific columns such as `read_id`, `chrom`, `ref_position`,
#' `ref_mod_strand`, `query_kmer`, `call_prob`, `call_code`, `base_qual`, and `flag`.
#' The function filters out unplaced '-' strand reads with position 0, adds new fields like `start` and `end`,
#' and writes the result as a compressed Arrow dataset with `zstd` compression.
#'
#' @importFrom dplyr filter mutate select if_else rename
#' @importFrom arrow write_dataset open_delim_dataset
#' @importFrom rlang .data
#' 
#' @examples
#' \dontrun{
#' make_ch3_archive("calls.tsv", "sample1", "output/", short_ids = TRUE)
#' }
#'
#' @export
make_ch3_archive <- function(file_name, 
                             sample_name,
                             out_path,
                             short_ids = TRUE) 
{
  start_time <- Sys.time()
  
  stopifnot("Invalid file_name" = 
              file.exists(file_name))
  
  # Read data as arrow table
  meth_data <- 
    open_delim_dataset(
      file_name, delim = "\t") |>
    rename(read_position = forward_read_position) |>
    select(
      read_id, chrom, read_position, ref_position, ref_mod_strand, read_length, query_kmer,
      call_prob, call_code, base_qual, flag) |>
    filter(!(ref_mod_strand == "-" & ref_position == 0)) |>
    mutate(
      sample_name = !!sample_name, 
      .before = read_id) |>
    mutate( # create shrot read IDs if TRUE, compresses file size.
      read_id = if_else(short_ids, gsub(".*-", "", read_id), read_id),
      .after = sample_name) |>
    mutate(
      start = if_else(ref_mod_strand == "-" & query_kmer == "CG", 
                      ref_position - 1,
                      ref_position),
      end = if_else(ref_mod_strand == "+" & query_kmer == "CG", 
                    ref_position + 2,
                    ref_position + 1),
      .after = chrom
    ) |>
    select(-ref_mod_strand) |>
    select(-ref_position) |>
    write_dataset(
      path = out_path, 
      basename_template = paste0(sample_name, "-{i}.ch3"),
      compression = "zstd"
    )
  
  elapsed_time <- round(as.numeric(Sys.time() - start_time, units = "secs"))
  
  # Report
  message("Wrote calls to the ", 
          sample_name, " file in ", out_path, 
          " (", elapsed_time, " secs)")
  
  # Return archive file name for piping/assignment only
  invisible(paste0(out_path, "/", sample_name, "-{i}.ch3"))
}