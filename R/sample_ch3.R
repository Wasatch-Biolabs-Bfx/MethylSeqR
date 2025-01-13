#' Downsample a CH3 File to a Desired Coverage
#'
#' This function downsamples a CH3 file to a specified coverage and writes the 
#' resulting data to a new CH3 file.
#'
#' @param ch3_file A string representing the file path to the input CH3 file.
#' @param out_path A string specifying the directory where the output CH3 file will be saved.
#' @param coverage A numeric value specifying the desired coverage for the downsampled CH3 file. Defaults to 5.
#' @param sample_name_suffix A string to append to the sample name in the output file. Defaults to an empty string.
#' @param genome_size A numeric value representing the size of the genome in base pairs. Defaults to 3.2e9 (human genome size).
#'
#' @details
#' This function calculates the current coverage of the input CH3 file based on the 
#' genome size and read lengths. If the current coverage is less than the desired coverage, 
#' the function copies the full file. Otherwise, it samples a fraction of reads to achieve 
#' the specified coverage. The sampled data is written to a new CH3 file with an updated 
#' sample name.
#'
#' The sample name in the output file is generated using the desired coverage, the 
#' original sample name extracted from the input file path, and an optional suffix.
#'
#' @return The function does not return a value but writes the downsampled CH3 file to the 
#' specified output path.
#'
#' @examples
#' # Downsample a CH3 file to 10x coverage and save it in the specified directory
#' 
#' # Set up the file path for the test data located in inst/test_data/
#' 
#' ch3_file_ex <- system.file("test_data/chr21_5xSample_BloodA-0.ch3", package = "MethylSeqR")
#' 
#' # Run the function with the example data
#' 
#' result = sample_ch3(ch3_file_ex, 
#'   out_path = system.file("downsampled, package = "MethylSeqR"), 
#'   coverage = 1, 
#'   sample_name_suffix = "_downsampled"
#' )
#'
#' @export

sample_ch3 <- function(ch3_file,
                       out_path,
                       coverage = 5, 
                       sample_name_suffix = "",
                       genome_size = 3.2e9)
{
  dat <- open_dataset(ch3_file)
  
  new_sample_name <- paste0(coverage, 
                            "xSample_", 
                            gsub(".*?/(.*)-..ch3", "\\1", ch3_file), 
                            sample_name_suffix)
  
  read_data <- 
    dat |>
    select(read_id, read_length) |>
    distinct() |>
    collect()
  
  full_coverage <- sum(read_data$read_length) / genome_size
  
  message("Original Coverage: ", full_coverage)
  
  fraction <- coverage / full_coverage
  
  if (fraction >= 1) {
    message("Warning: Original coverage less than desired coverage, copying full file.")
    dat |>
      mutate(sample_name = new_sample_name) |>
      write_dataset(
        path = out_path, 
        basename_template = paste0(new_sample_name, "-{i}.ch3"))
  } else {
    sampled_ids  <- sample(read_data$read_id, fraction * nrow(read_data), 
                           replace = FALSE)
    
    dat |>
      filter(
        read_id %in% sampled_ids) |>
      mutate(sample_name = new_sample_name) |>
      write_dataset(
        path = out_path, 
        basename_template = paste0(new_sample_name, "-{i}.ch3"))
  }
}