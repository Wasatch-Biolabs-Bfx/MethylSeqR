open_dat <- function(path,
                     samples = "all",
                     chrs = "all",
                     min_read_length = 100,
                     min_call_prob = .9,
                     min_base_qual = 10,
                     max_memory = "2GB",
                     max_threads = 4)
{
  # Perform checks and calcs for arguments
  stopifnot(
    "A single folder must be provided" =
      length(path) == 1)

  if (all(samples == "all"))
    samples <-
      list.dirs(path = path,
                full.names = FALSE, recursive = FALSE)

  if (all(chrs == "all"))
    chrs <-
      list.dirs(path = paste0(path, "/", samples),
                full.names = FALSE, recursive = FALSE) |>
      unique()

  stopifnot(
    "Could not find one or more samples in your path" =
      dir.exists(paste0(path, "/", samples)))

  # Set up db connection
  db_con <- dbConnect(duckdb(tempfile()))

  # Set memory and processor limits
  dbExecute(
    db_con,
    sprintf("SET memory_limit = '%s';
             SET threads TO %i;",
            max_memory, max_threads))

  # open (and return) data
  open_dataset(
    path,
    partitioning = c("sample_name", "chrom")) |>
    filter(
      sample_name %in% samples,
      chrom %in% chrs,
      call_prob >= min_call_prob,
      read_length >= min_read_length,
      base_qual >= min_base_qual) |>
    to_duckdb(con = db_con) 
}
