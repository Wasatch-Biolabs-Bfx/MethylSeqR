.compute_windows_chr <- function(cur_chrom,
                                 db_con,
                                 modseq_dat,
                                 window_size)
{
  chrom_dat <- 
    modseq_dat |>
    filter(chrom == cur_chrom)
  
  positions <- tibble(ref_position = 1:max(pull(chrom_dat, ref_position)))
  
  # Join tables
  dat <- 
    chrom_dat |>
    right_join(
      positions, 
      by = join_by(ref_position),
      copy = TRUE) |>
    tidier::mutate(
      across(c(cov, ends_with("_counts")), sum),
      across(ends_with("_frac"), ~ sum(.x * cov, na.rm = TRUE) / 
                                   sum(cov, na.rm = TRUE)),
      .by = sample_name,
      .order_by = ref_position,
      .frame = c(window_size, 0)) |>
    filter(!is.na(sample_name)) |>
    rename(end = ref_position) |>
    mutate(
      .before = end,
      start = end - window_size + 1) |>
    collect()
  
  message("Finished ", cur_chrom, ": ", nrow(dat), " windows have data.")
  
  if (nrow(dat) == 0) return()
  
  if (dbExistsTable(db_con, "sliding_windows")) {
    dbAppendTable(db_con, "sliding_windows", dat)
  } else {
    copy_to(dest = db_con, name = "sliding_windows", df = dat)
  }
}

compute_sliding_windows <- function(modseq_dat,
                                    window_size = 500)
{
  # Setup DB
  db_con <- dbConnect(duckdb(tempfile()))
  
  if(dbExistsTable(db_con, "sliding_windows")) {
    dbRemoveTable(db_con, "sliding_windows")
  }
  
  # Get vector of chromosome names
  chrs <- 
    modseq_dat |>
    select(chrom) |>
    distinct() |>
    arrange(chrom) |>
    pull()
  
  return_vals <- lapply(chrs, .compute_windows_chr, 
                        db_con = db_con, 
                        modseq_dat = modseq_dat, 
                        window_size = window_size)
  
  tbl(db_con, "sliding_windows")
}
