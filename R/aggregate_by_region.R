#' Aggregate positional methylation data by regions provided with an annotation file.
#'
#' @param modseq_dat A duckdb object of methylation data already processed with summarize_by_pos().
#' @param annot_file Path to an annotation file of regions to aggregate data into- TSV format.
#' This file is required to have three columns- chrom, start, end. Optional fourth column of region_name can be included.
#' @return A duckdb object of all regional methylation data to be processed.
#' @examples
#' regional_data = aggregate_regions(data, "CpG_islands.tsv")
aggregate_regions <- function(modseq_dat,
                              annot_file)
{
  # Read annotation
  annotation <-
    read_tsv(annot_file,
             col_names = c("chrom", "start", "end", "region_name"),
             show_col_types = FALSE)

  # make sure column names did not get included and mess up code...
  if (annotation[1,1] %in% c("chr", "Chr", "chrom", "Chrom")) {
    annotation = annotation[-1, ]
  }

  stopifnot("Invalid annotation format. File must have chr, start, end." =
              ncol(annotation) >= 3)

  if (ncol(annotation == 3)) {
    annotation <-
      annotation |>
      mutate(
        region_name = paste(chrom, start, end, sep = "_"))
  }

  annotation <-
    annotation |>
    reframe(
      .by = c(region_name, chrom),
      ref_position = start:end)

  # Create regional dataframe
  modseq_dat |>
    print(head()) |>
    right_join(
      annotation,
      by = join_by(chrom, ref_position),
      copy = TRUE) |>
    summarize(
      .by = c(sample_name, region_name),
      cov = sum(cov),
      across(ends_with("_counts"), sum),
      across(ends_with("_frac"), ~ sum(.x * cov) / sum(cov)))
}
