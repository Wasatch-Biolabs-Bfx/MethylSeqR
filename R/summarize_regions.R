library(readr)
library(tidyr)
library(dplyr)
library(dbplyr)
library(duckdb)
library(duckplyr)

# Function
summarize_regions <- function(ch3_db,
                              region_bed,
                              join_type = "inner")
{
  # Read annotation
  annotation <-
    read_csv(region_bed,
             col_names = c("chrom", "start", "end", "region_name"),
             show_col_types = FALSE)
  
  # make sure column names did not get included and mess up code...
  if (annotation[1,1] %in% c("chr", "Chr", "chrom", "Chrom")) {
    annotation = annotation[-1, ]
  }
  if (ncol(annotation) < 3 || ncol(annotation) > 4) {
    stop("Invalid annotation format. File must have 3 or 4 columns:
         chr, start, end, region_name (optional) annotation.")
  }
  if (!join_type %in% c("inner", "right", "left", "full")) {
    stop("Invalid join type")
  }
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
  
  # Open the database connection
  if (inherits(ch3_db, "character")) {
    ch3_db <- list(db_file = ch3_db)
    class(ch3_db) <- "ch3_db"
    db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
    ch3_db$tables <- dbListTables(db_con)
  } else if (inherits(ch3_db, "ch3_db")) {
    db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
  } else {
    stop("Invalid ch3_db class. Must be character or ch3_db.")
  }
  
  # Create regional data frame- offer a left, right or inner join
  
  # left join, if they wanna keep reads that are outside of the annotation table
  
  # right join, if they want to keep reads that are all in annotation table and
  # even those regions not included in data frame
  
  # inner join- regions in both annotation and data frame...
  my_join <- switch(join_type,
                    "inner" = inner_join,
                    "right" = right_join,
                    "left" = left_join,
                    "full" = full_join)
  
  # Initialize progress bar (for the number of regions to be processed)
  # pb <- progress_bar$new(
  #   format = "[:bar] :percent [Elapsed time: :elapsed]",
  #   total = nrow(annotation),  # Progress based on number of annotation regions
  #   complete = "=",   
  #   incomplete = "-", 
  #   current = ">",    
  #   clear = FALSE,    
  #   width = 100
  # )
  
  # make a loop through chr
  # "diff_meth" 
  
  # Begin summarizing by region
  # Perform the join and aggregation
  regional_data <- tbl(db_con, "positions") |>
    my_join(annotation, 
            by = join_by(chrom, 
                         ref_position), 
            copy = TRUE) |>
    summarize(
      .by = c(sample_name, region_name),
      cov = sum(cov),
      across(ends_with("_counts"), sum),
      across(ends_with("_frac"), ~ sum(.x * cov) / sum(cov))) |>
    compute(name = "temp_table", temporary = TRUE)
  
  # Update progress bar for each iteration (in this case, for each region)
  # pb$tick()
  
  # Write the aggregated data to the database (if needed)
  dbWriteTable(conn = db_con, name = "regions", value = regional_data, overwrite = TRUE)
  
  # Finish up: update table list and close the connection
  ch3_db$tables <- dbListTables(db_con)
  dbDisconnect(db_con, shutdown = TRUE)
  
  return(ch3_db)
}