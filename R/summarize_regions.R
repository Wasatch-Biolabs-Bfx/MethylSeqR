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
  
  # left join- keep reads that are outside of the annotation table
  
  # right join- keep reads in annotation table + regions not included in data frame
  
  # inner join- regions in both annotation and data frame...
  my_join <- switch(join_type,
                    "inner" = inner_join,
                    "right" = right_join,
                    "left" = left_join,
                    "full" = full_join)
  
  if (dbExistsTable(db_con, "regions"))
    dbRemoveTable(db_con, "regions")
  
  # Create regions table by unique chrs
  chroms <-
    annotation |>
    select(chrom) |>
    distinct() |>
    filter(nchar(chrom) < 6) |>
    pull()
  
  print(chroms)
  
  # Create Progress Bar
  pb <- progress_bar$new(
    format = "[:bar] :percent [Elapsed time: :elapsed]",
    total = length(chroms) + 1,
    complete = "=",
    incomplete = "-",
    current = ">",
    clear = FALSE,
    width = 100)
  
  pb$tick()
  
  for (chr in chroms) {
    # Begin summarizing by region- perform the join and aggregation
    tbl(db_con, "positions") |>
      filter(chrom == chr) |>
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
    
    # Create or append table
    dbExecute(db_con, 
              "CREATE TABLE IF NOT EXISTS regions AS 
              SELECT * FROM temp_table WHERE 1=0")
    
    dbExecute(db_con, 
              "INSERT INTO regions SELECT * FROM temp_table")
    
    dbRemoveTable(db_con, "temp_table")
    
    pb$tick()
  }
  
  # Close progress bar
  pb$terminate()
  
  # Finish up: update table list and close the connection
  ch3_db$tables <- dbListTables(db_con)
  dbDisconnect(db_con, shutdown = TRUE)
  
  return(ch3_db)
}

# library(readr)
# library(tidyr)
# library(dplyr)
# library(dbplyr)
# library(duckdb)
# library(duckplyr)
# 
# # Function
# summarize_regions <- function(ch3_db,
#                               region_bed,
#                               join_type = "inner")
# {
#   # Read annotation
#   annotation <-
#     read_csv(region_bed,
#              col_names = c("chrom", "start", "end", "region_name"),
#              show_col_types = FALSE)
#   
#   # make sure column names did not get included and mess up code...
#   if (annotation[1,1] %in% c("chr", "Chr", "chrom", "Chrom")) {
#     annotation = annotation[-1, ]
#   }
#   if (ncol(annotation) < 3 || ncol(annotation) > 4) {
#     stop("Invalid annotation format. File must have 3 or 4 columns:
#          chr, start, end, region_name (optional) annotation.")
#   }
#   if (!join_type %in% c("inner", "right", "left", "full")) {
#     stop("Invalid join type")
#   }
#   if (ncol(annotation == 3)) {
#     annotation <-
#       annotation |>
#       mutate(
#         region_name = paste(chrom, start, end, sep = "_"))
#   }
#   
#   annotation <-
#     annotation |>
#     reframe(
#       .by = c(region_name, chrom),
#       ref_position = start:end)
#   
#   # Open the database connection
#   if (inherits(ch3_db, "character")) {
#     ch3_db <- list(db_file = ch3_db)
#     class(ch3_db) <- "ch3_db"
#     db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
#     ch3_db$tables <- dbListTables(db_con)
#   } else if (inherits(ch3_db, "ch3_db")) {
#     db_con <- dbConnect(duckdb(ch3_db$db_file), read_only = FALSE)
#   } else {
#     stop("Invalid ch3_db class. Must be character or ch3_db.")
#   }
#   
#   # Create regional data frame- offer a left, right or inner join
#   
#   # left join- keep reads that are outside of the annotation table
#   
#   # right join- keep reads in annotation table + regions not included in data frame
#   
#   # inner join- regions in both annotation and data frame...
#   my_join <- switch(join_type,
#                     "inner" = inner_join,
#                     "right" = right_join,
#                     "left" = left_join,
#                     "full" = full_join)
#   
#   # Create regions table by unique chrs
#   # chroms <-
#   #   annotation |>
#   #   select(chrom) |>
#   #   distinct() |>
#   #   pull()
#   
#   # Create Progress Bar
#   # pb <- progress_bar$new(
#   #   format = "[:bar] :percent [Elapsed time: :elapsed]",
#   #   total = length(chroms) + 1,
#   #   complete = "=",
#   #   incomplete = "-",
#   #   current = ">",
#   #   clear = FALSE,
#   #   width = 100)
#   # 
#   # pb$tick()
#   # 
#   # for (chr in chroms) {
#     # Begin summarizing by region- perform the join and aggregation
#   if (dbExistsTable(db_con, "regions"))
#     dbRemoveTable(db_con, "regions")
#   
#     regions = tbl(db_con, "positions") |>
#       my_join(annotation, 
#               by = join_by(chrom, 
#                            ref_position), 
#               copy = TRUE) |>
#       summarize(
#         .by = c(sample_name, region_name),
#         cov = sum(cov),
#         across(ends_with("_counts"), sum),
#         across(ends_with("_frac"), ~ sum(.x * cov) / sum(cov))) |>
#       compute()
#     
#     dbWriteTable(
#       conn = db_con, 
#       name = "regions", 
#       value = regions
#       append = TRUE)
#     
#     
#     # compute(name = "regions", temporary = FALSE) |>
#     
#     # Create or append table
#     # dbExecute(db_con, 
#     #           "CREATE TABLE IF NOT EXISTS regions AS 
#     #           SELECT * FROM temp_table WHERE 1=0")
#     # 
#     # dbExecute(db_con, 
#     #           "INSERT INTO regions SELECT * FROM temp_table")
#     # 
#     # dbRemoveTable(db_con, "temp_table")
#     # 
#     # pb$tick()
#   # }
#   
#   # Close progress bar
#   # pb$terminate()
#   
#   # Finish up: update table list and close the connection
#   ch3_db$tables <- dbListTables(db_con)
#   dbDisconnect(db_con, shutdown = TRUE)
#   
#   return(ch3_db)
# }