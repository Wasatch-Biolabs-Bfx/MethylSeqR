## ----include=FALSE------------------------------------------------------------
library(knitr, quietly=TRUE)
library(MethylSeqR, quietly=TRUE)
# library(Biostrings, quietly=TRUE)
opts_chunk$set(tidy=FALSE, tidy.opts=list(width.cutoff=70))

## ----eval=FALSE---------------------------------------------------------------
# # Convert a modkit TSV to a compressed .ch3 archive
# make_ch3_archive(
#   file_name   = "/path/to/modkit_calls.tsv",  # tab-delimited input
#   sample_name = "Astrocytes",                 # used in output filenames
#   out_path    = "/path/to/ch3_out",           # directory to write .ch3 files
#   short_ids   = TRUE                          # shorten read_id to reduce size
# )
# 
# # The directory now contains e.g.: Astrocytes-0.ch3, Astrocytes-1.ch3, ...
# # Build a DB from those .ch3 files:
# make_ch3_db(
#   ch3_files = "/path/to/ch3_out",   # a directory OR a vector of file paths
#   db_name   = "astro_db"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Example: directory of CH3 files; auto-derived sample names
# ch3_files <- system.file("TEST_DATA", package = "MethylSeqR")   # adjust as needed
# ch3_db    <- tempfile("example_db")
# make_ch3_db(ch3_files, ch3_db)
# 
# # Example: explicit sample names override the fallback
# make_ch3_db(
#   ch3_files = c(
#     sample1  = "/data/sample1-0.ch3",
#     sample2     = "/data/sample2.ch3"
#   ),
#   db_name   = "named_samples"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Minimal usage: default m/h classes
# summarize_ch3_positions(ch3_db)
# 
# # Custom classes, including a novel code 'a', with stricter filtering
# summarize_ch3_positions(
#   ch3_db,
#   mod_code      = c("a", "m + h"),
#   unmod_code    = "-",
#   unmod_label   = "c",
#   min_num_calls = 5
# )

## ----eval=FALSE---------------------------------------------------------------
# region_bed = system.file("Islands_hg38_test.csv", package = "MethylSeqR")
# # Minimal call
# summarize_ch3_regions(ch3_db, region_file = region_bed)
# 
# # Custom call
# summarize_ch3_regions(
#   ch3_db,
#   region_file   = region_bed,
#   join          = "inner",
#   mod_code      = c("m","h","m + h"),   # supports novel codes too, e.g., "a"
#   min_num_calls = 5
# )

## ----eval=FALSE---------------------------------------------------------------
# # 100 bp windows, non-overlapping (step_size = window_size)
# summarize_ch3_windows(ch3_db, window_size = 100, step_size = 100)
# 
# # 2 kb windows with offsets; include a novel 'a' class
# summarize_ch3_windows(
#   ch3_db,
#   window_size   = 2000,
#   step_size     = 20,
#   mod_code      = c("a", "m + h"),
#   min_num_calls = 25,
#   overwrite     = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Calculate windows
# summarize_ch3_reads(ch3_db)

## ----eval=FALSE---------------------------------------------------------------
# # Path or ch3_db object both work
# get_ch3_dbinfo("my_data.ch3.db")

## ----eval=FALSE---------------------------------------------------------------
# get_ch3_tableinfo("my_data.ch3.db", table_name = "positions")

## ----eval=FALSE---------------------------------------------------------------
# get_ch3_cols("my_data.ch3.db", "windows")

## ----eval=FALSE---------------------------------------------------------------
# n_cpg <- get_ch3_cpg_count("my_data.ch3.db", table_name = "calls")

## ----eval=FALSE---------------------------------------------------------------
# # Rename in the 'positions' table
# rename_ch3_samples(
#   "my_data.ch3.db",
#   table       = "positions",
#   samples_map = c("Cortical_Neurons" = "Cortex", "Blood" = "PBMC"),
#   preview     = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# remove_ch3_table("my_data.ch3.db", "tmp_debug_table")

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Get methylation statistics for the 'positions' call type without plotting
# calc_ch3_diff(ch3_db = ch3_db,
#               call_type = "positions",
#               cases = "blood",
#               controls = "sperm")

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Collapse significant windows
# collapse_ch3_windows("my_data.ch3.db")

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Collapse significant windows
# run_ch3_analysis(ch3_db,
#              out_path = "/Users/analysis",
#              call_type = "windows",
#              cases = c("sperm"),
#              controls = c("blood"))

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Get coverage statistics for the 'positions' call type
# plot_ch3_cov(ch3_db = ch3_db, call_type = "positions")

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Get methylation statistics for the 'positions' call type without plotting
# plot_ch3_modfrac(ch3_db = ch3_db, call_type = "positions")

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Run the correlation matrix function using the 'positions' call type and plot the results
# calc_ch3_samplecor(ch3_db = ch3_db, call_type = "positions", plot = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Calculate PCA
# plot_ch3_pca(ch3_db)

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Calculate windows
# run_ch3_qc(ch3_db)

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Get any table in your database and create a variable in R to work with
# positions = get_ch3_table("my_data.ch3.db", "positions")
# regions = get_ch3_table("my_data.ch3.db", "regions")
# 
# mod_diff_regions = get_ch3_table("my_data.ch3.db", "mod_diff_regions")
# collapsed_windows = get_ch3_table("my_data.ch3.db", "collapsed_windows")

## ----eval=FALSE---------------------------------------------------------------
# # Specify the path to the database
# ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
# 
# # Export the table to your computer
# export_ch3_table("my_data.ch3.db",
#              table = "windows",
#              out_path = "/Desktop/My_Folder/windows.csv")
# 
# export_ch3_table("my_data.ch3.db",
#              table = "mod_diff_windows",
#              out_path = "/Desktop/My_Folder/mod_diff_windows.csv")

