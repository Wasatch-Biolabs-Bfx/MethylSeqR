<div style="display: flex; align-items: center; justify-content: center;">
  <img src="inst/WBL_METHYLSEQR.png" alt="MethylSeqR Logo" style="width: 275px;">
</div>


# MethylSeqR

## Version 0.8.1 
**(Updated December 3rd 2025)**

***Note***: *This is an early release - changes may occur that significantly change the functionality and structure of the data and functions. The user should be aware that subsequent releases may break code written using earlier releases.*

MethylSeqR is an R package managing Direct Whole Methylome Sequencing (dWMS) data. It creates a database, and processes it with unique options. Data can be summarized by positions, windows, or provided an annotation bed file, by unique genomic regions. The package also offers quality control functions, differential methylation, and a sliding window analysis.



## Installation

For easy visualization and data management, we encourage you to download and use RStudio.
RStudio download instructions can be found [here](https://posit.co/download/rstudio-desktop/).

```{r, eval = FALSE}
# Install the devtools package if necessary
install.packages("devtools")

# Install MethylSeqR from GitHub
devtools::install_github("Wasatch-Biolabs-Bfx/MethylSeqR", build_vignettes = FALSE)

# Access Package
library(MethylSeqR)
```

***For Linux Users:*** *System packages may need to be intalled in order to use devtools. Instructions can be found online. An example guide for this can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-16-04).*



## Creating .ch3 Files

If your sequencing data is coming from [Wasatch Biolabs](https://www.wasatchbiolabs.com/), 
`.ch3` files can be delivered directly with your batch upon request. 
These files contain methylation calls from a third-generation sequencing run. 
Typically, one .ch3 file is created per sample, and each file compresses large amounts of raw data into a usable intermediate format.

If .ch3 files were not provided, you can build them yourself using the make_ch3_archive() function in this package. 
This takes a tab-delimited file of methylation calls (such as the output of modkit extract-calls) and compresses it into .ch3 format for downstream use.

```{r, eval=FALSE}
# Convert a calls.tsv file to compressed .ch3 format
make_ch3_archive(
  file_name   = "calls.tsv",   # input modkit calls file
  sample_name = "sample1",     # will be embedded in output filenames
  out_path    = "output_dir/", # where to write .ch3 files
  short_ids   = TRUE           # optionally shorten read_id to reduce size
)
```

### Important:

* .ch3 files use 0-based, half-open genome coordinates (like BED files).
* Example: a CpG at base 1000 (1-based) will appear as start=999, end=1001.
* Each output archive is written in compressed Parquet format (zstd), typically producing multiple .ch3 files per sample (e.g., sample1-0.ch3, sample1-1.ch3).

## Example Data

If you’d like to test the package without generating your own .ch3 files, small example data are included with the package in:

`MethylSeqR/inst/extdata/ch3_files/`

You can download these example .ch3 files directly from the GitHub repository under the inst/extdata/ch3_files folder.

These test files can be used to practice building a database and running the full analysis workflow.

```{r, eval = FALSE}
# Example: build a test database using included example files
example_path <- system.file("extdata/ch3_files", package = "MethylSeqR")
ch3_db <- make_ch3_db(example_path, db_name = "example_db")
```


## Instructions
Begin with .ch3 files created by Wasatch Biolabs or individually, and build a database using `make_ch3_db()`. 
This will at first hold a calls table. If you would like to see key stats on your CH3 file, call `get_ch3_stats()`. This shows information like CpG coverage, calls by flag value, high confidence calls, high quality calls, and average read length.

### Summarize Data
After a database is created, a user can summarize their data by position (`summarize_ch3_positions()`), by regions (`summarize_ch3_regions()`), by windows (`summarize_ch3_windows()`), or by reads ((`summarize_ch3_reads()`). 



### Differential Methylation
A differential methylation analysis can be conducted on positional, regional, or window data using `calc_ch3_diff()`. After calculating methylation differences between windows, use `collapse_ch3_windows()` to collapse significant windows in a methylation dataset. This merges contiguous regions that meet the specified criteria.



### Get Database Stats
If you would like to see key stats on your database at any time, including what unique sample names are in the data for a differential analysis, call `get_ch3_dbinfo()`. 
To see what columns are in a table in your database and how many records (rows) there are, call `get_ch3_tableinfo()` with your database and desired table name.



### Quality Control
`run_ch3_qc()` can be called to visually assess any data. Running a QC can take a long time on large data, so set the argument `max_rows` to a reasonable value (ex. 1000) to assess data faster. To view and extract a table, call `export_table()` to export any data table from the database to a file, or use `get_table()` to import as a tibble into your local environment. Similarily, use `max_calls` if you are fine with a smaller, randomized set of data.

If you would like to run everything in one command, call `run_ch3_analysis()`.



## Paradigm
You can pipe your functions together, or feel free to call each function one line at a time. Below are two examples of this.
```{r, eval = FALSE}
setwd("/home/directory/analysis")

# Build database and run analysis in a pipe
ch3_db <- make_ch3_db(
  ch3_files = "../ch3_files_directory", 
             db_name = "my_data") |>
  summarize_ch3_windows() |>
  calc_ch3_diff(call_type = "windows",
              cases =
                c("sperm"),
              controls =
                c("blood")) |>
  collapse_ch3_windows() 
  
# Build and analyze through separate lines
ch3_db <- make_ch3_db(ch3_files = "../ch3_files_directory", db_name = "my_data")
ch3_db <- summarize_ch3_windows(ch3_db)
ch3_db <- calc_ch3_diff(ch3_db, call_type = "windows", cases = c("sperm"), controls = c("blood"))
ch3_db <- collapse_ch3_windows(ch3_db) 

# Run entire differential methylation analysis in one function
run_ch3_analysis(ch3_db, 
             out_path = "/analysis",
             call_type = "windows",
             cases = c("sperm"),
             controls = c("blood"))

# Check to see what's in your database at any time
get_ch3_dbinfo(ch3_db)

# Check to see the columns in any table 
get_ch3_tableinfo(ch3_db, "windows")

# Export differentially methylated data to your computer
export_ch3_table(ch3_db, "collapsed_windows", "../results_directory")

# OR, to work with your data locally in your R environment
positions <- get_ch3_table(ch3_db, "positions")
windows <- get_ch3_table(ch3_db, "windows")
regions <- get_ch3_table(ch3_db, "regions")

#DONE! Data has been analyzed and exported!

```

***Warning*** *- If using samples other than human or with unique chromosome names, remember to adjust the chrs argument in each function!*




### Convenience Functions

MethylSeqR also provides a few helper utilities to make it easier to inspect and manage your database:

```{r, eval=FALSE}
# View all column names in a given table
get_ch3_cols(ch3_db, "calls")

# Count unique CpG sites (based on start/end)
get_ch3_cpg_count(ch3_db, table_name = "calls")

# Show high-level database statistics (size, tables, sample names)
get_ch3_dbinfo(ch3_db)

# Get detailed information about a specific table
get_ch3_tableinfo(ch3_db, "positions")

# Rename sample names inside any table
rename_ch3_samples(ch3_db, "positions",
                   samples_map = c("old_name" = "new_name"))

# Remove a table from the database
remove_ch3_table(ch3_db, "temp_table")
```

## Getting Help
To see detailed documentation on a specific function in R, call `?{function}`. Example:
```{r, eval = FALSE}
?make_ch3_db()
```
This will render development documentation for that function in the Help tab in Rstudio

### Vignette
To get detailed instructions and help working through the package, download and view our vignette in this github repo at `docs/MethylSeqRWalkthrough.html`.

Or, build and follow along the vignette by calling:
```{r, eval = FALSE}
browseVignettes("MethylSeqR")
```
in R and click on HTML in your browser. Or, to browse the vignette in your R environment, call

```{r, eval = FALSE}
vignette("MethylSeqRWalkthrough")
```

## License

This package is distributed under the **Personal and Internal Research License (v1.1)**, based on the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0), with an additional grant of rights for internal research use provided by Wasatch Biolabs.

- ✅ Free for personal use, academic research, educational, and other noncommercial purposes.  
- ✅ Commercial organizations may also use it internally for research and evaluation.  
- ❌ Not permitted for commercial services, fee-for-service analysis, redistribution, or production use.  

See the full license text in the [LICENSE](./LICENSE) file included with this repository.

If you have any suggestions or requested features, please email Jonathon Hill at jhill@byu.edu.

#### Developed by Wasatch Biolabs.
#### Visit us on [our website](https://www.wasatchbiolabs.com/) for more details.

<div style="margin-top: 40px; text-align: center;"> <img src="inst/wbl_main_logo.png" alt="Wasatch Biolabs Logo" style="width: 200px;"> </div> 
