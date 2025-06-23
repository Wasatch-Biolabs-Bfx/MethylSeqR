<div style="display: flex; align-items: center; justify-content: center;">
  <img src="inst/WBL_METHYLSEQR.png" alt="MethylSeqR Logo" style="width: 275px;">
</div>


# MethylSeqR

## Version 0.8.1 
**(Updated June 23 2025)**

***Note***: *This is an early release - changes may occur that significantly change the functionality and structure of the data and functions. The user should be aware that subsequent releases may break code written using earlier releases.*

MethylSeqR is an R package managing Direct Whole Methylome Sequencing (dWMS) data. It creates a database, and processes it with unique options. Data can be summarized by positions, windows, or provided an annotation bed file, by unique genomic regions. The package also offers quality control functions, differential methylation, and a sliding window analysis.

## Installation

```{r, eval = FALSE}
# Install the devtools package if necessary
install.packages("devtools")

# Install MethylSeqR from GitHub
devtools::install_github("Wasatch-Biolabs-Bfx/MethylSeqR", build_vignettes = FALSE)

# Access Package
library(MethylSeqR)
```

***For Linux Users:*** *System packages may need to be intalled in order to use devtools. Instructions can be found online. An example guide for this can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-16-04).*

## Instructions

Begin with CH3 files created by Wasatch Biolabs and build a database using the `make_ch3_db()`. This will at first hold a calls table. If you would like to see key stats on your CH3 file, call `get_ch3_stats()`. This shows information like CpG coverage, calls by flag value, high confidence calls, high quality calls, and average read length.

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
  ch3_files = "/ch3_files_directory", 
             db_name = "my_data") |>
  summarize_ch3_windows() |>
  calc_ch3_diff(call_type = "windows",
              cases =
                c("sperm"),
              controls =
                c("blood")) |>
  collapse_ch3_windows() 
  
# Build and analyze through separate lines
ch3_db <- make_ch3_db(ch3_files = "/ch3_files_directory", db_name = "my_data")
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
export_ch3_table(ch3_db, "collapsed_windows", "/results_directory")

#DONE! Data has been analyzed and exported!

```

***Warning*** *- If using samples other than human or with unique chromosome names, remember to adjust the chrs argument in each function!*

## Vignette
To get detailed instructions and help working through the package, build and follow along the vignette by calling:
```{r, eval = FALSE}
browseVignettes("MethylSeqR")
```
in R and click on HTML in your browser. Or, to browse the vignette in your R environment, call

```{r, eval = FALSE}
vignette("MethylSeqRWalkthrough")
```

#### Developed by Wasatch Biolabs.
#### Visit us on [our website](https://www.wasatchbiolabs.com/) for more details.

<div style="margin-top: 40px; text-align: center;"> <img src="inst/wbl_main_logo.png" alt="Wasatch Biolabs Logo" style="width: 200px;"> </div> 
