<div style="display: flex; align-items: center; justify-content: center;">
  <img src="inst/WBL_METHYLSEQR.png" alt="MethylSeqR Logo" style="width: 275px;">
</div>


# MethylSeqR

## Version 0.8.0 
**(Updated April 3 2025)**

***Note***: *This is an early release - changes may occur that significantly change the functionality and structure of the data and functions. The user should be aware that subsequent releases may break code written using earlier releases.*

MethylSeqR is an R package managing Direct Whole Methylome Sequencing (dWMS) data. It creates a database, and processes it with unique options. Data can be summarized by positions, windows, or provided an annotation bed file, by unique genomic regions. The package also offers quality control functions, differential methylation, and a sliding window analysis.

## Installation

```{r, eval = FALSE}
# Install the devtools package if necessary
install.packages("devtools")

# Install MethylSeqR from GitHub
devtools::install_github("Wasatch-Biolabs-Bfx/MethylSeqR", build_vignettes = TRUE)

# Access Package
library(MethylSeqR)
```

***For Linux Users:*** *System packages may need to be intalled in order to use devtools. Instructions can be found online. An example guide for this can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-16-04).*

## Paradigm 
Begin with CH3 files created by Wasatch Biolabs and build a database using the `make_ch3_db()`. This will at first hold a calls table. If you would like to see key stats on your CH3 file, call `get_ch3_stats()`. This shows information like CpG coverage, calls by flag value, high confidence calls, high quality calls, and average read length.

After a database is created, a user can summarize their data by position (`summarize_positions()`), by regions (`summarize_regions()`) or windows (`summarize_windows()`). A differential methylation analysis can be conducted on positional, regional, or window data using `calc_mod_diff()`. If you would like to see key stats on your database, including what unique sample names are in the data for a differential analysis, call `get_db_stats()`.

`run_qc()` can be called to visually assess any data. To view and extract a table, call `export_table()` to export any data table from the database to a file, or use `get_table()` to import as a tibble into your local environment.

***Warning*** *- If using samples other than human or with unique chromosome names, remember to adjust the chrs argument in each function!*

## Vignette
To get detailed instructions and help working through the package, build and follow along the vignette by calling:
```{r, eval = FALSE}
browseVignettes("MethylSeqR")
```
and click on HTML in your browser. Or, to browse the vignette in your R environment, call

```{r, eval = FALSE}
vignette("MethylSeqRWalkthrough")
```

#### Developed by Wasatch Biolabs.
#### Visit us on [our website](https://www.wasatchbiolabs.com/) for more details.

<div style="margin-top: 40px; text-align: center;"> <img src="inst/wbl_main_logo.png" alt="Wasatch Biolabs Logo" style="width: 200px;"> </div> 
