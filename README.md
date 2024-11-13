<img src="inst/wbl_main_logo.png" alt="Main Logo" width="500" height="150">

# MethylSeqR

## Version 0.0.5
*Note: This is an early release - changes may occur that significantly change the functionality and structure of the data and functions. The user should be aware that subsequent releases may break code written using earlier releases.*

MethylSeqR is an R package managing Direct Whole Methylome Sequencing (dWMS) data. It creates a database, and processes it with unique options. Data can be summarized by positions, windows, or provided an annotation bed file, by unique genomic regions. The package also offers quality control functions, differential methylation, and a sliding window analysis.

## Installation

```{r, eval = FALSE}
# Install the devtools package if necessary
# install.packages("devtools")
# library(devtools)

# Install MethylSeqR from GitHub
devtools::install_github("https://github.com/Wasatch-BioLabs/MethylSeqR.git")

# Access Package
library(MethylSeqR)

```

## Paradigm 
Begin with CH3 files using the `make_pos_db()` function to generate a database. This will at first hold a positions table. After that, you can use `qc_wrapper()` to visually assess positional data, or summarize by region (`summarize_regions()`) or windows (`summarize_windows`).

