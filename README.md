<div style="display: flex; align-items: center; justify-content: center;">
  <img src="inst/WBL_METHYLSEQR.png" alt="MethylSeqR Logo" style="width: 275px;">
</div>


# MethylSeqR

## Version 0.6.0
*Note: This is an early release - changes may occur that significantly change the functionality and structure of the data and functions. The user should be aware that subsequent releases may break code written using earlier releases.*

MethylSeqR is an R package managing Direct Whole Methylome Sequencing (dWMS) data. It creates a database, and processes it with unique options. Data can be summarized by positions, windows, or provided an annotation bed file, by unique genomic regions. The package also offers quality control functions, differential methylation, and a sliding window analysis.

## Installation

```{r, eval = FALSE}
# Install the devtools package if necessary
install.packages("devtools")

# Install MethylSeqR from GitHub
devtools::install_github("Wasatch-BioLabs/MethylSeqR", build_vignettes = TRUE)

# Access Package
library(MethylSeqR)
```

***For Linux Users:*** *System packages may need to be intalled in order to use devtools. Instructions can be found online. An example guide for this can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-16-04).*

## Paradigm 
Begin with CH3 files using the `make_pos_db()` function to generate a database. This will at first hold a positions table. 

After a positional table is made, a user can summarize by region (`summarize_regions()`) or windows (`summarize_windows()`). A differential methylation analysis can be conducted on positional, regional, or window data using `calc_mod_diff()`.

`qc_wrapper()` can be called to visually assess data. Call `export_tables()` to export any data table from the database to a file, or use `get_table()` to import as a tibble into your local environment.

#### Developed by Wasatch Biolabs.
#### Visit us on [our website](https://www.wasatchbiolabs.com/) for more details.

<div style="margin-top: 40px; text-align: center;"> <img src="inst/wbl_main_logo.png" alt="Wasatch Biolabs Logo" style="width: 200px;"> </div> 
