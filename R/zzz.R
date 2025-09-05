#' MethylSeqR: Methylation Data Analysis from Nanopore Sequencing
#'
#' MethylSeqR is an R package developed by Wasatch Biolabs for efficient
#' preprocessing, summarization, and visualization of methylation data.
#' It supports native DNA methylation sequencing workflows and is
#' optimized for high-throughput pipelines.
#'
#' @section Features:
#' - Database-backed storage using DuckDB
#' - Sliding window, positional, regional, and read-level summarization of methylation levels
#' - Support for multiple methylation contexts
#' - Quality Control Visualization
#' - Efficient Differential Methylation Analysis 
#' 
#' @author
#' Jonathon Hill \email{jonathon@wasatchbiolabs.com} (aut, cre)
#' Hailey Zimmerman \email{hailey@renewbt.com} (aut)
#' 
#' @name MethylSeqR
#'
#' @references
#' Wasatch Biolabs (2025). *MethylSeqR: Tools for Methylation Analysis in Clinical and Research Settings.*
#' https://www.wasatchbiolabs.com/
#' 
#' For bug reports and feature requests:  
#' https://github.com/Wasatch-Biolabs-Bfx/MethylSeqR
"_PACKAGE"

.onAttach <- function(lib, pkg)
{
  pv <- packageVersion(pkg)
  msg1 <- 
    glue::glue("========================================",    
               "  ╔╦╗╔═╗╔╦╗╦ ╦╦ ╦╦  ╔═╗╔═╗╔═╗ ╦═╗       ",
               "  ║║║║╣  ║ ╠═╣╚╦╝║  ╚═╗║╣ ║═╬╗╠╦╝       ",
               "  ╩ ╩╚═╝ ╩ ╩ ╩ ╩ ╩═╝╚═╝╚═╝╚═╝╚╩╚═{pv}   ",
               "========================================",
               .sep = "\n")
  
  msg2 <- paste("MethylSeqR v", packageVersion(pkg))
  
  if (interactive()) {
    packageStartupMessage(msg1)
  } else {
    packageStartupMessage(msg2)
  }
  
  packageStartupMessage("Created by Wasatch Biolabs")
  packageStartupMessage("Research & Clinical Nanopore Sequencing")
  packageStartupMessage("www.wasatchbiolabs.com")
  packageStartupMessage("")
  packageStartupMessage("This package is licensed for personal") 
  packageStartupMessage("and internal research use only.")
  packageStartupMessage("See the LICENSE file or visit")
  packageStartupMessage("https://your-website-or-repo/LICENSE.")
}
