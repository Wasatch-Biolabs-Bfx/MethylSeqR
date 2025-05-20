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
#' @docType Package
#' @name MethylSeqR
#'
#' @references
#' Wasatch Biolabs (2025). *MethylSeqR: Tools for Methylation Analysis in Clinical and Research Settings.*
#' https://www.wasatchbiolabs.com/
#' 
#' For bug reports and feature requests:  
#' https://github.com/Wasatch-Biolabs-Bfx/MethylSeqR
NULL

.onAttach <- function(lib, pkg)
{
  msg1 <- paste0(
    "========================================\n",    
    "  ╔╦╗╔═╗╔╦╗╦ ╦╦ ╦╦  ╔═╗╔═╗╔═╗ ╦═╗\n",
    "  ║║║║╣  ║ ╠═╣╚╦╝║  ╚═╗║╣ ║═╬╗╠╦╝\n",
    "  ╩ ╩╚═╝ ╩ ╩ ╩ ╩ ╩═╝╚═╝╚═╝╚═╝╚╩╚═", packageVersion(pkg), "\n",
    "========================================\n",
    "Created by Wasatch Biolabs\n",
    "Research & Clinical Nanopore Sequencing\n",
    "www.wasatchbiolabs.com\n")
  
  msg2 <- paste("MethylSeqR", packageVersion(pkg), "www.wasatchbiolabs.com")
  
  if (interactive()) {
    packageStartupMessage(msg1)
  } else {
    packageStartupMessage(msg2)
  }
}
