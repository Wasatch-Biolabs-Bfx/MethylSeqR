.onAttach <- function(lib, pkg)
{
  msg1 <- paste0(
    "========================================\n",    
    "  ╔╦╗╔═╗╔╦╗╦ ╦╦ ╦╦  ╔═╗╔═╗╔═╗ ╦═╗\n",
    "  ║║║║╣  ║ ╠═╣╚╦╝║  ╚═╗║╣ ║═╬╗╠╦╝\n",
    "  ╩ ╩╚═╝ ╩ ╩ ╩ ╩ ╩═╝╚═╝╚═╝╚═╝╚╩╚═ ", packageVersion(pkg), "\n",
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
