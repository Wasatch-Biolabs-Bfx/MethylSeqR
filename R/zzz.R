.onAttach <- function(lib, pkg)
{
  pv <- packageVersion(pkg)
  msg1 <- 
    glue::glue("========================================",    
               "  ╔╦╗╔═╗╔╦╗╦ ╦╦ ╦╦  ╔═╗╔═╗╔═╗ ╦═╗",
               "  ║║║║╣  ║ ╠═╣╚╦╝║  ╚═╗║╣ ║═╬╗╠╦╝",
               "  ╩ ╩╚═╝ ╩ ╩ ╩ ╩ ╩═╝╚═╝╚═╝╚═╝╚╩╚═{pv}",
               "========================================")
  
  msg2 <- paste("MethylSeqR", packageVersion(pkg))
  
  if (interactive()) {
    packageStartupMessage(msg1)
  } else {
    packageStartupMessage(msg2)
  }
  
  packageStartupMessage("Created by Wasatch Biolabs")
  packageStartupMessage("Research & Clinical Nanopore Sequencing")
  packageStartupMessage("www.wasatchbiolabs.com\n")
  packageStartupMessage("This package is licensed for personal") 
  packageStartupMessage("and internal research use only.")
  packageStartupMessage("See the LICENSE file or visit")
  packageStartupMessage("https://your-website-or-repo/LICENSE.")
}
