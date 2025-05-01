.onLoad <- function(libname, pkgname) {
  # Check if 'duckdb' is installed
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    message("The 'duckdb' package is required but not installed. Installing it now...")
    install.packages("duckdb")
  }
}