% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_positions.R
\name{summarize_positions}
\alias{summarize_positions}
\title{Summarize Methylation Positions from Calls Table}
\usage{
summarize_positions(ch3_db, mod_type = c("c", "m", "h", "mh"), min_cov = 1)
}
\arguments{
\item{ch3_db}{A DuckDB database connection or file path (character) to the `.ch3.db` file.}

\item{mod_type}{A character vector specifying which DNA modification types to include in the summarized data. 
The vector can include values like `"m"`, `"h"`, or `"mh"`.}

\item{min_cov}{Minimum coverage required to include a position in the summary. Default is 1.}
}
\value{
The modified `ch3_db` object with the updated `positions` table.
}
\description{
This function processes methylation call data from a DuckDB database, summarizing 
coverage and call counts into a `positions` table.
}
\details{
The function:
- Connects to the DuckDB database.
- Summarizes methylation calls by `sample_name`, `chrom`, `start`, and `end`.
- Computes coverage and call counts (`m`, `h`, `mh`, `c`).
- Filters based on `min_cov`.
- Creates a `positions` table in the database.

A progress bar is displayed during execution.
}
\examples{
summarize_positions("example_ch3.ch3.db", min_cov = 1)

}
