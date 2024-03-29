% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{SelectInformativeGenes}
\alias{SelectInformativeGenes}
\title{Select informative genes from aggregated pseudo-bulk profiles.}
\usage{
SelectInformativeGenes(
  reference,
  regress_out = NULL,
  K = 1000,
  mean_cutoff = 1,
  pseudocount = 0.01
)
}
\arguments{
\item{reference}{Tall and skinny matrix of cell type profiles.}

\item{regress_out}{factors to remove when computing gene variances}

\item{K}{how many genes to return}

\item{mean_cutoff}{Genes with mean expression below this value are excluded.}

\item{pseudocount}{cv is computed as sd / (mean + pseudocount).}
}
\description{
Select informative genes from aggregated pseudo-bulk profiles.
}
