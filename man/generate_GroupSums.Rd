% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LSI_functions.R
\name{generate_GroupSums}
\alias{generate_GroupSums}
\title{Function to form row and column sums and means on count matrix from Seurat object.}
\usage{
generate_GroupSums(mat, groups = NULL, na.rm = TRUE, sparse = FALSE)
}
\arguments{
\item{mat}{raw count matrix.}

\item{groups}{which groups to summarize.}

\item{na.rm}{whether to remove any NA from matrix.}

\item{sparse}{whether to use sparse matrix conversion.}
}
\description{
Summarized group of count matrix for LSI
}
