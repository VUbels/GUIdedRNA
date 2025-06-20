% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LSI_functions.R
\name{generate_VariableGenes}
\alias{generate_VariableGenes}
\title{Variable gene function selection enabling.}
\usage{
generate_VariableGenes(mat, nvar = 2000, blacklist = NULL)
}
\arguments{
\item{mat}{raw count matrix from Seurat object.}

\item{nvar}{number of variable genes to be found.}

\item{blacklist}{list of genes to exclude from variable gene analysis.}
}
\value{
Variable list of genes
}
\description{
Generates variable genes from matrix
}
