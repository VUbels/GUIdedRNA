% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plotUMAP}
\alias{plotUMAP}
\title{Function to run UMAP on LSI results.}
\usage{
plotUMAP(
  df,
  dataType = "qualitative",
  cmap = NULL,
  covarLabel = "",
  point_size = 0.5,
  namedColors = FALSE,
  plotTitle = NULL,
  colorLims = NULL,
  na.value = "grey35",
  useRaster = TRUE
)
}
\arguments{
\item{df}{Seurat object.}

\item{cmap}{colour map to pass to function.}

\item{namedColors}{vector colour map.}

\item{plotTitle}{title to plot.}

\item{colorLims}{colour limits for minimal and maximal percentage based on expression data.}

\item{na.value}{colour for cells not expression gene.}

\item{useRaster}{raster plot.}

\item{rawCounts}{count matrix from Seurat object.}

\item{coverLabel}{cover label.}

\item{point.size}{individual point size per cell.}
}
\description{
Function to run UMAP on LSI results.
}
