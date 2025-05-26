# inst/shiny-app/global.R
library(GUIdedRNA)
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(Matrix)
library(dplyr)
library(irlba)
library(edgeR)
library(DoubletFinder)
library(celda)
library(decontX)
library(Seurat)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(cowplot)
library(plyr)
library(SeuratObject)
library(fs)
library(R6)
library(harmony)
library(GenomicFeatures)
library(GenomicRanges)
library(AnnotationDbi)
library(sparseMatrixStats)
library(matrixStats)

# Global message queue - using a hidden environment variable for storage
.message_env <- new.env()
.message_env$queue <- character(0)

# This function will be replaced when the app initializes
send_message <- function(msg, tab_id = NULL) {
  message(paste("Message queued:", msg))
}

get_messages <- function() {
  msgs <- .message_env$queue
  .message_env$queue <- character(0)
  return(msgs)
}

message("GUIdedRNA global.R loaded successfully")