# global.R
library(R6)
library(Seurat)
library(ggplot2)
library(DT)
library(Matrix)
library(dplyr)
library(DoubletFinder)
library(decontX)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(irlba)
library(edgeR)

# Global message queue
message("Loading global.R file...")
source("../R/preprocessing_functions.R")
source("../R/LSI_functions.R")

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
