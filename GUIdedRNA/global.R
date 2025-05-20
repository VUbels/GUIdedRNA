# global.R
library(R6)
library(DT)


library(ggplot2)
library(Matrix)
library(dplyr)
library(irlba)
library(edgeR)

library(DoubletFinder)
library(decontX)
library(Seurat)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Global message queue
message("Loading global.R file...")
source("../R/preprocessing_functions.R")
source("../R/LSI_functions.R")
source("../R/plotting_functions.R")

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
