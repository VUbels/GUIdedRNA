# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(Matrix)
library(dplyr)
library(cowplot)
library(Seurat)
library(harmony)
library(edgeR)
library(celda)

# Load problematic packages with error handling
safe_library <- function(package) {
  tryCatch({
    library(package, character.only = TRUE)
    message(paste("✓ Loaded:", package))
    return(TRUE)
  }, error = function(e) {
    warning(paste("✗ Failed to load:", package, "- Error:", e$message))
    return(FALSE)
  })
}

# Try to load the problematic packages
GENOMIC_FEATURES_AVAILABLE <- safe_library("GenomicFeatures")
TXDB_AVAILABLE <- safe_library("TxDb.Hsapiens.UCSC.hg38.knownGene")
DOUBLETFINDER_AVAILABLE <- safe_library("DoubletFinder")

# Load GUIdedRNA package
library(GUIdedRNA)

# Function to check if a required package is available
check_package_availability <- function(package_name) {
  switch(package_name,
         "GenomicFeatures" = GENOMIC_FEATURES_AVAILABLE,
         "TxDb.Hsapiens.UCSC.hg38.knownGene" = TXDB_AVAILABLE,
         "DoubletFinder" = DOUBLETFINDER_AVAILABLE,
         TRUE  # Default to TRUE for other packages
  )
}

# Setup volumes for file browser
volumes <- c()

# Add home directory
volumes <- c(volumes, Home = fs::path_home())

# Add root for Linux/Mac
if (.Platform$OS.type == "unix") {
  volumes <- c(volumes, "Root" = "/")
  
  # Add common WSL paths if running in WSL
  wsl_check <- system("grep -q Microsoft /proc/version", ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (wsl_check == 0) {  # Running in WSL
    # Add Windows drives that might be mounted
    if (dir.exists("/mnt/c")) volumes <- c(volumes, "C:" = "/mnt/c")
    if (dir.exists("/mnt/d")) volumes <- c(volumes, "D:" = "/mnt/d")
    if (dir.exists("/mnt/e")) volumes <- c(volumes, "E:" = "/mnt/e")
    if (dir.exists("/mnt/f")) volumes <- c(volumes, "F:" = "/mnt/f")
    if (dir.exists("/mnt/g")) volumes <- c(volumes, "G:" = "/mnt/g")
    
    if (dir.exists("/mnt")) volumes <- c(volumes, "Windows Drives" = "/mnt")
  }
  
  # Add other common Unix locations
  if (dir.exists("/media")) volumes <- c(volumes, "Media" = "/media")
  if (dir.exists("/home")) volumes <- c(volumes, "Home" = "/home")
}

# Add Windows drives if on Windows
if (.Platform$OS.type == "windows") {
  for (letter in LETTERS) {
    drive <- paste0(letter, ":")
    if (dir.exists(drive)) {
      volumes <- c(volumes, drive)
      names(volumes)[length(volumes)] <- paste0(letter, ":")
    }
  }
}

# Docker volume support
if (dir.exists("/data")) {
  volumes <- c(volumes, "Data" = "/data")
}
if (dir.exists("/output")) {
  volumes <- c(volumes, "Output" = "/output")
}

message("Global.R loaded successfully")
message(paste("Available volumes:", paste(names(volumes), collapse = ", ")))
message(paste("GenomicFeatures available:", GENOMIC_FEATURES_AVAILABLE))
message(paste("TxDb.Hsapiens.UCSC.hg38.knownGene available:", TXDB_AVAILABLE))
message(paste("DoubletFinder available:", DOUBLETFINDER_AVAILABLE))