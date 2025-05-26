library(GUIdedRNA)
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(Matrix)
library(plyr)
library(dplyr)
library(irlba)

# Load packages with error handling
safe_library <- function(package) {
  tryCatch({
    library(package, character.only = TRUE)
    message(paste("✓ Loaded:", package))
  }, error = function(e) {
    warning(paste("✗ Failed to load:", package, "- Error:", e$message))
  })
}

# Essential packages that must load
essential_packages <- c("edgeR", "Seurat", "SeuratObject", "cowplot", "fs", "R6", "harmony")
for(pkg in essential_packages) {
  safe_library(pkg)
}

# Optional packages that can fail gracefully
optional_packages <- c("celda", "decontX", "GenomicFeatures", "GenomicRanges", 
                      "AnnotationDbi", "sparseMatrixStats", "matrixStats",
                      "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
for(pkg in optional_packages) {
  safe_library(pkg)
}

# Docker-compatible volume setup
setup_docker_volumes <- function() {
  volumes <- c()
  
  # Check if we're running in Docker (look for mounted volumes)
  if (dir.exists("/data")) {
    volumes <- c(volumes, "Data Volume" = "/data")
  }
  
  if (dir.exists("/output")) {
    volumes <- c(volumes, "Output Volume" = "/output")
  }
  
  # Always include root and common directories
  volumes <- c(volumes, "Root" = "/")
  
  # Check for host system mounts (from docker volume mapping)
  if (dir.exists("/host_root")) {
    volumes <- c(volumes, "Host Root" = "/host_root")
  }
  
  if (dir.exists("/host_home")) {
    volumes <- c(volumes, "Host Home" = "/host_home")
  }
  
  if (dir.exists("/host_mnt")) {
    volumes <- c(volumes, "Host Mounted Drives" = "/host_mnt")
    
    # Add individual drives if they exist in host_mnt
    for (letter in c("c", "d", "e", "f", "g", "h")) {
      drive_path <- paste0("/host_mnt/", letter)
      if (dir.exists(drive_path)) {
        drive_name <- paste0("Windows ", toupper(letter), ":")
        volumes <- c(volumes, drive_name)
        names(volumes)[length(volumes)] <- drive_name
        volumes[length(volumes)] <- drive_path
      }
    }
  }
  
  # Check for direct mounted Windows drives
  if (dir.exists("/mnt")) {
    volumes <- c(volumes, "Direct Mounted Drives" = "/mnt")
    
    # Add individual drives if they exist
    for (letter in c("c", "d", "e", "f", "g", "h")) {
      drive_path <- paste0("/mnt/", letter)
      if (dir.exists(drive_path)) {
        drive_name <- paste0("Drive ", toupper(letter), ":")
        volumes <- c(volumes, drive_name)
        names(volumes)[length(volumes)] <- drive_name
        volumes[length(volumes)] <- drive_path
      }
    }
  }
  
  # Add home directory if it exists
  if (dir.exists("/home")) {
    volumes <- c(volumes, "Home" = "/home")
  }
  
  # Add tmp directory
  if (dir.exists("/tmp")) {
    volumes <- c(volumes, "Temp" = "/tmp")
  }
  
  # Check for media directories
  if (dir.exists("/media")) {
    volumes <- c(volumes, "Media" = "/media")
  }
  
  return(volumes)
}

# Set up volumes globally
.docker_volumes <- setup_docker_volumes()

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
message(paste("Available volumes:", paste(names(.docker_volumes), collapse = ", ")))