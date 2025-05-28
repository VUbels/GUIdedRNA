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

# Enhanced volume detection for Docker environments
setup_volumes <- function() {
  volumes <- c()
  
  # Always include home directory if accessible
  tryCatch({
    home_path <- fs::path_home()
    if (dir.exists(home_path)) {
      volumes <- c(volumes, Home = home_path)
    }
  }, error = function(e) {
    message("Could not access home directory")
  })
  
  # Docker-specific volumes (mounted by docker-compose)
  docker_mounts <- list(
    "Data" = "/data",
    "Output" = "/output", 
    "Host System" = "/host_system",
    "Host Home" = "/host_home",
    "Host Media" = "/host_media",
    "Host Mnt" = "/host_mnt"
  )
  
  for (name in names(docker_mounts)) {
    path <- docker_mounts[[name]]
    if (dir.exists(path)) {
      volumes <- c(volumes, path)
      names(volumes)[length(volumes)] <- name
      message(paste("✓ Found Docker mount:", name, "at", path))
    }
  }
  
  # Windows drives mounted through WSL/Docker
  windows_drives <- c("C", "D", "E", "F", "G", "H", "I", "J")
  for (drive in windows_drives) {
    # Check multiple possible mount points
    mount_points <- c(
      paste0("/host_drives/", drive),
      paste0("/mnt/", tolower(drive)),
      paste0("/host_mnt/", tolower(drive))
    )
    
    for (mount_point in mount_points) {
      if (dir.exists(mount_point)) {
        volumes <- c(volumes, mount_point)
        names(volumes)[length(volumes)] <- paste0(drive, ":")
        message(paste("✓ Found Windows drive:", drive, "at", mount_point))
        break  # Use first available mount point
      }
    }
  }
  
  # Standard Unix/Linux paths
  if (.Platform$OS.type == "unix") {
    unix_paths <- list(
      "Root" = "/",
      "Media" = "/media",
      "Mnt" = "/mnt",
      "Volumes" = "/Volumes"  # macOS
    )
    
    for (name in names(unix_paths)) {
      path <- unix_paths[[name]]
      if (dir.exists(path) && !name %in% names(volumes)) {
        volumes <- c(volumes, path)
        names(volumes)[length(volumes)] <- name
      }
    }
  }
  
  # Windows native paths (if running on Windows directly)
  if (.Platform$OS.type == "windows") {
    for (letter in LETTERS) {
      drive <- paste0(letter, ":")
      if (dir.exists(drive)) {
        volumes <- c(volumes, drive)
        names(volumes)[length(volumes)] <- paste0(letter, ":")
      }
    }
  }
  
  # Remove duplicates and invalid paths
  valid_volumes <- c()
  for (i in seq_along(volumes)) {
    path <- volumes[i]
    name <- names(volumes)[i]
    
    if (dir.exists(path) && !path %in% valid_volumes) {
      valid_volumes <- c(valid_volumes, path)
      names(valid_volumes)[length(valid_volumes)] <- name
    }
  }
  
  return(valid_volumes)
}

# Setup volumes
volumes <- setup_volumes()

message("Global.R loaded successfully")
message(paste("Available volumes:", paste(names(volumes), collapse = ", ")))
message(paste("Volume paths:", paste(volumes, collapse = ", "))) 
message(paste("GenomicFeatures available:", GENOMIC_FEATURES_AVAILABLE))
message(paste("TxDb.Hsapiens.UCSC.hg38.knownGene available:", TXDB_AVAILABLE))
message(paste("DoubletFinder available:", DOUBLETFINDER_AVAILABLE))