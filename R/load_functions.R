#' Variable gene function selection enabling.
#'
#' @description Reads folder, or folder of folder or 10X files to create Seurat objects
#' @returns List of seurat objects
#' @export

read_folder <- function(folder_path) {
  # Initialize a list to store Seurat objects
  seurat_list <- list()
  
  # 1. Process files in the main directory with pattern: name_matrix.mtx.gz, etc.
  main_files <- list.files(folder_path, pattern = "(_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz)$", 
                           full.names = TRUE, recursive = FALSE)
  
  # Find base names for files in the main directory
  main_base_names <- unique(sapply(strsplit(basename(main_files), 
                                            "_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz"), 
                                   `[`, 1))
  
  # Process datasets in the main directory
  for (base in main_base_names) {
    
    mtx_file <- file.path(folder_path, paste0(base, "_matrix.mtx.gz"))
    features_file <- file.path(folder_path, paste0(base, "_features.tsv.gz"))
    cells_file <- file.path(folder_path, paste0(base, "_barcodes.tsv.gz"))
    
    # Check if all files exist
    if (file.exists(mtx_file) && file.exists(features_file) && file.exists(cells_file)) {
      # Read the matrix
      tryCatch({
        # Read the matrix
        counts <- Seurat::ReadMtx(
          mtx = mtx_file,
          features = features_file,
          cells = cells_file
        )
        
        # Create Seurat object
        seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = base)
        
        # Calculate percent mitochondrial
        seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
        
        # Add to list
        seurat_list[[base]] <- seurat_obj
        
        message(paste("Successfully processed dataset:", base))
      }, error = function(e) {
        message(paste("Error processing", base, ":", e$message))
      })
    }
  }

  # Return the complete list of Seurat objects
  if (length(seurat_list) == 0) {
    message("No valid datasets found in the specified directory.")
    return(NULL)
  } else if (length(seurat_list) == 1) {
    message("Returning a single Seurat object.")
    return(seurat_list[[1]])  # Return the single object directly
  } else {
    message(paste("Returning a list of", length(seurat_list), "Seurat objects."))
    return(seurat_list)
  }
}