#' Assigns expected doublet ground variable
#'
#' Function returns the expected doublet rate based on 10X data dependent on cell count rather than flat rate.
#' Increasing cell count in smaller samples.
#' @export

assign_ExpectedDoublet <- function(val) {
  # Initialize expected_doublet value
  expected_doublet <- NA
  
  # Check the range of value_A and assign the corresponding value to expected_doublet value based on 10x information
  if (val >= 0 & val <= 500) {
    expected_doublet <- 0.004
  } else if (val >= 501 & val <= 1000) {
    expected_doublet <- 0.008
  } else if (val >= 1001 & val <= 2000) {
    expected_doublet <- 0.016
  } else if (val >= 2001 & val <= 3000) {
    expected_doublet <- 0.023
  } else if (val >= 3001 & val <= 4000) {
    expected_doublet <- 0.031
  } else if (val >= 4001 & val <= 5000) {
    expected_doublet <- 0.039
  } else if (val >= 5001 & val <= 6000) {
    expected_doublet <- 0.046
  } else if (val >= 6001 & val <= 7000) {
    expected_doublet <- 0.054
  } else if (val >= 6001 & val <= 8000) {
    expected_doublet <- 0.061
  } else if (val >= 6001 & val <= 9000) {
    expected_doublet <- 0.069
  } else if (val >= 9001) {
    expected_doublet <- 0.076
  }
  
  # Return the expected doublet value
  return(expected_doublet)
}

#' Doublet Removal function
#'
#' Fully runs the doublet removal protocol through DoubletFinder() based on variable from assign_ExpectedDoublet().
#' Returns Seurat object with TRUE/FALSE assignment to duplicate cells in META data.
#' 
#' @param seurat_list List of Seurat objects.
#' @export

preprocess_DoubletRemoval <- function(seurat_list) {
  
  processed_list <- list()
  object_list <- names(seurat_list)
  total_datasets <- length(object_list)
  
  send_message(paste("Starting doublet removal for", total_datasets, "datasets..."))
  
  for (i in seq_along(object_list)) {
    dataset <- object_list[i]
    
    send_message(paste("Processing dataset", i, "of", total_datasets, ":", dataset))
    
    seurat_obj <- seurat_list[[dataset]]
    
    if(!("seurat_clusters" %in% colnames(seurat_obj@meta.data))) {
      send_message(paste("  - Clustering", dataset, "for doublet detection..."))
      seurat_obj <- Seurat::NormalizeData(seurat_obj)
      seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
      seurat_obj <- Seurat::ScaleData(seurat_obj)
      seurat_obj <- Seurat::RunPCA(seurat_obj)
      seurat_obj <- Seurat::FindNeighbors(seurat_obj)
      seurat_obj <- Seurat::FindClusters(seurat_obj)
      seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)
    }
    
    send_message(paste("  - Running parameter sweep for", dataset, "..."))
    sweep.res.list_obj <- DoubletFinder::paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
    sweep.stats_obj <- DoubletFinder::summarizeSweep(sweep.res.list_obj, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats_obj)
    
    val <- length(seurat_obj@assays$RNA@layers$counts@p)
    expected_doublet <- assign_ExpectedDoublet(val)
    send_message(paste("  - Expected doublet rate for", dataset, ":", expected_doublet))
    
    send_message(paste("  - Computing homotypic proportion for", dataset, "..."))
    annotations <- seurat_obj@meta.data$seurat_clusters
    homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
    nExp_poi <- round(expected_doublet*nrow(seurat_obj@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    send_message(paste("  - Running DoubletFinder on", dataset, "..."))
    seurat_obj <- DoubletFinder::doubletFinder(seurat_obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
    
    # Create the plot
    send_message(paste("  - Creating diagnostic plot for", dataset, "..."))
    bcmvn$pK <- as.numeric(bcmvn$pK)
    range_x <- range(bcmvn$pK)
    positions <- seq(range_x[1], range_x[2], length.out = 6)
    x_at_max <- which.max(bcmvn$BCmetric)
    
    plot_obj <- ggplot2::ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
      ggtitle(dataset) +
      geom_point() +
      geom_line() +
      geom_vline(aes(xintercept = x_at_max), linetype="dashed", color = "black") +
      geom_vline(xintercept = positions, color="grey50", linetype="solid") +
      scale_x_continuous(breaks = c(min(bcmvn$pK), x_at_max, max(bcmvn$pK))) +
      theme_minimal() + 
      theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    
    tryCatch({
      output_filename <- paste0(output_folder, "/Preprocessing/DoubletRemoval_", dataset, ".png")
      dir.create(dirname(output_filename), showWarnings = FALSE, recursive = TRUE)
      ggsave(filename = output_filename, plot = plot_obj, width = 6, height = 4, bg = 'white')
      send_message(paste("  - Saved plot to:", output_filename))
    }, error = function(e) {
      send_message(paste("  - Warning: Could not save plot:", e$message))
    })
    
    # Get dimensions before doublet removal
    nPreDoublet <- ncol(seurat_obj)
    
    # Find the DoubletFinder classification column
    df_class_col <- grep("^DF\\.classification", colnames(seurat_obj@meta.data), value = TRUE)
    
    if (length(df_class_col) == 0) {
      send_message(paste("  - Warning: No DoubletFinder classification column found for", dataset))
      processed_list[[dataset]] <- seurat_obj
      next
    }
    
    # Filter doublets
    is_singlet <- seurat_obj@meta.data[[df_class_col]] == "Singlet"
    cells_to_keep <- colnames(seurat_obj)[is_singlet]
    
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    
    # Get dimensions after doublet removal
    nPostDoublet <- ncol(seurat_obj)
    doublets_removed <- nPreDoublet - nPostDoublet
    doublet_percentage <- round(100 * doublets_removed / nPreDoublet, 2)
    
    send_message(sprintf("  - Dataset %s: Removed %s doublet cells (%.2f%%). %s cells remaining.", 
                         dataset, doublets_removed, doublet_percentage, nPostDoublet))
    
    processed_list[[dataset]] <- seurat_obj
  }
  
  send_message(paste("Doublet removal complete for all", total_datasets, "datasets!"))
  return(processed_list)
}

#' Removes ambient RNA from cells using DecontX pipeline
#'
#' @param seurat_list List of Seurat objects
#' @export

preprocess_AmbientRNA <- function(seurat_list) {
  
  processed_list <- list()
  object_list <- names(seurat_list)
  total_datasets <- length(object_list)
  
  send_message(paste("Starting ambient RNA removal for", total_datasets, "datasets..."))
  
  for (i in seq_along(object_list)) {
    dataset <- object_list[i]
    seurat_obj <- seurat_list[[dataset]]
    
    send_message(paste("Processing dataset", i, "of", total_datasets, ":", dataset))
    
    if(!("seurat_clusters" %in% colnames(seurat_obj@meta.data))) {
      send_message(paste("  - Clustering", dataset, "for ambient RNA detection..."))
      seurat_obj <- Seurat::NormalizeData(seurat_obj)
      seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
      seurat_obj <- Seurat::ScaleData(seurat_obj)
      seurat_obj <- Seurat::RunPCA(seurat_obj)
      seurat_obj <- Seurat::FindNeighbors(seurat_obj)
      seurat_obj <- Seurat::FindClusters(seurat_obj)
      seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)
    }
    
    counts <- Seurat::GetAssayData(object = seurat_obj, layer = "counts")
    clusters <- Idents(seurat_obj) %>% as.numeric()
    
    # Add debugging information
    send_message(sprintf("  - Original matrix: %d genes x %d cells, class: %s", 
                         nrow(counts), ncol(counts), class(counts)[1]))
    
    # SPECIFIC FIX: Check for and handle the matrix dimension issue
    # Use Matrix::rowSums for sparse matrices to avoid dimension issues
    if (is(counts, "sparseMatrix") || is(counts, "dgCMatrix")) {
      expressed_genes <- Matrix::rowSums(counts) > 0
    } else {
      expressed_genes <- rowSums(counts) > 0
    }
    
    send_message(sprintf("  - Found %d genes with expression > 0", sum(expressed_genes)))
    
    # CRITICAL: Use drop = FALSE to prevent matrix->vector conversion
    x <- counts[expressed_genes, , drop = FALSE]
    
    # Additional safety check
    if (is.null(dim(x)) || length(dim(x)) < 2) {
      send_message(sprintf("  - ERROR: Matrix lost dimensions after filtering. Class: %s, Dimensions: %s", 
                           class(x)[1], paste(dim(x), collapse = "x")))
      send_message(paste("  - Skipping", dataset, "due to dimension error"))
      processed_list[[dataset]] <- seurat_obj
      next
    }
    
    cell_count <- ncol(x)
    gene_count <- nrow(x)
    send_message(sprintf("  - Running decontX on %s cells with %s non-zero genes...", cell_count, gene_count))
    
    # Additional debugging - let's see exactly what we're passing to decontX
    send_message(sprintf("  - Matrix x: class=%s, dims=%s, sparse=%s", 
                         class(x)[1], paste(dim(x), collapse = "x"), is(x, "sparseMatrix")))
    send_message(sprintf("  - Clusters: length=%d, class=%s, range=%s", 
                         length(clusters), class(clusters), paste(range(clusters), collapse = "-")))
    
    decon <- tryCatch({
      decontX::decontX(x, z=clusters, verbose=FALSE)
    }, error = function(e) {
      send_message(paste("  - Error in decontX:", e$message))
      send_message(paste("  - Full error details:", toString(e)))
      
      # Try without clusters as fallback
      send_message(paste("  - Attempting decontX without predefined clusters..."))
      tryCatch({
        decontX::decontX(x, z=NULL, verbose=FALSE)
      }, error = function(e2) {
        send_message(paste("  - Fallback also failed:", e2$message))
        return(NULL)
      })
    })
    
    if (is.null(decon)) {
      send_message(paste("  - Skipping", dataset, "due to decontX error"))
      processed_list[[dataset]] <- seurat_obj
      next
    }
    
    send_message(paste("  - Applying decontX results to Seurat object for", dataset))
    # Save desired information back to Seurat Object
    newCounts <- decon$decontXcounts
    
    # Add back unexpressed genes and sort according to original counts
    # CRITICAL: Use drop = FALSE here too
    unexpressed_counts <- counts[!expressed_genes, , drop = FALSE]
    newCounts <- rbind(newCounts, unexpressed_counts)[rownames(counts), , drop = FALSE]
    
    seurat_obj[["RNA"]]@layers$counts <- as(round(newCounts), "sparseMatrix")
    seurat_obj$estConp <- decon$contamination
    
    seurat_obj <- subset(seurat_obj)
    
    mean_contamination <- round(mean(decon$contamination) * 100, 2)
    send_message(sprintf("  - Average contamination in %s: %.2f%%", dataset, mean_contamination))
    
    processed_list[[dataset]] <- seurat_obj
  }
  
  send_message(paste("Ambient RNA removal complete for all", total_datasets, "datasets!"))
  return(processed_list)
}

#' Removes cells with deficient RNA after Ambient RNA removal
#'
#' @param seurat_list List of Seurat objects
#' @export

remove_lowRNA <- function(seurat_list) {
  
  processed_list <- list()
  object_list <- names(seurat_list)
  total_datasets <- length(object_list)
  
  send_message(paste("Starting low RNA cell filtering for", total_datasets, "datasets..."))
  
  for (i in seq_along(object_list)) {
    dataset <- object_list[i]
    
    # Access the dataset using the name
    seurat_obj <- seurat_list[[dataset]]
    
    # Get dimensions before filtering
    nPreDecon <- ncol(seurat_obj)
    
    send_message(paste("Filtering dataset", i, "of", total_datasets, ":", dataset))
    
    # Filter cells with low RNA
    seurat_obj <- subset(seurat_obj, subset = c(nFeature_RNA > 100 & nCount_RNA > 500))
    
    # Get dimensions after filtering
    nPostDecon <- ncol(seurat_obj)
    cells_removed <- nPreDecon - nPostDecon
    cells_removed_pct <- round(100 * cells_removed / nPreDecon, 2)
    
    send_message(sprintf("  - Dataset %s: Removed %s cells with too little RNA (%.2f%%). %s cells remaining.", 
                         dataset, cells_removed, cells_removed_pct, nPostDecon))
    
    # Store in processed list
    processed_list[[dataset]] <- seurat_obj
  }
  
  send_message(paste("Low RNA cell filtering complete for all", total_datasets, "datasets!"))
  return(processed_list)
}

#' Function to completely remove genes from a list of Seurat objects
#'
#' @param seurat_list List of Seurat objects to filter.
#' @param genes_to_remove vector of gene names to remove completely. If "remove_ensembl" is passed, will automatically detect ENSG genes.
#' @param send_message function to send messages to the UI (optional).
#' @return filtered list of Seurat objects.
#' @export
#'
remove_GenesFromSeuratList <- function(seurat_list, genes_to_remove, send_message = NULL) {
  
  processed_list <- list()
  object_list <- names(seurat_list)
  total_datasets <- length(object_list)
  
  send_message(sprintf("Removing genes from %d datasets...", total_datasets))
  
  for (i in seq_along(object_list)) {
    dataset_name <- object_list[i]
    seurat_obj <- seurat_list[[dataset_name]]
    
    # Get current gene names
    current_genes <- rownames(seurat_obj)
    
    # Check if we need to automatically detect ENSG genes
    if (length(genes_to_remove) == 1 && genes_to_remove == "remove_ensembl") {
      # Find all genes matching ENSG pattern (case insensitive)
      # Patterns like: ENSG00000304370.1, ensg00000351123.2, ENSG00000123456, etc.
      ensg_genes <- current_genes[grepl("^ensg[0-9]+", current_genes, ignore.case = TRUE)]
      genes_to_remove_existing <- ensg_genes
      
      # send_message(sprintf("  - Dataset %s: Found %d ENSG genes to remove", dataset_name, length(ensg_genes)))
      
      # if (length(ensg_genes) > 0) {
      #   # Show some examples of ENSG genes found
      #   example_genes <- head(ensg_genes, 5)
      #   send_message(sprintf("  - Examples: %s", paste(example_genes, collapse = ", ")))
      # }
      
    } else {
      # Original behavior: use provided gene list
      genes_to_remove_existing <- intersect(genes_to_remove, current_genes)
      genes_not_found <- setdiff(genes_to_remove, current_genes)
      
      if (length(genes_not_found) > 0) {
        send_message(sprintf("  - Dataset %s: %d genes not found in object.", dataset_name, length(genes_not_found)))
      }
    }
    
    if (length(genes_to_remove_existing) == 0) {
      send_message(sprintf("  - Dataset %s: No specified genes found. Keeping original object.", dataset_name))
      processed_list[[dataset_name]] <- seurat_obj
      next
    }
    
    # Get genes to keep
    genes_to_keep <- setdiff(current_genes, genes_to_remove_existing)
    
    # send_message(sprintf("  - Dataset %s: Removing %d genes, keeping %d genes out of original %d genes.", 
    #                      dataset_name, length(genes_to_remove_existing), length(genes_to_keep), length(current_genes)))
    
    # Filter the Seurat object
    tryCatch({
      # Subset the Seurat object to keep only the desired genes
      filtered_seurat <- seurat_obj[genes_to_keep, ]
      processed_list[[dataset_name]] <- filtered_seurat
      
    }, error = function(e) {
      send_message(sprintf("  - Error removing genes from %s: %s", dataset_name, e$message))
      send_message(sprintf("  - Keeping original %s object.", dataset_name))
      processed_list[[dataset_name]] <- seurat_obj
    })
  }
  
  send_message("Gene removal completed for all datasets.")

  


  return(processed_list)
}