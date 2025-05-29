#' Variable gene function selection enabling.
#'
#' Generates variable genes from matrix
#'
#' @return Variable list of genes
#' @param mat raw count matrix from Seurat object.
#' @param nvar number of variable genes to be found.
#' @param blacklist list of genes to exclude from variable gene analysis.
#' @export
#' 
generate_VariableGenes <- function(mat, nvar = 2000, blacklist = NULL){
  # Safety checks
  if (is.null(mat) || length(dim(mat)) != 2) {
    stop("Input matrix must be a 2-dimensional matrix or sparse matrix")
  }
  
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    warning("Input matrix is empty")
    return(character(0))
  }
  
  # Ensure rownames exist
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("Gene_", 1:nrow(mat))
    message("Matrix had no rownames. Generated generic names.")
  }
  
  # Get the top nvar variable genes present in mat (a gene x sample/cell matrix)
  # If blacklist is present, do not return any genes in the blacklist
  if(!is.null(blacklist) && length(blacklist) > 0){
    # Ensure blacklist is character vector
    blacklist <- as.character(blacklist)
    
    ncount <- nrow(mat)
    # Use proper logical indexing with %in%
    genes_to_keep <- !(rownames(mat) %in% blacklist)
    mat <- mat[genes_to_keep, , drop = FALSE]
    
    message(sprintf("Removed %s genes overlapping blacklist prior to selecting variable genes...", ncount - nrow(mat)))
    
    # Check if we still have genes left
    if (nrow(mat) == 0) {
      warning("All genes were filtered out by blacklist")
      return(character(0))
    }
  }
  
  # Calculate variance based on matrix type
  tryCatch({
    if(is(mat, "sparseMatrix")){
      if (!requireNamespace("sparseMatrixStats", quietly = TRUE)) {
        message("sparseMatrixStats not available, converting to dense matrix")
        gene_vars <- matrixStats::rowVars(as.matrix(mat))
      } else {
        gene_vars <- sparseMatrixStats::rowVars(mat)
      }
    } else {
      if (!requireNamespace("matrixStats", quietly = TRUE)) {
        gene_vars <- apply(mat, 1, var)
      } else {
        gene_vars <- matrixStats::rowVars(mat)
      }
    }
    
    # Handle potential NA values in variance calculation
    gene_vars[is.na(gene_vars)] <- 0
    
    # Get top variable genes
    n_genes_to_select <- min(nvar, length(gene_vars))
    top_indices <- head(order(gene_vars, decreasing = TRUE), n_genes_to_select)
    varGenes <- rownames(mat)[top_indices]
    
    message(sprintf("Selected %d variable genes out of %d total genes", 
                    length(varGenes), nrow(mat)))
    
    return(varGenes)
    
  }, error = function(e) {
    message(paste("Error in variance calculation:", e$message))
    message("Returning first n genes as fallback")
    
    # Fallback: return first n genes
    n_genes_to_select <- min(nvar, nrow(mat))
    return(rownames(mat)[1:n_genes_to_select])
  })
}

#' Sparse Log of X function. Log normalization on sparse matrix, returns sparse matrix.
#'
#' Work around for staying with sparse matrix format in Seurat object to reduce object size.
#' 
#' @return Logarithmic sparse matrix
#' @param spmat raw count matrix from Seurat object in sparse format.
#' @param logtype c("log", "log2", "log10"), standardized to log2.
#' @param scaleFactor scaling factor, standardized to 10^4.
#' @export
#' 
run_SparseLogX <- function(spmat, logtype="log2", scale=FALSE, scaleFactor=10^4){
  stopifnot(any(logtype == c("log", "log2", "log10")))
  
  # Check if input is a sparse matrix
  if(!is(spmat, "sparseMatrix")){
    message("Warning: Converting dense matrix to sparse format for efficiency")
    spmat <- as(spmat, "sparseMatrix")
  }
  
  # Handle empty matrix case
  if(length(spmat@x) == 0) {
    message("Empty matrix detected, returning as is")
    return(spmat)
  }
  
  # Scale matrix if requested - with safeguards against zero column sums
  if(scale == TRUE){
    colSums_vals <- Matrix::colSums(spmat)
    
    # Check for zero column sums to avoid division by zero
    if(any(colSums_vals == 0)) {
      message(sprintf("Note: %d columns with zero sums detected; adding small pseudocount", sum(colSums_vals == 0)))
      colSums_vals[colSums_vals == 0] <- 1e-10
    }
    
    # Scale columns more efficiently with direct manipulation of sparse matrix slots
    # This avoids creating intermediate matrices
    i <- spmat@i
    p <- spmat@p
    x <- spmat@x
    
    # Apply scaling to each column based on column sums
    for(col in 1:ncol(spmat)) {
      idx_start <- p[col] + 1
      idx_end <- p[col + 1]
      
      if(idx_start <= idx_end) {
        col_indices <- idx_start:idx_end
        x[col_indices] <- (x[col_indices] / colSums_vals[col]) * scaleFactor
      }
    }
    
    # Reconstruct sparse matrix with scaled values
    spmat@x <- x
    
    # Force garbage collection to free memory
    invisible(gc())
  }
  
  # Log transform the sparse matrix in place
  # Operating directly on the non-zero values (x slot)
  if(logtype == "log") {
    spmat@x <- log(spmat@x + 1)
  } else if(logtype == "log2") {
    spmat@x <- log2(spmat@x + 1) 
  } else {
    spmat@x <- log10(spmat@x + 1)
  }
  
  # Force garbage collection again
  invisible(gc())
  
  return(spmat)
}

#' Function to generate a gene blacklist based on user selections
#'
#' @param count_matrix raw count matrix from Seurat object.
#' @param selected_blacklist vector of selected blacklist options from UI.
#' @param send_message function to send messages to the UI (optional).
#' @return vector of gene names to blacklist.
#' @export
#' 
generate_GeneBlacklist <- function(count_matrix, selected_blacklist, send_message = NULL) {
  # Initialize empty blacklist
  blacklist.genes <- character(0)
  
  # If no message function provided, create a dummy one
  if (is.null(send_message)) {
    send_message <- function(msg) {
      message(msg)
    }
  }
  
  pkg_available <- c(
    "TxDb" = requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE),
    "GenomicFeatures" = requireNamespace("GenomicFeatures", quietly = TRUE),
    "org.Hs.eg.db" = requireNamespace("org.Hs.eg.db", quietly = TRUE),
    "AnnotationDbi" = requireNamespace("AnnotationDbi", quietly = TRUE),
    "GenomicRanges" = requireNamespace("GenomicRanges", quietly = TRUE)
  )
  
  send_message(paste("Package availability:", paste(names(pkg_available), pkg_available, sep="=", collapse=", ")))
  
  
  # Safety check: ensure count_matrix has rownames
  if (is.null(rownames(count_matrix))) {
    send_message("Warning: Count matrix has no rownames. Cannot generate gene blacklist.")
    return(character(0))
  }
  
  # Convert rownames to character vector to ensure proper matching
  gene_names <- as.character(rownames(count_matrix))
  
  # Ensure selected_blacklist is properly formatted
  if (is.null(selected_blacklist) || length(selected_blacklist) == 0) {
    send_message("No genes selected for blacklisting.")
    return(character(0))
  }
  
  # Convert to character vector if not already
  selected_blacklist <- as.character(selected_blacklist)
  
  send_message(paste("Processing blacklist options:", paste(selected_blacklist, collapse = ", ")))
  
  # Add mitochondrial genes if selected
  if("blacklist_mitogenes" %in% selected_blacklist) {
    tryCatch({
      mt.genes <- grep(pattern = "^MT-", x = gene_names, value = TRUE)
      if (length(mt.genes) > 0) {
        blacklist.genes <- c(blacklist.genes, mt.genes)
        send_message(sprintf("Added %d mitochondrial genes to blacklist.", length(mt.genes)))
      } else {
        send_message("No mitochondrial genes found with pattern ^MT-")
      }
    }, error = function(e) {
      send_message(paste("Error finding mitochondrial genes:", e$message))
    })
  }
  
  # Add sex chromosome genes if selected
  if("blacklist_sexgenes" %in% selected_blacklist) {
    tryCatch({
      # Load required packages explicitly
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
        send_message("TxDb.Hsapiens.UCSC.hg38.knownGene package not available. Installing...")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
      }
      
      if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
        send_message("GenomicFeatures package not available. Installing...")
        BiocManager::install("GenomicFeatures")
      }
      
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        send_message("org.Hs.eg.db package not available. Installing...")
        BiocManager::install("org.Hs.eg.db")
      }
      
      # Load the packages
      library(TxDb.Hsapiens.UCSC.hg38.knownGene)
      library(GenomicFeatures)
      library(org.Hs.eg.db)
      library(AnnotationDbi)
      
      # Extract sex chromosome genes
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      geneGR <- GenomicFeatures::genes(txdb)
      
      seqnames_clean <- as.vector(GenomicRanges::seqnames(geneGR))
      sex_mask <- seqnames_clean %in% c("chrY", "chrX")
      sexGenesGR <- geneGR[sex_mask]
      
      
      if (length(sexGenesGR) > 0) {
        matchedGeneSymbols <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                                    keys = sexGenesGR$gene_id,
                                                    columns = c("ENTREZID", "SYMBOL"),
                                                    keytype = "ENTREZID")
        
        # Filter out NA symbols and ensure character vector
        sexChr.genes <- as.character(matchedGeneSymbols$SYMBOL)
        sexChr.genes <- sexChr.genes[!is.na(sexChr.genes)]
        
        if (length(sexChr.genes) > 0) {
          blacklist.genes <- c(blacklist.genes, sexChr.genes)
          send_message(sprintf("Added %d X/Y chromosomal genes to blacklist.", length(sexChr.genes)))
        } else {
          send_message("No valid sex chromosome gene symbols found.")
        }
      } else {
        send_message("No sex chromosome genes found in annotation.")
      }
    }, error = function(e) {
      send_message(paste("Error finding sex chromosome genes:", e$message))
      send_message("This may be due to missing annotation packages. Try installing: BiocManager::install(c('TxDb.Hsapiens.UCSC.hg38.knownGene', 'GenomicFeatures', 'org.Hs.eg.db'))")
    })
  }
  
  # Add ribosomal genes if selected
  if("blacklist_rbgenes" %in% selected_blacklist) {
    tryCatch({
      RPS.genes <- grep(pattern = "^RPS", x = gene_names, value = TRUE)
      RPL.genes <- grep(pattern = "^RPL", x = gene_names, value = TRUE)
      
      ribo_genes <- c(RPS.genes, RPL.genes)
      if (length(ribo_genes) > 0) {
        blacklist.genes <- c(blacklist.genes, ribo_genes)
        send_message(sprintf("Added %d ribosomal genes to blacklist (RPS: %d, RPL: %d).", 
                             length(ribo_genes), length(RPS.genes), length(RPL.genes)))
      } else {
        send_message("No ribosomal genes found with patterns ^RPS or ^RPL")
      }
    }, error = function(e) {
      send_message(paste("Error finding ribosomal genes:", e$message))
    })
  }
  
  if("blacklist_ensembl" %in% selected_blacklist) {
    tryCatch({
      # Pattern to match Ensembl gene IDs (ENSG followed by numbers and optional version)
      ensembl.genes <- grep(pattern = "^ENSG[0-9]+", x = gene_names, value = TRUE)
      
      if (length(ensembl.genes) > 0) {
        blacklist.genes <- c(blacklist.genes, ensembl.genes)
        send_message(sprintf("Added %d Ensembl genes to blacklist.", length(ensembl.genes)))
      } else {
        send_message("No Ensembl genes found with pattern ^ENSG")
      }
    }, error = function(e) {
      send_message(paste("Error finding Ensembl genes:", e$message))
    })
  }
  
  # Remove duplicates and ensure character vector
  if (length(blacklist.genes) > 0) {
    blacklist.genes <- unique(as.character(blacklist.genes))
    
    # Filter to only genes that actually exist in the count matrix
    existing_genes <- blacklist.genes[blacklist.genes %in% gene_names]
    missing_genes <- setdiff(blacklist.genes, existing_genes)
    
    if (length(missing_genes) > 0) {
      send_message(sprintf("Warning: %d blacklisted genes not found in count matrix.", length(missing_genes)))
    }
    
    blacklist.genes <- existing_genes
    send_message(sprintf("Final blacklist contains %d genes present in the dataset.", length(blacklist.genes)))
  } else {
    send_message("No valid genes added to blacklist.")
  }
  
  return(blacklist.genes)
}

#' Function to form row and column sums and means on count matrix from Seurat object.
#'
#' Summarized group of count matrix for LSI
#' 
#' @param mat raw count matrix.
#' @param groups which groups to summarize.
#' @param na.rm whether to remove any NA from matrix.
#' @param sparse whether to use sparse matrix conversion.
#' @export
#'
generate_GroupSums <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  # Input validation
  if (is.null(groups)) {
    stop("groups parameter cannot be NULL")
  }
  
  if (length(groups) != ncol(mat)) {
    stop(sprintf("Length of groups (%d) must equal number of columns in matrix (%d)", 
                 length(groups), ncol(mat)))
  }
  
  # Ensure groups is a proper vector
  groups <- as.character(groups)
  unique_groups <- unique(groups)
  
  if (length(unique_groups) == 0) {
    stop("No unique groups found")
  }
  
  # Calculate group sums
  tryCatch({
    gm <- lapply(unique_groups, function(x) {
      group_indices <- which(groups == x)
      if (length(group_indices) == 0) {
        return(rep(0, nrow(mat)))
      }
      
      # Subset matrix for this group
      group_mat <- mat[, group_indices, drop = FALSE]
      
      if (sparse || is(mat, "sparseMatrix")) {
        Matrix::rowSums(group_mat, na.rm = na.rm)
      } else {
        rowSums(group_mat, na.rm = na.rm)
      }
    })
    
    # Convert list to matrix
    gm <- do.call(cbind, gm)
    colnames(gm) <- unique_groups
    
    # Ensure rownames are preserved
    if (!is.null(rownames(mat))) {
      rownames(gm) <- rownames(mat)
    }
    
    return(gm)
    
  }, error = function(e) {
    message(paste("Error in generate_GroupSums:", e$message))
    # Return a fallback matrix
    fallback_mat <- matrix(0, nrow = nrow(mat), ncol = length(unique_groups))
    colnames(fallback_mat) <- unique_groups
    if (!is.null(rownames(mat))) {
      rownames(fallback_mat) <- rownames(mat)
    }
    return(fallback_mat)
  })
}

#' Function to run iterative LSI on a Seurat object.
#'
#' Runs the entire interative LSI pipeline on a seurat object.
#'
#' @param seurat_obj Seurat object.
#' @param rawCounts count matrix from Seurat object.
#' @param blacklist.genes vector of gene names to blacklist.
#' @param nVarGenes number of variable genes to use.
#' @param resolution vector of resolution values for clustering.
#' @param harmonize which dimensions to harmonize.
#' @param nPCs number of principal components to use.
#' @param covariates vector of covariates to use for harmony integration.
#' @param umapNeighbors number of neighbors for UMAP.
#' @param umapMinDist minimum distance for UMAP.
#' @param umapDistMetric distance metric for UMAP.
#' @param send_message function to send messages to the UI.
#' @return list containing updated Seurat object and LSI results.
#' @export
#'
process_LSI <- function(seurat_obj, 
                        rawCounts = NULL,
                        blacklist.genes = NULL,
                        nVarGenes = 2000, 
                        resolution = c(0.2, 0.4, 0.8), 
                        nPCs = 30, 
                        harmonize = c(1, 2), 
                        covariates = c("orig.ident"), 
                        umapNeighbors = 50,
                        umapMinDist = 0.5,
                        umapDistMetric = "cosine",
                        send_message = NULL) {
  
  if (is.null(send_message)) {
    send_message <- function(msg) { message(msg) }
  }
  
  # Perform log normalization
  send_message("Performing log normalization...")
  log2CP10k <- run_SparseLogX(spmat = rawCounts, logtype="log2", scale=TRUE, scaleFactor=10^4)
  seurat_obj <- SeuratObject::SetAssayData(object=seurat_obj, layer="data", new.data=log2CP10k)
  
  # Create list to store LSI results
  lsiOut <- list()
  
  # Initialize clusters here, just like in the working script
  clusters <- NULL
  
  # Run iterative LSI
  send_message("Running iterative LSI...")
  set.seed(1)
  
  for(i in seq_along(resolution)){
    send_message(sprintf("Starting LSI round %d with resolution %s...", i, resolution[i]))
    
    # If first round, compute variable genes on raw data first
    if(i == 1){
      send_message(sprintf("Identifying top %s variable genes among all cells...", nVarGenes))
      varGenes <- generate_VariableGenes(log2CP10k, nvar=nVarGenes, blacklist=blacklist.genes)
    } else {
      # For remaining rounds, calculate variable genes using previous clusters
      # Try-catch to handle potential errors
      tryCatch({
        clusterMat <- edgeR::cpm(generate_GroupSums(rawCounts, clusters, sparse=TRUE), log=TRUE, prior.count=3)
        send_message(sprintf("Identifying top %s variable genes from round %s LSI...", nVarGenes, i-1))
        varGenes <- generate_VariableGenes(clusterMat, nvar=nVarGenes, blacklist=blacklist.genes)
      }, error = function(e) {
        send_message(paste("Error in cluster-based gene selection:", e$message))
        send_message("Falling back to log-normalized data for gene selection")
        varGenes <- generate_VariableGenes(log2CP10k, nvar=nVarGenes, blacklist=blacklist.genes)
      })
    }
    
    # Run LSI
    send_message(sprintf("Running LSI on %d variable genes...", length(varGenes)))
    LSIi <- LSI_function(rawCounts[varGenes,], nComponents = max(nPCs), binarize = FALSE)
    
    # Harmonize if necessary - using same condition as working script
    if(i %in% harmonize && length(covariates) > 0){
      send_message(sprintf("Harmonizing LSI SVD PCs for round %s...", i))
      tryCatch({
        harmonized_pcs <- harmony::RunHarmony(
          data_mat  = LSIi$matSVD,
          meta_data = seurat_obj@meta.data,
          vars_use  = covariates,
          do_pca    = FALSE,
          # Add some conservative parameters for stability
          nclust = min(20, max(3, floor(ncol(LSIi$matSVD)/5)))
        )
        LSIi$matSVD <- harmonized_pcs
      }, error = function(e) {
        send_message(paste("Harmony integration failed:", e$message))
        send_message("Using original SVD results without Harmony integration.")
      })
    }
    
    # Create dimension reduction
    reducName <- paste0("LSI_iter", i)
    send_message(sprintf("Creating dimension reduction object for %s...", reducName))
    seurat_obj[[reducName]] <- SeuratObject::CreateDimReducObject(
      embeddings = LSIi$matSVD, 
      key = sprintf("LSI%s_", i), 
      assay = "RNA"
    )
    
    # Find neighbors and clusters
    send_message("Finding nearest neighbors...")
    seurat_obj <- Seurat::FindNeighbors(
      object = seurat_obj, 
      reduction = reducName, 
      dims = 1:nPCs, 
      force.recalc = TRUE
    )
    
    send_message(sprintf("Clustering with resolution %s...", resolution[i]))
    seurat_obj <- Seurat::FindClusters(
      object = seurat_obj, 
      resolution = resolution[i]
    )
    
    # Store cluster identities
    clusters <- Idents(seurat_obj)
    send_message(sprintf("Found %d clusters in iteration %d", length(unique(clusters)), i))
    
    # Store information
    lsiOut[[reducName]] <- list(
      lsiMat = LSIi$matSVD,
      svd = LSIi$svd,
      varFeatures = varGenes, 
      clusters = clusters
    )
    
    # Force garbage collection
    invisible(gc())
    
    send_message(sprintf("Completed LSI round %d", i))
  }
  
  # Run UMAP using the last LSI iteration
  lastIteration <- paste0("LSI_iter", length(resolution))
  send_message(sprintf("Running UMAP on final LSI iteration (%s)...", lastIteration))
  
  # Run UMAP on the final iteration
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = lastIteration, # Use the variable with the last iteration name
    dims = 1:nPCs,
    n.neighbors = umapNeighbors,
    min.dist = umapMinDist,
    metric = umapDistMetric
  )
  
  # Return results
  return(list(
    seurat_obj = seurat_obj,
    lsiOut = lsiOut
  ))
}

#' Function to run LSI on a Seurat object.
#'
#' @param mat raw count matrix from Seurat object.
#' @param nComponents number of right singular vectors to estimate.
#' @param binarize whether to binarize dgCMatrix class.
#' @return list containing LSI results.
#' @export
LSI_function <- function(mat, nComponents, binarize = FALSE){
  
  if(binarize){
    message("Binarizing matrix...")
    # The 'x' slot of the dgCMatrix class contains the non-zero elements of the matrix
    mat@x[mat@x > 0] <- 1
  }
  
  #Calculate RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  
  # Check for zero row sums and remove those rows
  if(any(rowSm == 0)) {
    message(sprintf("Removing %d features with zero counts...", sum(rowSm == 0)))
    mat <- mat[rowSm > 0, , drop = FALSE]
    rowSm <- rowSm[rowSm > 0]
  }
  
  # Calculate TF-IDF
  message(sprintf("Calculating TF-IDF with %s features (terms) and %s cells (documents)...", nrow(mat), ncol(mat)))
  start <- Sys.time()
  scaleTo <- 10^4
  
  message("Calculating TF using Matrix transpose...")
  
  # Add a small epsilon to avoid division by zero
  # Use Matrix::t() instead of base::t()
  tf <- Matrix::t(Matrix::t(mat) / (colSm + 1e-10))
  
  message("TF calculation successful")
  
  idf <- as(ncol(mat) / rowSm, "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% tf
  
  # Log transform TF-IDF
  tfidf <- run_SparseLogX(tfidf, logtype="log", scale=TRUE, scaleFactor=scaleTo)
  
  # Replace any NaN or Inf values
  if(any(is.na(tfidf@x) | is.infinite(tfidf@x))) {
    message("Warning: NaN or Inf values detected in TF-IDF matrix. Replacing with zeros...")
    tfidf@x[is.na(tfidf@x) | is.infinite(tfidf@x)] <- 0
  }
  
  # Clean up
  rm(tf)
  rm(idf)
  invisible(gc())
  
  # Calculate SVD for LSI with error handling
  message("Calculating SVD for LSI...")
  tryCatch({
    # Check if matrix is suitable for SVD
    if(any(is.na(tfidf@x)) || any(is.infinite(tfidf@x))) {
      stop("Matrix contains NA or Inf values")
    }
    
    svd <- irlba::irlba(tfidf, nv=nComponents, nu=nComponents)
    svdDiag <- Matrix::diag(x=svd$d)
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(svd$v))  # Use Matrix::t() here too
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC", seq_len(ncol(matSVD)))
    
    message(sprintf("LSI complete: %s minutes", round(difftime(Sys.time(), start, units="mins"), 3)))
  }, error = function(e) {
    message(paste("Error in SVD calculation:", e$message))
    message("Using fallback SVD method...")
    
    # Convert sparse matrix to regular matrix for prcomp
    dense_tfidf <- as.matrix(tfidf)
    
    # Replace any remaining problematic values
    dense_tfidf[is.na(dense_tfidf)] <- 0
    dense_tfidf[is.infinite(dense_tfidf)] <- 0
    
    # Use standard PCA as fallback with reduced rank
    pca_result <- prcomp(t(dense_tfidf), center = TRUE, scale. = FALSE, 
                         rank. = min(nComponents, min(dim(dense_tfidf))-1))
    
    # Create SVD-like output
    svd <- list(
      d = pca_result$sdev[1:min(nComponents, length(pca_result$sdev))],
      u = matrix(0, nrow=nrow(dense_tfidf), ncol=length(pca_result$sdev)),  # Placeholder
      v = pca_result$rotation
    )
    
    # Create output matrix similar to original function
    matSVD <- pca_result$x
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC", seq_len(ncol(matSVD)))
    
    message("Fallback SVD completed successfully")
  })
  
  # Return result
  if(is.null(rownames(mat))){
    rownames(mat) <- 1:nrow(mat)
  }
  
  return(
    list(
      matSVD = matSVD, 
      rowSm = rowSm, 
      colSm = colSm, 
      svd = svd, 
      binarize = binarize
    )
  )
}