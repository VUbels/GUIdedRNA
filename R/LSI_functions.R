#' Variable gene function selection enabling.
#'
#' @description Work around for staying with sparse matrix format in seurat object to reduce object size.
#' @returns Variable list of genes
#' @param mat raw count matrix from Seurat object.
#' @param nvar number of variable genes to be found.
#' @param blacklist list of genes to exclude from variable gene analysis.
#' @export
#' 

generate_VariableGenes <- function(mat, nvar = 2000, blacklist = NULL){
  # Get the top nvar variable genes present in mat (a gene x sample/cell matrix)
  # If blacklist is present, do not return any genes in the blacklist
  if(!is.null(blacklist)){
    ncount <- nrow(mat)
    mat <- mat[!rownames(mat) %in% blacklist,]
    message(sprintf("Removed %s genes overlapping blacklist prior to selecting variable genes...", ncount - nrow(mat)))
  }
  if(is(mat, "sparseMatrix")){
    varGenes <- rownames(mat)[head(order(sparseMatrixStats::rowVars(mat), decreasing = TRUE), nvar)]
  }else{
    varGenes <- rownames(mat)[head(order(matrixStats::rowVars(mat), decreasing = TRUE), nvar)]
  }
  return(varGenes)
}

#' Sparse Log of X function. Log normalization on sparse matrix, returns sparse matrix.
#'
#' @description Work around for staying with sparse matrix format in Seurat object to reduce object size.
#' @returns Logarithmic sparse matrix
#' @param spmat raw count matrix from Seurat object in sparse format.
#' @param logtype c("log", "log2", "log10"), standardized to log2.
#' @param scaleFactor scaling factor, standardized to 10^4.
#' @export

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

#' Function to generate a gene blacklist based on user selections.
#'
#' @param count_matrix raw count matrix from Seurat object.
#' @param selected_blacklist vector of selected blacklist options from UI.
#' @param send_message function to send messages to the UI (optional).
#' @return vector of gene names to blacklist.
#' @export
#' 
generate_GeneBlacklist <- function(count_matrix, selected_blacklist, send_message = NULL) {
  # Initialize empty blacklist
  blacklist.genes <- c()
  
  # If no message function provided, create a dummy one
  if (is.null(send_message)) {
    send_message <- function(msg) {
      message(msg)
    }
  }
  
  # If blacklist options are selected
  if(length(selected_blacklist) > 0) {
    # Add mitochondrial genes if selected
    if("blacklist_mitogenes" %in% selected_blacklist) {
      mt.genes <- grep(pattern = "^MT-", x = rownames(count_matrix), value = TRUE)
      blacklist.genes <- c(blacklist.genes, mt.genes)
      send_message("Added mitochondrial genes to blacklist.")
    }
    
    # Add sex chromosome genes if selected
    if("blacklist_sexgenes" %in% selected_blacklist) {
      # Extract sex chromosome genes
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      geneGR <- GenomicFeatures::genes(txdb)
      sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
      matchedGeneSymbols <- select(org.Hs.eg.db,
                                   keys = sexGenesGR$gene_id,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "ENTREZID")
      sexChr.genes <- matchedGeneSymbols$SYMBOL
      blacklist.genes <- c(blacklist.genes, sexChr.genes)
      send_message("Added X/Y chromosomal genes to blacklist.")
    }
    
    # Add ribosomal genes if selected
    if("blacklist_rbgenes" %in% selected_blacklist) {
      RPS.genes <- grep(pattern = "^RPS", x = rownames(count_matrix), value = TRUE)
      RPL.genes <- grep(pattern = "^RPL", x = rownames(count_matrix), value = TRUE)
      blacklist.genes <- c(blacklist.genes, RPS.genes, RPL.genes)
      send_message("Added ribosomal genes to blacklist.")
    }
    
    # Remove duplicates
    blacklist.genes <- unique(blacklist.genes)
    send_message(sprintf("Final blacklist contains %d genes.", length(blacklist.genes)))
  } else {
    send_message("No genes selected for blacklisting.")
  }
  
  return(blacklist.genes)
}

#' Function to form row and column sums and means on count matrix from Seurat object.
#'
#' @param mat raw count matrix.
#' @param groups which groups to summarize.
#' @param na.rm whether to remove any NA from matrix.
#' @param sparse whether to use sparse matrix conversion.
#' @export

generate_GroupSums <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}


#' Function to run iterative LSI on a Seurat object.
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
    mat <- mat[rowSm > 0, ]
    rowSm <- rowSm[rowSm > 0]
  }
  
  # Calculate TF-IDF
  message(sprintf("Calculating TF-IDF with %s features (terms) and %s cells (documents)...", nrow(mat), ncol(mat)))
  start <- Sys.time()
  scaleTo <- 10^4
  
  # Add a small epsilon to avoid division by zero
  tf <- t(t(mat) / (colSm + 1e-10))
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
    matSVD <- t(svdDiag %*% t(svd$v))
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

