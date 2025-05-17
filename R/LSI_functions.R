#' Variable gene function selection enabling.
#'
#' Work around for staying with sparse matrix format in seurat object to reduce object size.
#' 
#' @param mat raw count matrix from Seurat object.
#' @param nvar number of variable genes to be found.
#' @param blacklist list of genes to exclude from variable gene analysis.

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
#' Work around for staying with sparse matrix format in Seurat object to reduce object size.
#' 
#' @param spmat raw count matrix from Seurat object in sparse format.
#' @param logtype c("log", "log2", "log10"), standardized to log2.
#' @param scaleFactor scaling factor, standardized to 10^4.

run_SparseLogX <- function(spmat, logtype="log2", scale=FALSE, scaleFactor=10^4){
  stopifnot(any(logtype == c("log", "log2", "log10")))
  
  if(scale == TRUE){
    spmat <- t(t(spmat)/Matrix::colSums(spmat)) * scaleFactor
  }
  
  if(is(spmat, "sparseMatrix")){
    matsum <- summary(spmat) # Get the sparse matrix summary
    if(logtype == "log"){
      logx <- log(matsum$x + 1) 
    }else if(logtype == "log2"){
      logx <- log2(matsum$x + 1) 
    }else{
      logx <- log10(matsum$x + 1)
    }
    logmat <- sparseMatrix(i = matsum$i, j = matsum$j, # convert back to sparse matrix
                           x = logx, dims = dim(spmat),
                           dimnames = dimnames(spmat))
  }else{
    if(logtype == "log"){
      logmat <- log(spmat + 1) 
    }else if(logtype == "log2"){
      logmat <- log2(spmat + 1) 
    }else{
      logmat <- log10(spmat + 1) 
    }
  }
  return(logmat)
}

#' Function to generate a gene blacklist based on user selections
#'
#' @param count_matrix raw count matrix from Seurat object
#' @param selected_blacklist vector of selected blacklist options from UI
#' @param send_message function to send messages to the UI (optional)
#' @return vector of gene names to blacklist
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

#' Function to run iterative LSI on a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param blacklist.genes vector of gene names to blacklist
#' @param nVarGenes number of variable genes to use
#' @param resolution vector of resolution values for clustering
#' @param nPCs number of principal components to use
#' @param covariates vector of covariates to use for harmony integration
#' @param umapNeighbors number of neighbors for UMAP
#' @param umapMinDist minimum distance for UMAP
#' @param umapDistMetric distance metric for UMAP
#' @param send_message function to send messages to the UI
#' @return list containing updated Seurat object and LSI results
#'
process_LSI <- function(seurat_obj, 
                              blacklist.genes = NULL,
                              nVarGenes = 2000, 
                              resolution = c(0.2, 0.4, 0.8), 
                              nPCs = 30, 
                              covariates = c("orig.ident"),
                              umapNeighbors = 50,
                              umapMinDist = 0.5,
                              umapDistMetric = "cosine",
                              send_message = NULL) {
  
  # If no message function provided, create a dummy one
  if (is.null(send_message)) {
    send_message <- function(msg) {
      message(msg)
    }
  }

  # Perform log normalization
  send_message("Performing log normalization...")
  log2CP10k <- run_SparseLogX(rawCounts, logtype="log2", scale=TRUE, scaleFactor=10^4)
  seurat_obj <- SeuratObject::SetAssayData(object=seurat_obj, layer="data", new.data=log2CP10k)
  
  # Create list to store LSI results
  lsiOut <- list()
  
  # Run iterative LSI
  send_message("Running iterative LSI...")
  set.seed(1)
  
  for(i in seq_along(resolution)){
    send_message(sprintf("Starting LSI round %d with resolution %s...", i, resolution[i]))
    
    # If first round, compute variable genes on raw data first
    if(i == 1){
      send_message(sprintf("Identifying top %s variable genes among all cells...", nVarGenes))
      varGenes <- generate_VariableGenes(log2CP10k, nvar=nVarGenes, blacklist=blacklist.genes)
    }else{
      # For remaining rounds, calculate variable genes using previous clusters
      clusterMat <- edgeR::cpm(groupSums(rawCounts, clusters, sparse=TRUE), log=TRUE, prior.count=3)
      send_message(sprintf("Identifying top %s variable genes from round %s LSI...", nVarGenes, i-1))
      varGenes <- generate_VariableGenes(clusterMat, nvar=nVarGenes, blacklist=blacklist.genes)
    }
    
    # Run LSI and find clusters
    send_message(sprintf("Running LSI on %d variable genes...", length(varGenes)))
    LSIi <- LSI_function(rawCounts[varGenes,], nComponents = max(nPCs), binarize = FALSE)
    
    # Harmonize if necessary
    send_message(sprintf("Harmonizing LSI SVD PCs for round %s...", i))
    harmonized_pcs <- harmony::HarmonyMatrix(
      data_mat  = LSIi$matSVD,
      meta_data = seurat_obj@meta.data,
      vars_use  = covariates, # Covariates to 'harmonize'
      do_pca    = FALSE
    )
    LSIi$matSVD <- harmonized_pcs
    
    # Create and store dimension reduction
    reducName <- paste0("LSI_iter", i)
    send_message(sprintf("Creating dimension reduction object for %s...", reducName))
    seurat_obj[[reducName]] <- CreateDimReducObject(embeddings = LSIi$matSVD, key = sprintf("LSI%s_", i), assay = "RNA")
    
    # Find neighbors and clusters
    send_message("Finding nearest neighbors...")
    seurat_obj <- Seurat::FindNeighbors(object = seurat_obj, reduction = reducName, dims = 1:nPCs, force.recalc = TRUE)
    
    send_message(sprintf("Clustering with resolution %s...", resolution[i]))
    seurat_obj <- Seurat::FindClusters(object = seurat_obj, resolution = resolution[i], algorithm = 4, method = "igraph")
    
    # Store cluster identities
    clusters <- Idents(seurat_obj)
    
    # Store information
    lsiOut[[reducName]] <- list(
      lsiMat = LSIi$matSVD,
      svd = LSIi$svd,
      varFeatures = varGenes, 
      clusters = clusters
    )
    
    send_message(sprintf("Completed LSI round %d", i))
  }
  
  # Return results
  return(list(
    seurat_obj = seurat_obj,
    lsiOut = lsiOut
  ))
}

#' Function to run LSI on a Seurat object
#'
#' @param mat raw count matrix from Seurat object.
#' @param nComponents number of right singular vectors to estimate.
#' @param binarize whether to binarize dgCMatrix class
#' @return list containing LSI results


LSI_function <- function(mat, nComponents, binarize = FALSE){
  
  if(binarize){
    message("Binarizing matrix...")
    # The 'x' slot of the dgCMatrix class contains the non-zero elements of the matrix
    mat@x[mat@x > 0] <- 1
  }
  
  #Calculate RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  
  # Calculate TF-IDF
  message(sprintf("Calculating TF-IDF with %s features (terms) and %s cells (documents)...", nrow(mat), ncol(mat)))
  start <- Sys.time()
  scaleTo <- 10^4
  tf <- t(t(mat) / colSm)
  idf <- as(ncol(mat) / rowSm, "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% tf
  # Log transform TF-IDF
  tfidf <- run_SparseLogX(tfidf, logtype="log", scale=TRUE, scaleFactor=scaleTo)
  
  # Clean up
  rm(tf)
  rm(idf)
  invisible(gc())
  
  # Calculate SVD for LSI
  message("Calculating SVD for LSI...")
  svd <- irlba::irlba(tfidf, nv=nComponents, nu=nComponents)
  svdDiag <- Matrix::diag(x=svd$d)
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC", seq_len(ncol(matSVD)))
  
  # Return matSVD and svd
  message(sprintf("LSI complete: %s minutes", round(difftime(Sys.time(), start, units="mins"), 3)))
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

