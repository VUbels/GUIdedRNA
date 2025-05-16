#' Variable gene function selection enabling.
#'
#' Work around for staying with sparse matrix format in seurat object to reduce object size.
#' 
#' @param mat raw count matrix from Seurat object.
#' @param nvar number of variable genes to be found.
#' @param blacklist list of genes to exclude from variable gene analysis.

getVarGenes <- function(mat, nvar = 2000, blacklist = NULL){
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

sparseLogX <- function(spmat, logtype="log2", scale=FALSE, scaleFactor=10^4){
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

#' TF-IDF LSI adapted from Jeff Granja, who adapted from flyATAC (i.e. Cusanovich et al. 2018)
#'
#' Calculate the Term Frequency - Inverse Document Frequency (TF-IDF) for a feature x cell counts
#' matrix, then calculate the Singular Value Decomposition of that matrix, which is then used as
#' input for Seurat's SNN clustering
#' 
#' @param mat raw count matrix from Seurat object in sparse format.
#' @param nComponents number of right singular vectors to estimate.
#' @param binarize whether to binarize matrix, standardized to FALSE.


runLSI <- function(mat, nComponents, binarize = FALSE){
  
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
  tfidf <- sparseLogX(tfidf, logtype="log", scale=TRUE, scaleFactor=scaleTo)
  
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