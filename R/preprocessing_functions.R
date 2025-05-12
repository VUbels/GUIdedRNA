#' Assigns expected doublet ground variable
#'
#' Function returns the expected doublet rate based on 10X data dependent on cell count rather than flat rate. Increasing cell count in smaller samples.
#' Assign_ExpectedDoublet_variable()


assign_ExpectedDoublet <- function(val) {
  
  # Initialize expected_doublet value
  expected_doublet<- NA
  
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
#' Fully runs the doublet removal protocol through DoubletFinder() based on variable from Assign_ExpectedDoublet_variable(), Returns Seurat object with TRUE/FALSE assignment to duplicate cells in META data


preprocess_DoubletRemoval <- function(seurat_obj)
{
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)
  
  sweep.res.list_obj <- paramSweep(seu_obj, PCs = 1:10, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_obj)
  
  val <- length(seurat_obj@assays[["RNA"]]@counts@p)
  expected_doublet <- assign_ExpectedDoublet(val)

  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(expected_doublet*nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seurat_obj <- doubletFinder(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
  
  # Create the plot
  
  range_x <- range(bcmvn$pK)
  positions <- seq(range_x[1], range_x[2], length.out = 6)
  
  plot_obj <- ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    ggtitle(seurat_names) +
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
  
  output_filename <- paste0(output_directory, "Clustering_01/DoubletFinderFunction_", seurat_names, ".png")
  ggsave(filename = output_filename, plot = plot_obj, width = 6, height = 4, bg = 'white')
  
  paste0(output_directory, output_filename)
  
  return(seurat_obj)
  
}

#####################################
# Removal of decontaminated RNA
#####################################

runDecontX <- function(seurat_obj, seed=1){
  counts <- GetAssayData(object = seurat_obj, slot = "counts")
  clusters <- Idents(seurat_obj) %>% as.numeric()
  
  # Run on only expressed genes
  x <- counts[rowSums(counts)>0,]
  message(sprintf("Running decontX on %s cells with %s non-zero genes...", dim(x)[2], dim(x)[1]))
  decon <- decontX(x, z=clusters, verbose=TRUE, seed=seed)
  
  # Save desired information back to Seurat Object
  # We will place the estimated 'decontaminated counts' in place of the original counts ('RNA')
  # and keep the original counts as a separate assay called 'origCounts'
  newCounts <- decon$decontXcounts
  # Add back unexpressed genes and sort according to original counts
  newCounts <- rbind(newCounts, counts[rowSums(counts)==0,])[rownames(counts),]
  seurat_obj[["RNA"]]@counts <- as(round(newCounts), "sparseMatrix")
  seurat_obj$estConp <- decon$contamination # Estimated 'contamination proportion, 0 to 1'
  
  return(seurat_obj)
}

Remove_lowRNA <- function(seurat_obj){
  nPreDecon <- dim(seurat_obj)[2]
  seurat_obj <- subset(seurat_obj, subset = (nFeature_RNA > minFeatures & nCount_RNA > minCounts))
  nPostDecon <- dim(seurat_obj)[2]
  
  message(sprintf("Removed %s cells with too little RNA after decontamination. %s cells remaining.", 
                  nPreDecon - nPostDecon, nPostDecon))
  
  return(seurat_obj)
  
}

#####################################
# Cluster scRNA using iterative LSI
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(Matrix)
})