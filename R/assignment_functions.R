#' 
#'
#' Calculate the Term Frequency - Inverse Document Frequency (TF-IDF) for a feature x cell counts
#' matrix, then calculate the Singular Value Decomposition of that matrix, which is then used as
#' input for Seurat's SNN clustering
#' 
#' @param mat raw count matrix from Seurat object in sparse format.
#' @param nComponents number of right singular vectors to estimate.
#' @param binarize whether to binarize matrix, standardized to FALSE.

# Identify genes we want to blacklist during clustering

mt.genes <- grep(pattern = "^MT-", x = rownames(rawCounts), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(rawCounts), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(rawCounts), value = TRUE)

# X/Y chromosome genes:
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
geneGR <- GenomicFeatures::genes(txdb)
sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
matchedGeneSymbols <- select(org.Hs.eg.db,
                             keys = sexGenesGR$gene_id,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "ENTREZID")
sexChr.genes <- matchedGeneSymbols$SYMBOL


# Genes to ignore (just for clustering purposes)
blacklist.genes <- c(
  mt.genes,
  sexChr.genes,
  s.genes,
  g2m.genes,
  RPS.genes,
  RPL.genes
)