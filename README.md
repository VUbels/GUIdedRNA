# GUIdedRNA

A comprehensive Shiny application for guided RNA-seq analysis including quality control, preprocessing, dimensionality reduction, clustering, and cell type annotation.

## Features

- **Quality Control**: Interactive QC metrics visualization and filtering
- **Preprocessing**: Doublet removal and ambient RNA correction
- **LSI Analysis**: Two-round iterative Latent Semantic Indexing
- **Cell Type Annotation**: Guided clustering and annotation workflow
- **Integration**: Combine results from multiple analysis rounds
- **Export**: Comprehensive results export in multiple formats

## Installation

### From GitHub
```r
# Install devtools if you haven't already
install.packages("devtools")

# Install BiocManager dependencies

BiocManager::install(c("sparseMatrixStats", "AnnotationDbi", "edgeR", "GenomicRanges", "GenomicFeatures", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "celda", "decontX"))

# Install GUIdedRNA
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
devtools::install_github("VUbels/GUIdedRNA")
```

### Using Docker
```bash
# Clone the repository
git clone https://github.com/VUbels/GUIdedRNA.git
cd GUIdedRNA

# Quick deployment
chmod +x deploy.sh
./deploy.sh

# Then open http://localhost:3838
```

## Usage

### R Package
```r
library(GUIdedRNA)

# Launch the application
launch_GUIdedRNA()

# Launch with custom settings
launch_GUIdedRNA(port = 8080, host = "localhost")
```

### Docker
```bash
# Using docker-compose (recommended)
docker-compose up

# Or manual Docker
docker build -t guidedrna .
docker run -p 3838:3838 -v ./data:/data -v ./output:/output guidedrna
```

## Workflow

1. **Setup**: Upload 10X Genomics data or folder containing multiple datasets
2. **Sample Information**: Add metadata and sample attributes
3. **Quality Control**: Set filtering parameters and visualize QC metrics
4. **Preprocessing**: Run doublet removal and ambient RNA correction
5. **LSI Round 1**: Initial dimensionality reduction and broad clustering
6. **Initial Clustering**: Assign broad cell types using marker genes
7. **LSI Round 2**: Refined analysis on cell type subsets
8. **Final Clustering**: Detailed cell type annotation
9. **Integration**: Combine all results into final annotations
10. **Download**: Export processed data and results

## Requirements

- R >= 4.0.0
- Seurat >= 5.0.0
- Docker (for containerized deployment)

## Support

For issues and questions, please visit our [GitHub Issues](https://github.com/VUbels/GUIdedRNA/issues) page.
