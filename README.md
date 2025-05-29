# GUIdedRNA

A comprehensive Shiny application for guided RNA-seq analysis including quality control, preprocessing, dimensionality reduction, clustering, and cell type annotation.

## Features

- **Quality Control**: Interactive QC metrics visualization and filtering
- **Preprocessing**: Doublet removal and ambient RNA correction
- **LSI Analysis**: Two-round iterative Latent Semantic Indexing
- **Cell Type Annotation**: Guided clustering and annotation workflow
- **Integration**: Combine results from multiple analysis rounds
- **Export**: Comprehensive results export in multiple formats

## Installation and use for regular use

**For Regular Users (Recommended)**
No technical knowledge required - just double-click and wait!  

### Windows/IOS/Linux Users
1) Download this repository (green "Code" button → "Download ZIP")  
2) Extract to your Desktop  
3) Download and install https://www.docker.com/products/docker-desktop/  
4) Run Docker Desktop  
5a) Double-click INSTALL_AND_RUN.bat (WINDOWS)
5b) Double-click INSTALL_AND_RUN.sh (IOS/Linux)
6) Let Docker build container (5-15 minutes)
7) If browser does not open simply go to http://localhost:3838

After building docker container is completed once, subsequently running app is instantaneous.

### For Advanced Users & Developers
WSL2 Users (Windows Subsystem for Linux)  
If you need Windows drive access and advanced volume mounting:  
If working through command line simply open in any browser at http://localhost:3838  

### Installation of GUIdedRNA library and use through IDE

```r
# Install devtools
install.packages("devtools")

# Install BiocManager dependencies
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("sparseMatrixStats", "AnnotationDbi", "edgeR", "GenomicRanges", "GenomicFeatures", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "celda", "decontX"))

# Install GUIdedRNA
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
devtools::install_github("VUbels/GUIdedRNA")
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
