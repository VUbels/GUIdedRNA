# Use Bioconductor base image with proper R/Bioconductor version alignment
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgsl-dev \
    libglpk-dev \
    libhdf5-dev \
    libfftw3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install core Shiny packages using BiocManager for version consistency
RUN R -e "BiocManager::install(c('shiny', 'shinydashboard', 'shinyFiles', 'shinyjs', 'shinycssloaders', 'DT'), ask=FALSE)"

# Install required Bioconductor packages (all at once for version consistency)
RUN R -e "BiocManager::install(c('sparseMatrixStats', 'AnnotationDbi', 'edgeR', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'celda', 'decontX'), ask=FALSE)"

# Install CRAN packages through BiocManager for version consistency
RUN R -e "BiocManager::install(c('Matrix', 'irlba', 'Rcpp', 'RcppArmadillo', 'fields', 'KernSmooth', 'ROCR', 'parallel', 'plyr', 'matrixStats', 'dplyr', 'ggplot2', 'cowplot', 'fs', 'R6', 'remotes', 'devtools'), ask=FALSE)"

# Install Seurat ecosystem through BiocManager
RUN R -e "BiocManager::install(c('SeuratObject', 'Seurat', 'harmony'), ask=FALSE)"

# Install DoubletFinder with error handling (optional dependency)
RUN R -e "tryCatch({ remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade='never'); cat('DoubletFinder installed successfully\\n') }, error=function(e) { cat('DoubletFinder failed, will be installed by GUIdedRNA if needed\\n') })"
RUN R -e "tryCatch({ remotes::install_github('immunogenomics/presto', upgrade='never'); cat('Presto installed successfully\\n') }, error=function(e) { cat('Presto failed, finding marker genes will use standard Wilcoxen\\n') })"

# Install GUIdedRNA - let it handle its own dependencies
RUN R -e "remotes::install_github('VUbels/GUIdedRNA', dependencies=TRUE, upgrade='never')"

# Verify installation
RUN R -e "library(GUIdedRNA); cat('GUIdedRNA loaded successfully\\n')"

# Create data and output directories
RUN mkdir -p /data /output
WORKDIR /app

# Create simple launch script
RUN echo 'library(GUIdedRNA)\n\
\n\
# Launch GUIdedRNA application\n\
if(exists("launch_GUIdedRNA")) {\n\
  cat("Using launch_GUIdedRNA function\\n")\n\
  launch_GUIdedRNA(host="0.0.0.0", port=3838, launch.browser=FALSE)\n\
} else {\n\
  cat("Using fallback launch method\\n")\n\
  app_dir <- system.file("shiny-app", package="GUIdedRNA")\n\
  if(app_dir != "" && dir.exists(app_dir)) {\n\
    shiny::runApp(appDir=app_dir, host="0.0.0.0", port=3838, launch.browser=FALSE)\n\
  } else {\n\
    stop("GUIdedRNA app directory not found")\n\
  }\n\
}' > /app/launch.R

# Expose port
EXPOSE 3838

# Simple launch command
CMD ["R", "--no-restore", "--file=/app/launch.R"]
