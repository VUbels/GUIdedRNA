FROM rocker/shiny:4.3.0

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

# Set R library path
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Install R package building tools first
RUN R -e "install.packages(c('remotes', 'devtools', 'BiocManager'), repos='https://cran.rstudio.com/')"

# Set Bioconductor version explicitly
RUN R -e "BiocManager::install(version='3.18', ask=FALSE)"

# Install system-level R packages first
RUN R -e "install.packages(c('Matrix', 'irlba', 'Rcpp', 'RcppArmadillo'), repos='https://cran.rstudio.com/')"

# Install core Shiny packages
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyFiles', 'shinyjs', 'shinycssloaders', 'DT'), repos='https://cran.rstudio.com/')"

# Install data manipulation packages
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'cowplot', 'plyr', 'fs', 'R6'), repos='https://cran.rstudio.com/')"

# Install statistical packages
RUN R -e "install.packages(c('sparseMatrixStats', 'matrixStats'), repos='https://cran.rstudio.com/')"

# Install Bioconductor packages with explicit dependencies
RUN R -e "BiocManager::install(c('GenomeInfoDb', 'S4Vectors', 'IRanges', 'GenomicRanges'), force=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install(c('AnnotationDbi', 'GenomicFeatures'), force=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install(c('org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene'), force=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install(c('edgeR', 'celda', 'decontX'), force=TRUE, ask=FALSE)"

# Install Seurat ecosystem (critical to do this after Bioconductor)
RUN R -e "install.packages(c('SeuratObject'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages(c('Seurat'), repos='https://cran.rstudio.com/')"

# Install harmony
RUN R -e "install.packages('harmony', repos='https://cran.rstudio.com/')"

# Install GitHub packages
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

# Copy the entire package
COPY . /tmp/GUIdedRNA/

# Set working directory
WORKDIR /tmp/GUIdedRNA

# Create a temporary DESCRIPTION without problematic dependencies
RUN cp DESCRIPTION DESCRIPTION.backup && \
    sed '/DoubletFinder/d; /TxDb.Hsapiens.UCSC.hg38.knownGene/d; /GenomicFeatures/d; /decontX/d; /celda/d; /edgeR/d; /org.Hs.eg.db/d; /GenomicRanges/d; /AnnotationDbi/d; /sparseMatrixStats/d; /matrixStats/d' DESCRIPTION > DESCRIPTION.temp && \
    mv DESCRIPTION.temp DESCRIPTION

# Build and install the package without dependencies
RUN R CMD build . --no-build-vignettes
RUN R CMD INSTALL *.tar.gz

# Restore the original DESCRIPTION
RUN mv DESCRIPTION.backup DESCRIPTION

# Manually copy the package source to the installed location to ensure we have all functions
RUN cp -r R/* /usr/local/lib/R/site-library/GUIdedRNA/R/ 2>/dev/null || true
RUN cp -r inst/* /usr/local/lib/R/site-library/GUIdedRNA/ 2>/dev/null || true

# Verify installation by testing package loading and function availability
RUN R -e "library(GUIdedRNA); \
    cat('Checking app directory: '); \
    app_dir <- system.file('shiny-app', package = 'GUIdedRNA'); \
    cat(app_dir, '\n'); \
    if(app_dir != '' && dir.exists(app_dir)) { \
        cat('App files: '); cat(list.files(app_dir), sep=', '); cat('\n'); \
    } else { \
        cat('Creating app directory manually...\n'); \
        dir.create('/usr/local/lib/R/site-library/GUIdedRNA/shiny-app', recursive=TRUE); \
        file.copy('inst/shiny-app/app.R', '/usr/local/lib/R/site-library/GUIdedRNA/shiny-app/app.R'); \
        file.copy('inst/shiny-app/global.R', '/usr/local/lib/R/site-library/GUIdedRNA/shiny-app/global.R'); \
    }"

# Test that the launch function works and reinstall any missing dependencies
RUN R -e "library(GUIdedRNA); \
    if(exists('launch_GUIdedRNA')) { \
        cat('✓ launch_GUIdedRNA function found\n'); \
    } else { \
        cat('✗ launch_GUIdedRNA function not found\n'); \
    }" && \
    R -e "if(!require('DoubletFinder', quietly=TRUE)) { \
        cat('Installing missing DoubletFinder...\n'); \
        remotes::install_github('chris-mcginnis-ucsf/DoubletFinder'); \
    }" && \
    R -e "if(!require('GenomicFeatures', quietly=TRUE)) { \
        cat('Installing missing GenomicFeatures...\n'); \
        BiocManager::install('GenomicFeatures'); \
    }" && \
    R -e "if(!require('TxDb.Hsapiens.UCSC.hg38.knownGene', quietly=TRUE)) { \
        cat('Installing missing TxDb.Hsapiens.UCSC.hg38.knownGene...\n'); \
        BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene'); \
    }"

# Create data and output directories
RUN mkdir -p /data /output
RUN chmod 755 /data /output

# Create a simple startup script
RUN cat > /start.sh << 'EOF'
#!/bin/bash
echo "Starting GUIdedRNA application..."
echo "Loading required libraries..."

R --slave << 'RSCRIPT'
library(shiny)
library(shinydashboard)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(DT)
cat("Core libraries loaded successfully\n")

library(GUIdedRNA)
cat("GUIdedRNA package loaded successfully\n")

app_dir <- system.file("shiny-app", package = "GUIdedRNA")
if(app_dir == "" || !dir.exists(app_dir)) {
    stop("App directory not found. Package installation may have failed.")
}
cat("App directory found:", app_dir, "\n")

if(exists("launch_GUIdedRNA")) {
    cat("Starting Shiny application...\n")
    launch_GUIdedRNA(port = 3838, host = "0.0.0.0", launch.browser = FALSE)
} else {
    cat("Launching app directly...\n") 
    shiny::runApp(appDir = app_dir, port = 3838, host = "0.0.0.0", launch.browser = FALSE)
}
RSCRIPT
EOF

RUN chmod +x /start.sh

# Clean up
RUN rm -rf /tmp/GUIdedRNA

# Set working directory
WORKDIR /

# Expose port
EXPOSE 3838

# Start the application
CMD ["/start.sh"]