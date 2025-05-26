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

# Copy the entire package (your structure is already correct!)
COPY . /tmp/GUIdedRNA/

# Verify the structure
RUN echo "=== Package structure ===" && \
    ls -la /tmp/GUIdedRNA/ && \
    echo "=== Shiny app files ===" && \
    ls -la /tmp/GUIdedRNA/inst/shiny-app/ && \
    echo "=== R package files ===" && \
    ls -la /tmp/GUIdedRNA/R/

# Build and install the package
WORKDIR /tmp/GUIdedRNA
RUN R CMD build . --no-build-vignettes
RUN R CMD INSTALL --build *.tar.gz

# Verify installation
RUN R -e "library(GUIdedRNA); print(system.file('shiny-app', package = 'GUIdedRNA')); list.files(system.file('shiny-app', package = 'GUIdedRNA'))"

# Create data and output directories
RUN mkdir -p /data /output
RUN chmod 755 /data /output

# Create startup script
RUN echo '#!/bin/bash\n\
echo "Starting GUIdedRNA application..."\n\
echo "Checking package installation..."\n\
R --slave -e "library(GUIdedRNA); cat(\"Package loaded successfully\\n\")"\n\
echo "Starting Shiny app..."\n\
R --slave -e "GUIdedRNA::launch_GUIdedRNA(port = 3838, host = \"0.0.0.0\", launch.browser = FALSE)"' > /start.sh

RUN chmod +x /start.sh

# Clean up
RUN rm -rf /tmp/GUIdedRNA

# Set working directory
WORKDIR /

# Expose port
EXPOSE 3838

# Start the application
CMD ["/start.sh"]