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

# Install system-level R packages first
RUN R -e "install.packages(c('Matrix', 'irlba', 'Rcpp', 'RcppArmadillo', 'fields', 'KernSmooth', 'ROCR', 'parallel', 'plyr', 'matrixStats'), repos='https://cran.rstudio.com/')"

# Install core Shiny packages
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyFiles', 'shinyjs', 'shinycssloaders', 'DT'), repos='https://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('sparseMatrixStats', 'AnnotationDbi', 'edgeR', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'celda', 'decontX'))"

# Install data manipulation packages
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'cowplot', 'fs', 'R6'), repos='https://cran.rstudio.com/')"

# Install Seurat ecosystem (critical to do this after Bioconductor)
RUN R -e "install.packages(c('Seurat'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages(c('SeuratObject'), repos='https://cran.rstudio.com/')"

# Install harmony
RUN R -e "install.packages('harmony', repos='https://cran.rstudio.com/')"

# Install GitHub packages
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"
RUN R -e "remotes::install_github('VUbels/GUIdedRNA')"

# Verify installation
RUN R -e "library(GUIdedRNA); cat('GUIdedRNA installed successfully from GitHub\n')"

# Create data and output directories
RUN mkdir -p /data /output /app

# Set working directory
WORKDIR /app

# Create startup script as R file
RUN echo 'library(shiny)\n\
library(GUIdedRNA)\n\
\n\
cat("Loading GUIdedRNA package...\\n")\n\
\n\
# Check if launch function exists\n\
if(exists("launch_GUIdedRNA")) {\n\
    cat("Using launch_GUIdedRNA function\\n")\n\
    launch_GUIdedRNA(host = "0.0.0.0", port = 3838, launch.browser = FALSE)\n\
} else {\n\
    cat("Looking for app directory...\\n")\n\
    app_dir <- system.file("shiny-app", package = "GUIdedRNA")\n\
    cat("App directory:", app_dir, "\\n")\n\
    if(app_dir != "" && dir.exists(app_dir)) {\n\
        cat("Starting Shiny app from directory\\n")\n\
        shiny::runApp(appDir = app_dir, host = "0.0.0.0", port = 3838, launch.browser = FALSE)\n\
    } else {\n\
        stop("GUIdedRNA app directory not found")\n\
    }\n\
}' > /app/start.R

# Expose port
EXPOSE 3838

# Run the application
CMD ["R", "--no-restore", "--file=/app/start.R"]