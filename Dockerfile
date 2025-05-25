# Dockerfile
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
    && rm -rf /var/lib/apt/lists/*

# Install R packages in order of dependencies
RUN R -e "install.packages(c('remotes', 'devtools', 'BiocManager'), repos='https://cran.rstudio.com/')"

# Install Bioconductor packages first
RUN R -e "BiocManager::install(c('edgeR', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'celda', 'GenomicFeatures', 'GenomicRanges', 'AnnotationDbi'), force = TRUE)"

# Install CRAN packages in batches to avoid timeouts
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyFiles', 'shinyjs', 'shinycssloaders'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages(c('DT', 'ggplot2', 'Matrix', 'dplyr', 'irlba'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages(c('cowplot', 'plyr', 'fs', 'R6', 'harmony'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages(c('sparseMatrixStats', 'matrixStats'), repos='https://cran.rstudio.com/')"

# Install Seurat and related packages
RUN R -e "install.packages(c('SeuratObject', 'Seurat'), repos='https://cran.rstudio.com/')"

# Install GitHub packages
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

# Set working directory
WORKDIR /app

# Copy package source
COPY . /app/

# Install the GUIdedRNA package
RUN R -e "devtools::install('.', dependencies = FALSE, upgrade = 'never')"

# Create directories and set permissions
RUN mkdir -p /data /output /srv/shiny-server
RUN groupadd -r shiny && useradd -r -g shiny shiny
RUN chown -R shiny:shiny /data /output /srv/shiny-server

# Copy the app to shiny server directory  
RUN cp -r /app/inst/shiny-app/* /srv/shiny-server/

# Expose port
EXPOSE 3838

# Create startup script
RUN echo '#!/bin/bash\nR -e "GUIdedRNA::launch_GUIdedRNA(port = 3838, host = \"0.0.0.0\", launch.browser = FALSE)"' > /start.sh && \
    chmod +x /start.sh

# Run as shiny user
USER shiny

# Start the application
CMD ["/start.sh"]
