# scRNA-seq Shiny App

This repository contains a Docker-based Shiny application that can be used for single-cell RNA analysis. It's primarily designed to use adjusted Latent Semantic Indexing for clusterization, an alternative to K-means clustering, with the intent to provide a higher annotation accuracy specifically in higher cell count datasets.

## Quick Start

### Prerequisites

- [Git](httpsgit-scm.comdownloads)
- [Docker](httpswww.docker.comproductsdocker-desktop)

### Running the App

On LinuxMac

```bash
git clone httpsgithub.comyourusernamescRNA-shiny-app.git
cd scRNA-shiny-app
chmod +x run.sh
.run.sh
```

On Windows

```cmd
git clone httpsgithub.comyourusernamescRNA-shiny-app.git
cd scRNA-shiny-app
setup.bat
```

Once running, open your browser and go to httplocalhost3838

## Future Extensions

This demo provides a foundation that can be extended to create a full single-cell RNA analysis pipeline with features like

- Data uploading and preprocessing
- Quality control and filtering
- Normalization and scaling
- Dimension reduction (PCA, t-SNE, UMAP)
- Clustering and marker gene identification
- Interactive visualization of results

## Customizing for scRNA-seq Analysis

To extend this app for real scRNA-seq analysis

1. Add Bioconductor packages to the Dockerfile
   ```
   RUN R -e install.packages('BiocManager'); BiocManagerinstall(c('Seurat', 'SingleCellExperiment', 'scater', 'scran'))
   ```

2. Modify the `app.R` file to include your analysis pipeline
3. Add example datasets in the `data` directory

## Repository Structure

```
scRNA-shiny-app
├── README.md          # This file
├── Dockerfile         # Docker configuration
├── run.sh             # Script to build and run on LinuxMac
├── app               # Shiny application files
│   ├── app.R          # Main Shiny app code
│   ├── data          # Example data files
│   └── www           # Static assets (CSS, images)
└── setup             # Setup scripts and utilities
    ├── setup.sh       # LinuxMac setup
    └── setup.bat      # Windows setup
```

## System Requirements

The Docker container is configured to use the resources allocated to Docker on your system. For analyzing larger scRNA-seq datasets, you may need to increase the resources allocated to Docker in Docker Desktop settings.

## Troubleshooting

- Port conflicts If port 3838 is already in use, modify the port number in `run.sh` or `setup.bat`
- Memory issues Increase memory allocation to Docker in Docker Desktop settings
- Container not starting Check Docker logs for more detailed error messages
