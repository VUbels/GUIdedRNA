# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)

library(Seurat)
library(ggplot2)
library(DT)
library(Matrix)
library(dplyr)

#library(GUIdedRNA)


options(shiny.maxRequestSize = 100 * 1024^3)

# Define UI
ui <- dashboardPage(
  
  skin = c("black"),
  dashboardHeader(title = "GUIdedRNA - scRNA analysis interface"),
  
  # Sidebar with analysis steps
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Setup", tabName = "upload", icon = icon("upload"), selected = TRUE),
      menuItem("Quality Control", tabName = "qc", icon = icon("check-circle")),
      menuItem("Preprocessing", tabName = "preprocess", icon = icon("filter")),
      menuItem("Normalization", tabName = "normalize", icon = icon("chart-line")),
      menuItem("Feature Selection", tabName = "features", icon = icon("list")),
      menuItem("Dimensionality Reduction", tabName = "dimreduce", icon = icon("project-diagram")),
      menuItem("Clustering", tabName = "cluster", icon = icon("object-group")),
      menuItem("Cell Type Annotation", tabName = "annotate", icon = icon("tags")),
      menuItem("Download Results", tabName = "download", icon = icon("download"))
    )
  ),
  
  # Main panel with tabbed content
  dashboardBody(
    shinyjs::useShinyjs(),
    
    tabItems(
      # Upload Data tab
      tabItem(tabName = "upload",
              fluidRow(
                box(
                  title = "Upload 10X Genomics files",
                  width = 12,
                  fileInput("matrixFile", "Matrix (.mtx)"),
                  fileInput("featuresFile", "Features/Genes (.tsv/.txt)"),
                  fileInput("barcodesFile", "Barcodes (.tsv/.txt)"),
                  actionButton("loadData_file", "Load Data")
                ),
              ),
              
              fluidRow(
                box(
                  title = "Upload Folder or list of Folders",
                  width = 12,
                  shinyDirButton("folderBtn_single", "Select Directory", "Choose a directory with 10X data"),
                  helpText("Select directory containing files, or directory containing dataset folders"),
                  textOutput("selectedFolder"),
                  actionButton("loadData_folder", "Load Data")
                ),
              ),
              
              fluidRow(
                box(
                  title = "Set Output Directory",
                  width = 12,
                  shinyDirButton("folderBtn_output", "Select Directory", "Choose an output directory"),
                  helpText("Select directory that will contain all output"),
                  textOutput("selectedFolder"),
                  actionButton("outputData_folder", "Set Output Folder")
                ),
              ),
              
              
              fluidRow(
                box(
                  title = "Data Summary",
                  width = 12,
                  verbatimTextOutput("dataSummary")
                )
              )
      ),
      
      # Quality Control tab
      tabItem(tabName = "qc",
              fluidRow(
                box(
                  title = "QC Parameters",
                  width = 4,
                  
                  # Side by side gene count inputs
                  fluidRow(
                    column(6, numericInput("minGenes", "Min Genes", min = 0, value = 200, step = 100)),
                    column(6, numericInput("maxGenes", "Max Genes", min = 0, value = 5000, step = 100))
                  ),
                  
                  # Side by side feature inputs
                  fluidRow(
                    column(6, numericInput("minFeatures", "Min Features", min = 0, value = 100, step = 100)),
                    column(6, numericInput("maxFeatures", "Max Features", min = 0, value = 3000, step = 100))
                  ),
                  
                  # Single mitochondrial percentage input
                  numericInput("maxMito", "Max Mitochondrial %", min = 0, max = 100, value = 20),
                  
                  actionButton("runQC", "Run QC")
                ),
                box(
                  title = "QC Metrics",
                  width = 8,
                  plotOutput("qcPlot") %>% withSpinner()
                )
              ),
              fluidRow(
                box(
                  title = "Filtered Data Summary",
                  width = 12,
                  verbatimTextOutput("qcSummary")
                )
              )
      ),
      
      
      tabItem(tabName = "preprocess",
              fluidRow(
                box(
                  title = "Preprocessing Options",
                  width = 4,
                  checkboxGroupInput("preprocessMethods", "Select Preprocessing Methods", 
                                     choices = c("Doublet Removal" = "doublet", 
                                                 "Ambient RNA Removal" = "ambient"),
                                     selected = c("doublet", "ambient")
                                     
                  ),
                  actionButton("commitPreprocessing", "Run Selected Methods", 
                               class = "btn-success", 
                               style = "width: 100%;")
                ),
                box(
                  title = "Preprocessing Results",
                  width = 8,
                  plotOutput("preprocessingPlot") %>% withSpinner()
                )
              )
      ),
      
      # Normalization tab
      tabItem(tabName = "normalize",
              fluidRow(
                box(
                  title = "Normalization Method",
                  width = 4,
                  selectInput("normMethod", "Method", 
                              choices = c("LogNormalize", "SCTransform"), 
                              selected = "LogNormalize"),
                  conditionalPanel(
                    condition = "input.normMethod == 'LogNormalize'",
                    numericInput("scaleFactor", "Scale Factor", value = 10000)
                  ),
                  actionButton("runNorm", "Run Normalization")
                ),
                box(
                  title = "Normalization Results",
                  width = 8,
                  plotOutput("normPlot") %>% withSpinner()
                )
              )
      ),
      
      # Feature Selection tab
      tabItem(tabName = "features",
              fluidRow(
                box(
                  title = "Variable Features Parameters",
                  width = 4,
                  numericInput("nFeatures", "Number of Features", value = 2000),
                  actionButton("findVarFeatures", "Find Variable Features")
                ),
                box(
                  title = "Variable Features Plot",
                  width = 8,
                  plotOutput("varFeaturesPlot") %>% withSpinner()
                )
              ),
              fluidRow(
                box(
                  title = "Top Variable Features",
                  width = 12,
                  DTOutput("topFeaturesTable")
                )
              )
      ),
      
      # Dimensionality Reduction tab
      tabItem(tabName = "dimreduce",
              fluidRow(
                box(
                  title = "PCA Parameters",
                  width = 4,
                  numericInput("nPCs", "Number of PCs", value = 30),
                  actionButton("runPCA", "Run PCA")
                ),
                box(
                  title = "PCA Results",
                  width = 8,
                  plotOutput("pcaPlot") %>% withSpinner(),
                  plotOutput("elbowPlot") %>% withSpinner()
                )
              ),
              fluidRow(
                box(
                  title = "UMAP/t-SNE Parameters",
                  width = 4,
                  numericInput("pcDims", "PCs to Use", value = 20),
                  selectInput("reduction", "Reduction Method", 
                              choices = c("UMAP", "t-SNE"), 
                              selected = "UMAP"),
                  actionButton("runReduction", "Run Reduction")
                ),
                box(
                  title = "Reduction Results",
                  width = 8,
                  plotOutput("reductionPlot") %>% withSpinner()
                )
              )
      ),
      
      # Clustering tab
      tabItem(tabName = "cluster",
              fluidRow(
                box(
                  title = "Clustering Parameters",
                  width = 4,
                  selectInput("clusterAlgo", "Algorithm", 
                              choices = c("Louvain", "Leiden"), 
                              selected = "Louvain"),
                  sliderInput("resolution", "Resolution", min = 0.1, max = 2.0, value = 0.8, step = 0.1),
                  actionButton("runClustering", "Run Clustering")
                ),
                box(
                  title = "Clustering Results",
                  width = 8,
                  plotOutput("clusterPlot") %>% withSpinner()
                )
              ),
              fluidRow(
                box(
                  title = "Cluster Statistics",
                  width = 12,
                  DTOutput("clusterStats")
                )
              )
      ),
      
      # Cell Type Annotation tab
      tabItem(tabName = "annotate",
              fluidRow(
                box(
                  title = "Marker Gene Analysis",
                  width = 4,
                  numericInput("maxMarkers", "Max Markers per Cluster", value = 10),
                  actionButton("findMarkers", "Find Markers")
                ),
                box(
                  title = "Marker Gene Results",
                  width = 8,
                  DTOutput("markerTable") %>% withSpinner()
                )
              ),
              fluidRow(
                box(
                  title = "Cell Type Assignment",
                  width = 12,
                  uiOutput("cellTypeUI"),
                  actionButton("assignCellTypes", "Assign Cell Types")
                )
              ),
              fluidRow(
                box(
                  title = "Annotated Clusters",
                  width = 12,
                  plotOutput("annotatedPlot") %>% withSpinner()
                )
              )
      ),
      
      # Download Results tab
      tabItem(tabName = "download",
              fluidRow(
                box(
                  title = "Download Options",
                  width = 12,
                  checkboxGroupInput("downloadItems", "Select Items to Download:",
                                     choices = c("Processed Seurat Object (RDS)" = "rds",
                                                 "Cell Metadata (CSV)" = "meta",
                                                 "Marker Genes (CSV)" = "markers",
                                                 "Dimensionality Reduction Coordinates (CSV)" = "dimred",
                                                 "Analysis Report (HTML)" = "report")),
                  downloadButton("downloadData", "Download Selected Items")
                )
              )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Reactive values to store data and analysis state
  values <- reactiveValues(
    seurat = NULL,
    original_seurat = NULL,
    qc_done = FALSE,
    preprocess_done = FALSE,
    norm_done = FALSE,
    features_done = FALSE,
    pca_done = FALSE,
    reduction_done = FALSE,
    cluster_done = FALSE,
    markers_done = FALSE,
    annotation_done = FALSE
  )

  # Define volumes in a platform-agnostic way
  volumes <- c()
  
  # Add home directory (works on all platforms)
  volumes <- c(volumes, Home = fs::path_home())
  
  # Add root for Linux/Mac
  if (.Platform$OS.type == "unix") {
    volumes <- c(volumes, "Root" = "/")
    
    # Add common WSL paths if running in WSL
    wsl_check <- system("grep -q Microsoft /proc/version", ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (wsl_check == 0) {  # Running in WSL
      # Add Windows drives that might be mounted
      if (dir.exists("/mnt/c")) volumes <- c(volumes, "C:" = "/mnt/c")
      if (dir.exists("/mnt/d")) volumes <- c(volumes, "D:" = "/mnt/d")
      if (dir.exists("/mnt/e")) volumes <- c(volumes, "E:" = "/mnt/e")
      if (dir.exists("/mnt/f")) volumes <- c(volumes, "F:" = "/mnt/f")
      if (dir.exists("/mnt/g")) volumes <- c(volumes, "G:" = "/mnt/g")
      
      # Add other common WSL locations
      if (dir.exists("/mnt")) volumes <- c(volumes, "Windows Drives" = "/mnt")
    }
    
    # Add other common Unix locations
    if (dir.exists("/media")) volumes <- c(volumes, "Media" = "/media")
    if (dir.exists("/mnt")) volumes <- c(volumes, "Mnt" = "/mnt")
  }
  
  # Add Windows drives if on Windows
  if (.Platform$OS.type == "windows") {
    for (letter in LETTERS) {
      drive <- paste0(letter, ":")
      if (dir.exists(drive)) {
        volumes <- c(volumes, drive)
        names(volumes)[length(volumes)] <- paste0(letter, ":")
      }
    }
  }
  
  # Initialize directory chooser
  shinyDirChoose(input, "folderBtn_single", roots = volumes, session = session)
  
  # Display selected folder
  output$selectedFolder <- renderText({
    if (!is.null(input$folderBtn_single)) {
      parseDirPath(volumes, input$folderBtn_single)
    } else {
      "No folder selected"
    }
  })
    
  # Data loading logic single dataset
  observeEvent(input$loadData_file, {
    req(input$matrixFile, input$featuresFile, input$barcodesFile)
    
    withProgress(message = 'Loading data...', {
      # Create temporary directory
      temp_dir <- tempdir()
      matrix_path <- file.path(temp_dir, "matrix.mtx")
      features_path <- file.path(temp_dir, "features.tsv")
      barcodes_path <- file.path(temp_dir, "barcodes.tsv")
      
      # Save uploaded files
      file.copy(input$matrixFile$datapath, matrix_path)
      file.copy(input$featuresFile$datapath, features_path)
      file.copy(input$barcodesFile$datapath, barcodes_path)
      
      # Read data
      counts <- Seurat::ReadMtx(
        mtx = matrix_path,
        features = features_path,
        cells = barcodes_path
      )
      
      # Create Seurat object
      seurat_obj <- Seurat::CreateSeuratObject(counts = counts)
      
      # Calculate percent mitochondrial
      seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      
      # Store in reactive values
      values$seurat <- seurat_obj
      values$original_seurat <- seurat_obj
    })
    
    # Show data summary
    output$dataSummary <- renderPrint({
      req(values$seurat)
      cat("Dataset loaded successfully.\n")
      cat(paste("Number of cells:", ncol(values$seurat), "\n"))
      cat(paste("Number of genes:", nrow(values$seurat), "\n"))
      cat("\nSample of gene names:\n")
      print(head(rownames(values$seurat)))
      cat("\nSample of cell barcodes:\n")
      print(head(colnames(values$seurat)))
    })
    
    # Navigate to QC tab
    updateTabItems(session, "tabs", "qc")
  })
  
  # Data loading for folder containing multiple datasets
  shinyDirChoose(input, "folderBtn_single", roots = volumes, session = session)
  
  # Add display for selected folder
  output$selectedFolder <- renderText({
    if (!is.null(input$folderBtn_single)) {
      parseDirPath(volumes, input$folderBtn_single)
    } else {
      "No folder selected"
    }
  })
  

  # Data loading for folder containing multiple datasets
  observeEvent(input$loadData_folder, {
    # Check if folder is selected
    if (is.null(input$folderBtn_single)) {
      showNotification("Please select a folder first", type = "error")
      return(NULL)
    }
    
    folder_path <- parseDirPath(volumes, input$folderBtn_single)
    
    # Check if path is valid
    if (is.null(folder_path) || folder_path == "") {
      showNotification("Please select a valid folder", type = "error")
      return(NULL)
    }
    
    withProgress(message = 'Processing folder...', {
      # Initialize a list to store Seurat objects
      seurat_list <- list()
      
      # 1. Process files in the main directory with pattern: name_matrix.mtx.gz, etc.
      main_files <- list.files(folder_path, pattern = "(_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz)$", 
                               full.names = TRUE, recursive = FALSE)
      
      # Find base names for files in the main directory
      main_base_names <- unique(sapply(strsplit(basename(main_files), 
                                                "_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz"), 
                                       `[`, 1))
      
      # Process datasets in the main directory
      for (base in main_base_names) {
        setProgress(detail = paste("Processing", base))
        
        mtx_file <- file.path(folder_path, paste0(base, "_matrix.mtx.gz"))
        features_file <- file.path(folder_path, paste0(base, "_features.tsv.gz"))
        cells_file <- file.path(folder_path, paste0(base, "_barcodes.tsv.gz"))
        
        # Check if all files exist
        if (file.exists(mtx_file) && file.exists(features_file) && file.exists(cells_file)) {
          # Read the matrix
          tryCatch({
            # Read the matrix
            counts <- Seurat::ReadMtx(
              mtx = mtx_file,
              features = features_file,
              cells = cells_file
            )
            
            # Create Seurat object
            seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = base)
            
            # Calculate percent mitochondrial
            seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
            
            # Add to list
            seurat_list[[base]] <- seurat_obj
          }, error = function(e) {
            showNotification(paste("Error processing", base, ":", e$message), type = "warning")
          })
        }
      }
      
      # 2. Process subdirectories with standard 10X naming
      subdirs <- list.dirs(folder_path, recursive = FALSE, full.names = TRUE)
      
      for (subdir in subdirs) {
        subdir_name <- basename(subdir)
        setProgress(detail = paste("Processing folder", subdir_name))
        
        # Construct the file paths for matrix, features, and barcodes
        # Check for both .gz and non-gz versions
        mtx_file <- NULL
        for (pattern in c("matrix.mtx.gz", "matrix.mtx")) {
          path <- file.path(subdir, pattern)
          if (file.exists(path)) {
            mtx_file <- path
            break
          }
        }
        
        features_file <- NULL
        for (pattern in c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")) {
          path <- file.path(subdir, pattern)
          if (file.exists(path)) {
            features_file <- path
            break
          }
        }
        
        cells_file <- NULL
        for (pattern in c("barcodes.tsv.gz", "barcodes.tsv")) {
          path <- file.path(subdir, pattern)
          if (file.exists(path)) {
            cells_file <- path
            break
          }
        }
        
        # Check if all files exist in the subdirectory
        if (!is.null(mtx_file) && !is.null(features_file) && !is.null(cells_file)) {
          tryCatch({
            # Read the matrix
            counts <- Seurat::ReadMtx(
              mtx = mtx_file,
              features = features_file,
              cells = cells_file
            )
            
            # Create Seurat object
            seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = subdir_name)
            
            # Calculate percent mitochondrial
            seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
            
            # Add to list
            seurat_list[[subdir_name]] <- seurat_obj
          }, error = function(e) {
            showNotification(paste("Error processing folder", subdir_name, ":", e$message), type = "warning")
          })
        }
      }
      
      # Check if any datasets were found
      if (length(seurat_list) == 0) {
        # Try one last approach - using Read10X directly on the folder
        tryCatch({
          counts <- Seurat::Read10X(data.dir = folder_path)
          
          # Create Seurat object
          seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = basename(folder_path))
          
          # Calculate percent mitochondrial
          seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
          
          # Add to list
          seurat_list[[basename(folder_path)]] <- seurat_obj
        }, error = function(e) {
          # Finally, if nothing works, show error
          showNotification("No valid datasets found in the uploaded folder", type = "error")
          return(NULL)
        })
      }
      
      # Combine datasets using Seurat v5's layer approach
      if (length(seurat_list) == 1) {
        # Only one dataset found
        combined_seurat <- seurat_list[[1]]
        
        # Add origin column
        combined_seurat$dataset_origin <- names(seurat_list)[1]
      } else {
        # First create a merged object with unique cell names
        cell_ids <- names(seurat_list)
        combined_seurat <- merge(
          seurat_list[[1]], 
          y = seurat_list[-1],
          add.cell.ids = cell_ids,
          project = "combined"
        )
        
        # Add a dataset origin metadata column
        combined_seurat$dataset_origin <- NA
        for (i in seq_along(cell_ids)) {
          dataset_name <- cell_ids[i]
          # Match cells from this dataset by their prefix
          cells_from_dataset <- grep(paste0("^", dataset_name, "_"), colnames(combined_seurat), value = TRUE)
          combined_seurat$dataset_origin[colnames(combined_seurat) %in% cells_from_dataset] <- dataset_name
        }
      }
      
      # Store in reactive values
      values$seurat <- combined_seurat
      values$original_seurat <- combined_seurat
      
      # Show success notification
      showNotification(paste("Loaded", length(seurat_list), "datasets successfully!"), type = "message")
    })
    
    # Show data summary
    output$dataSummary <- renderPrint({
      req(values$seurat)
      cat("Dataset loaded successfully.\n")
      cat(paste("Number of cells:", ncol(values$seurat), "\n"))
      cat(paste("Number of genes:", nrow(values$seurat), "\n"))
      cat("\nSample of gene names:\n")
      print(head(rownames(values$seurat)))
      cat("\nSample of cell barcodes:\n")
      print(head(colnames(values$seurat)))
    })
    
    # Navigate to QC tab
    updateTabItems(session, "tabs", "qc")
  })

  # Quality Control logic
  output$qcPlot <- renderPlot({
    req(values$seurat)
    VlnPlot(values$seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  })
  
  observeEvent(input$runQC, {
    req(values$seurat)
    
    withProgress(message = 'Running QC filtering...', {
      # Apply QC filtering
      filtered_seurat <- subset(values$seurat, 
                                nFeature_RNA > input$minGenes & 
                                  nFeature_RNA < input$maxGenes & 
                                  nCount_RNA > input$minFeatures &
                                  nCount_RNA < input$maxFeatures &
                                  percent.mt < input$maxMito)
      
      # Update Seurat object
      values$seurat <- filtered_seurat
      values$qc_done <- TRUE
    })
    
    # Show QC summary
    output$qcSummary <- renderPrint({
      req(values$seurat)
      cat("QC filtering completed.\n")
      cat(paste("Cells before filtering:", ncol(values$original_seurat), "\n"))
      cat(paste("Cells after filtering:", ncol(values$seurat), "\n"))
      cat(paste("Cells removed:", ncol(values$original_seurat) - ncol(values$seurat), "\n"))
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "preprocess")
  })
  
  # Preprocessing logic triggered by commit button
  observeEvent(input$commitPreprocessing, {
    # Check which methods are selected
    selected_methods <- input$preprocessMethods
    
    # If no methods are selected, show a notification and return
    if (length(selected_methods) == 0) {
      showNotification("Please select at least one preprocessing method", 
                       type = "warning")
      return()
    }
    
    # Run preprocessing based on selected methods
    preprocessing_results <- list()
    
    # Progress indication
    withProgress(message = 'Processing...', value = 0.1, {
      if ("doublet" %in% selected_methods) {
        # Perform doublet removal
        incProgress(0.3)
        preprocessing_results$doublet <- tryCatch({
          perform_doublet_removal()
        }, error = function(e) {
          showNotification(paste("Error in Doublet Removal:", e$message), 
                           type = "error")
          NULL
        })
      }
      
      if ("otherMethod" %in% selected_methods) {
        # Perform other preprocessing method
        incProgress(0.6)
        preprocessing_results$otherMethod <- tryCatch({
          perform_other_method()
        }, error = function(e) {
          showNotification(paste("Error in Other Method:", e$message), 
                           type = "error")
          NULL
        })
      }
      
      # Update the plot or results
      output$preprocessingPlot <- renderPlot({
        # Create a plot or visualization of preprocessing results
        # This will depend on your specific preprocessing methods
        if (length(preprocessing_results) > 0) {
          # Example plot logic - adjust according to your needs
          # You might want to create different plots based on selected methods
          par(mfrow = c(1, length(preprocessing_results)))
          for (method in names(preprocessing_results)) {
            # Placeholder - replace with actual plotting logic
            plot(preprocessing_results[[method]], 
                 main = paste("Results:", method))
          }
        }
      })
    })
  })

  # Normalization logic
  observeEvent(input$runNorm, {
    req(values$seurat, values$qc_done)
    
    withProgress(message = 'Running normalization...', {
      if(input$normMethod == "LogNormalize") {
        values$seurat <- NormalizeData(values$seurat, 
                                       normalization.method = "LogNormalize", 
                                       scale.factor = input$scaleFactor)
      } else if(input$normMethod == "SCTransform") {
        values$seurat <- SCTransform(values$seurat)
      }
      
      values$norm_done <- TRUE
    })
    
    # Plot distribution of normalized data
    output$normPlot <- renderPlot({
      req(values$seurat, values$norm_done)
      # Sample a few genes to show normalized expression distribution
      genes_to_plot <- sample(rownames(values$seurat), 4)
      VlnPlot(values$seurat, features = genes_to_plot, ncol = 2)
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "features")
  })
  
  # Feature selection logic
  observeEvent(input$findVarFeatures, {
    req(values$seurat, values$norm_done)
    
    withProgress(message = 'Finding variable features...', {
      values$seurat <- FindVariableFeatures(values$seurat, 
                                            selection.method = "vst", 
                                            nfeatures = input$nFeatures)
      
      values$features_done <- TRUE
    })
    
    # Plot variable features
    output$varFeaturesPlot <- renderPlot({
      req(values$seurat, values$features_done)
      VariableFeaturePlot(values$seurat)
    })
    
    # Show top variable features
    output$topFeaturesTable <- renderDT({
      req(values$seurat, values$features_done)
      top_features <- head(VariableFeatures(values$seurat), 20)
      data.frame(
        Feature = top_features,
        Mean = rowMeans(values$seurat@assays$RNA@data[top_features, ])
      )
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "dimreduce")
  })
  
  # Dimensionality reduction logic
  observeEvent(input$runPCA, {
    req(values$seurat, values$features_done)
    
    withProgress(message = 'Running PCA...', {
      values$seurat <- ScaleData(values$seurat, features = VariableFeatures(values$seurat))
      values$seurat <- RunPCA(values$seurat, features = VariableFeatures(values$seurat), npcs = input$nPCs)
      
      values$pca_done <- TRUE
    })
    
    # Plot PCA results
    output$pcaPlot <- renderPlot({
      req(values$seurat, values$pca_done)
      DimPlot(values$seurat, reduction = "pca")
    })
    
    # Plot elbow plot
    output$elbowPlot <- renderPlot({
      req(values$seurat, values$pca_done)
      ElbowPlot(values$seurat, ndims = input$nPCs)
    })
  })
  
  observeEvent(input$runReduction, {
    req(values$seurat, values$pca_done)
    
    withProgress(message = paste('Running', input$reduction, '...'), {
      if(input$reduction == "UMAP") {
        values$seurat <- RunUMAP(values$seurat, dims = 1:input$pcDims)
      } else if(input$reduction == "t-SNE") {
        values$seurat <- RunTSNE(values$seurat, dims = 1:input$pcDims)
      }
      
      values$reduction_done <- TRUE
    })
    
    # Plot reduction results
    output$reductionPlot <- renderPlot({
      req(values$seurat, values$reduction_done)
      DimPlot(values$seurat, reduction = tolower(input$reduction))
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "cluster")
  })
  
  # Clustering logic
  observeEvent(input$runClustering, {
    req(values$seurat, values$reduction_done)
    
    withProgress(message = 'Running clustering...', {
      # Find neighbors
      values$seurat <- FindNeighbors(values$seurat, dims = 1:input$pcDims)
      
      # Find clusters
      algorithm <- ifelse(input$clusterAlgo == "Louvain", 1, 4) # 1 for Louvain, 4 for Leiden
      values$seurat <- FindClusters(values$seurat, resolution = input$resolution, algorithm = algorithm)
      
      values$cluster_done <- TRUE
    })
    
    # Plot clustering results
    output$clusterPlot <- renderPlot({
      req(values$seurat, values$cluster_done)
      DimPlot(values$seurat, reduction = tolower(input$reduction), group.by = "seurat_clusters", label = TRUE)
    })
    
    # Show cluster statistics
    output$clusterStats <- renderDT({
      req(values$seurat, values$cluster_done)
      table_data <- table(values$seurat$seurat_clusters)
      data.frame(
        Cluster = names(table_data),
        Cell_Count = as.vector(table_data),
        Percentage = round(as.vector(table_data) / sum(table_data) * 100, 2)
      )
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "annotate")
  })
  
  # Cell type annotation logic
  observeEvent(input$findMarkers, {
    req(values$seurat, values$cluster_done)
    
    withProgress(message = 'Finding marker genes...', {
      # Find markers for each cluster
      all_markers <- FindAllMarkers(values$seurat, 
                                    only.pos = TRUE, 
                                    min.pct = 0.3, 
                                    logfc.threshold = 0.5)
      
      # Store markers
      values$markers <- all_markers %>%
        group_by(cluster) %>%
        top_n(n = input$maxMarkers, wt = avg_log2FC)
      
      values$markers_done <- TRUE
    })
    
    # Display marker table
    output$markerTable <- renderDT({
      req(values$markers, values$markers_done)
      values$markers
    })
    
    # Create UI for cell type assignment
    output$cellTypeUI <- renderUI({
      req(values$seurat, values$markers_done)
      clusters <- levels(values$seurat$seurat_clusters)
      
      lapply(clusters, function(cluster) {
        fluidRow(
          column(4, paste0("Cluster ", cluster)),
          column(8, textInput(paste0("label_cluster_", cluster), 
                              label = NULL, 
                              value = paste0("Cluster ", cluster)))
        )
      })
    })
  })
  
  observeEvent(input$assignCellTypes, {
    req(values$seurat, values$markers_done)
    
    withProgress(message = 'Assigning cell types...', {
      # Get cluster IDs and user-defined labels
      clusters <- levels(values$seurat$seurat_clusters)
      labels <- sapply(clusters, function(cluster) {
        input[[paste0("label_cluster_", cluster)]]
      })
      
      # Create new identity column
      new_ids <- plyr::mapvalues(values$seurat$seurat_clusters, 
                                 from = clusters, 
                                 to = labels)
      values$seurat$cell_type <- new_ids
      
      values$annotation_done <- TRUE
    })
    
    # Plot annotated clusters
    output$annotatedPlot <- renderPlot({
      req(values$seurat, values$annotation_done)
      DimPlot(values$seurat, reduction = tolower(input$reduction), group.by = "cell_type", label = TRUE)
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "download")
  })
  
  # Download results logic
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("scRNA_results_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      # Create temporary directory
      temp_dir <- tempdir()
      
      # Prepare files based on user selection
      if("rds" %in% input$downloadItems) {
        saveRDS(values$seurat, file.path(temp_dir, "seurat_object.rds"))
      }
      
      if("meta" %in% input$downloadItems) {
        write.csv(values$seurat@meta.data, file.path(temp_dir, "cell_metadata.csv"))
      }
      
      if("markers" %in% input$downloadItems && !is.null(values$markers)) {
        write.csv(values$markers, file.path(temp_dir, "marker_genes.csv"))
      }
      
      if("dimred" %in% input$downloadItems && !is.null(values$reduction_done)) {
        # Get the last used reduction
        red_type <- tolower(input$reduction)
        if(red_type %in% names(values$seurat@reductions)) {
          coords <- as.data.frame(values$seurat@reductions[[red_type]]@cell.embeddings)
          write.csv(coords, file.path(temp_dir, paste0(red_type, "_coordinates.csv")))
        }
      }
      
      if("report" %in% input$downloadItems) {
        # Generate report (would need additional code for comprehensive report)
        # Simple placeholder
        cat("# scRNA-seq Analysis Report\n\nDate: ", Sys.Date(), 
            "\n\n## Summary\n\n", file = file.path(temp_dir, "report.md"))
      }
      
      # Zip all files
      files_to_zip <- list.files(temp_dir, pattern = "\\.csv$|\\.rds$|\\.md$", full.names = TRUE)
      zip(file, files_to_zip)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)