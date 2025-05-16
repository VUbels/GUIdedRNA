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
library(DoubletFinder)
library(decontX)
#library(GUIdedRNA)

source("~/GUIdedRNA/GUIdedRNA/global.R")
source("~/GUIdedRNA/R/preprocessing_functions.R")

options(shiny.maxRequestSize = 100 * 1024^3)

# Global message queue 
message_queue <- character(0)

# Global message handling functions
send_message <- function(msg) {
  # Add to queue
  message_queue <<- c(message_queue, msg)
  # Also log to console
  message(paste("Message queued:", msg))
}

get_messages <- function() {
  msgs <- message_queue
  message_queue <<- character(0)
  return(msgs)
}

ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = "GUIdedRNA"),
  
  # Sidebar with analysis steps
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Setup", tabName = "upload", icon = icon("upload"), selected = TRUE),
      menuItem("Sample Information", tabName = "information", icon = icon("list")),
      menuItem("Quality Control", tabName = "qc", icon = icon("check-circle")),
      menuItem("Preprocessing", tabName = "preprocess", icon = icon("filter")),
      menuItem("LSI Round 1", tabName = "LSI_1", icon = icon("project-diagram")),
      menuItem("Clustering", tabName = "cluster", icon = icon("object-group")),
      menuItem("Cell Type Annotation", tabName = "annotate", icon = icon("tags")),
      menuItem("Download Results", tabName = "download", icon = icon("download"))
    )
  ),
  
  # Main panel with tabbed content
  dashboardBody(
    useShinyjs(),
    extendShinyjs(
      text = "shinyjs.appendToConsole = function(message) {
    var consoleOutput = document.getElementById('consoleOutput');
    if(consoleOutput) {
      var p = document.createElement('p');
      p.innerHTML = message;
      p.style.fontSize = '0.75em';  // This makes the font smaller
      p.style.margin = '0';         // Remove default paragraph margin for compact display
      p.style.padding = '2px 0';    // Add a little vertical padding between lines
      consoleOutput.appendChild(p);
      consoleOutput.scrollTop = consoleOutput.scrollHeight;
    }
  }",
      functions = c("appendToConsole")
    ),
    
    tabItems(
      # Upload Data tab
      tabItem(tabName = "upload",
              
              fluidRow(
                box(
                  title = "Set Output Directory",
                  width = 12,
                  shinyDirButton("folderBtn_output", "Select Directory", "Choose an output directory"),
                  helpText("Select directory that will contain all output"),
                  textOutput("outputFolder"),
                  actionButton("outputData_folder", "Set Output Folder")
                ),
              ),
              
              fluidRow(
                box(
                  title = "Upload Folder or list of Folders for multiple 10X Genomic datasets",
                  width = 12,
                  shinyDirButton("folderBtn_single", "Select Directory", "Choose a directory with 10X data"),
                  helpText("Select directory containing files, or directory containing dataset folders"),
                  textOutput("selectedFolder"),
                  actionButton("loadData_folder", "Load Data")
                ),
              ),
              
              fluidRow(
                box(
                  title = "Upload 10X Genomics files for Single Dataset",
                  width = 12,
                  fileInput("matrixFile", "Matrix (.mtx)"),
                  fileInput("featuresFile", "Features/Genes (.tsv/.txt)"),
                  fileInput("barcodesFile", "Barcodes (.tsv/.txt)"),
                  actionButton("loadData_file", "Load Data")
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
      
      # Setting sample information
      tabItem(tabName = "information",
              fluidRow(
                box(
                  title = "Sample Information Manager",
                  width = 12,
                  solidHeader = TRUE,
                  
                  fluidRow(
                    # Column 1: Display original identities from all objects
                    column(
                      width = 4,
                      div(
                        class = "panel panel-default",
                        div(class = "panel-heading", h4("Original Sample IDs")),
                        div(
                          class = "panel-body", 
                          style = "max-height: 400px; overflow-y: auto;",
                          uiOutput("originalIdentities")
                        )
                      )
                    ),
                    
                    # Column 2: Add new sample attributes
                    column(
                      width = 8,
                      div(
                        class = "panel panel-default",
                        div(class = "panel-heading", h4("Create New Sample Attributes")),
                        div(
                          class = "panel-body",
                          
                          # Input for new column name
                          textInput("newColumnName", "New Column Name", value = "sample_type"),
                          
                          # Dynamic UI for assigning values to samples
                          uiOutput("sampleValuesInput"),
                          
                          # Button to create the new column
                          actionButton("createNewColumn", "Create New Column", 
                                       class = "btn-success", style = "margin-top: 15px;"),
                          
                          # Display added columns
                          tags$hr(),
                          h4("Added Columns"),
                          uiOutput("addedColumnsUI")
                        )
                      )
                    )
                  )
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
                    column(6, numericInput("minFeatures", "Min Genes", min = 0, value = 200, step = 100)),
                    column(6, numericInput("maxFeatures", "Max Genes", min = 0, value = 5000, step = 100)),
                    helpText("Determines minimum and maximum variety of genes in data")
                  ),
                  
                  # Side by side feature inputs
                  fluidRow(
                    column(6, numericInput("minCount", "Min Count", min = 0, value = 100, step = 100)),
                    column(6, numericInput("maxCount", "Max Count", min = 0, value = 3000, step = 100)),
                    helpText("Determines minimum and maximum total RNA count in data")
                  ),
                  
                  # Single mitochondrial percentage input
                  numericInput("maxMito", "Max Mitochondrial %", min = 0, max = 100, value = 20),
                  helpText("Determines minimum and maximum mitochondrial RNA in data"),
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
                  # Add these style attributes for text alignment
                  div(id = "consoleOutput", 
                      style = "white-space: pre-wrap; height: 600px; overflow-y: auto; background-color: #f5f5f5; padding: 1px; font-family: monospace; text-align: left; vertical-align: top;")
                )
              )
            ),
      
      tabItem(tabName = "LSI_1",
              fluidRow(
                box(
                  title = "LSI Round 1",
                  width = 4,
                  checkboxGroupInput("blacklist_genes", "Ignore Mitochondrial, X/Y, Ribosomal genes from variable genes", 
                                     choices = c("Ignore Mitochondrial Genes" = "blacklist_mitogenes", 
                                                 "Ignore X/Y chomosomal genes" = "blacklist_sexgenes",
                                                 "Ignore Ribosomal genes" = "blacklist_rbgenes"),
                                     
                                     
                                     
                                     selected = c("blacklist_mitogenes", "blacklist_sexgenes", "blacklist_rbgenes")
                                     
                  ),
                  actionButton("commitLSI_1", "Run LSI with selected Methods", 
                               class = "btn-success", 
                               style = "width: 100%;")
                ),
                box(
                  title = "LSI Results",
                  width = 8,
                  # Add these style attributes for text alignment
                  div(id = "consoleOutput", 
                      style = "white-space: pre-wrap; height: 600px; overflow-y: auto; background-color: #f5f5f5; padding: 1px; font-family: monospace; text-align: left; vertical-align: top;")
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


#Defining server logic
server <- function(input, output, session) {

  # Initialize log text reactive value
  log_text <- reactiveVal("")
  
  
  display_ui_message <- function(msg) {
    # Function to display messages in the UI console
    # Update the log text
    current <- isolate(log_text())
    log_text(paste0(current, msg, "\n"))
    
    # Send to JavaScript console
    js$appendToConsole(msg)
  }
  
  # Create an observer that checks for new messages
  observe({
    msgs <- get_messages()
    if (length(msgs) > 0) {
      # Process all queued messages
      for (msg in msgs) {
        display_ui_message(msg)
      }
    }
    invalidateLater(100) 
  })
  
  
  assign("send_message", function(msg) {
    display_ui_message(msg)
  }, envir = .GlobalEnv)
  
  # Render log text
  output$preprocessingLog <- renderText({
    log_text()
  })
  
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
  
  # Initialize directory chooser for input
  shinyDirChoose(input, "folderBtn_single", roots = volumes, session = session)
  
  # Display selected input folder
  output$selectedFolder <- renderText({
    if (!is.null(input$folderBtn_single)) {
      parseDirPath(volumes, input$folderBtn_single)
      showNotification(paste("Input folder set to:", input$folderBtn_single), type = "message")
    } else {
      "No folder selected"
    }
  })
  
  # Initialize directory chooser for output
  shinyDirChoose(input, "folderBtn_output", roots = volumes, session = session)
  
  output$outputFolder <- renderText({
    if (!is.null(output_directory())) {
      output_directory()
    } else {
      "No folder selected"
    }
  })
  
  # Display selected output folder
  output_directory <- reactive({
    req(input$folderBtn_output)
    parseDirPath(volumes, input$folderBtn_output)
  })
  
  observeEvent(input$outputData_folder, {
    req(output_directory())
    output_folder <<- output_directory()
    showNotification(paste("Output folder set to:", output_folder), type = "message")
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
      values$seurat_list <- seurat_obj
      values$original_seurat_list <- seurat_obj
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
            
            # Calculate percent mitochondrial either in human or other data (murine/rat)
            seurat_obj[["percent.mt"]] <- if(any(grepl("^MT-", rownames(seurat_obj)))) {
              Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
            } else if(any(grepl("^mt-", rownames(seurat_obj)))) {
              Seurat::PercentageFeatureSet(seurat_obj, pattern = "^mt-")
            } else {
              # If neither pattern exists, create a zero vector as fallback
              rep(0, ncol(seurat_obj))
            }
            
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
      
      # Store the list in reactive values - instead of combining objects
      values$seurat_list <- seurat_list
      values$original_seurat_list <- seurat_list
      
      # Show success notification
      showNotification(paste("Loaded", length(seurat_list), "datasets successfully!"), type = "message")
    })
    
    # Show data summary
    output$dataSummary <- renderPrint({
      req(values$seurat_list)
      cat("Datasets loaded successfully.\n")
      cat(paste("Number of datasets:", length(values$seurat_list), "\n\n"))
      
      # Print summary information for each dataset
      for (dataset_name in names(values$seurat_list)) {
        seurat_obj <- values$seurat_list[[dataset_name]]
        cat(paste0("Dataset: ", dataset_name, "\n"))
        cat(paste("  Number of cells:", ncol(seurat_obj), "\n"))
        cat(paste("  Number of genes:", nrow(seurat_obj), "\n"))
        cat("\n")
      }
      
      # Show sample of gene names from first dataset
      if (length(values$seurat_list) > 0) {
        first_dataset <- values$seurat_list[[1]]
        cat("Sample of gene names from first dataset:\n")
        print(head(rownames(first_dataset)))
        cat("\nSample of cell barcodes from first dataset:\n")
        print(head(colnames(first_dataset)))
      }
    })
    
    # Navigate to QC tab
    updateTabItems(session, "tabs", "information")
  })
  
  # Logic to assign sample information to file list
  output$originalIdentities <- renderUI({
    req(values$seurat)
    
    sample_ids <- levels(values$seurat$orig.ident)
    
    tagList(
      tags$div(
        class = "sample-list",
        lapply(sample_ids, function(id) {
          tags$div(
            class = "sample-id-item",
            style = "padding: 8px; border-bottom: 1px solid #eee;",
            tags$span(class = "label label-default", style = "margin-right: 10px;", id),
            id
          )
        })
      )
    )
  })
  
  getAllSampleIDs <- function() {
    req(values$seurat_list)
    
    # Collect all unique orig.ident values from all Seurat objects
    all_ids <- c()
    for (seurat_obj in values$seurat_list) {
      sample_ids <- levels(seurat_obj$orig.ident)
      all_ids <- c(all_ids, sample_ids)
    }
    
    # Return unique IDs
    return(unique(all_ids))
  }
  
  output$originalIdentities <- renderUI({
    req(values$seurat_list)
    
    sample_ids <- getAllSampleIDs()
    
    tagList(
      tags$div(
        class = "sample-list",
        lapply(sample_ids, function(id) {
          # Count how many Seurat objects contain this sample ID
          count <- sum(sapply(values$seurat_list, function(seurat_obj) {
            id %in% levels(seurat_obj$orig.ident)
          }))
          
          tags$div(
            class = "sample-id-item",
            style = "padding: 8px; border-bottom: 1px solid #eee;",
            tags$span(class = "label label-default", style = "margin-right: 10px;", id),
            span(id),
            tags$span(
              class = "pull-right",
              style = "color: #888; font-size: 0.9em;",
              paste0("(in ", count, " object", ifelse(count > 1, "s", ""), ")")
            )
          )
        })
      )
    )
  })
  
  # Create input fields for assigning values to each sample
  output$sampleValuesInput <- renderUI({
    req(values$seurat_list)
    
    sample_ids <- getAllSampleIDs()
    
    tagList(
      tags$div(
        style = "max-height: 300px; overflow-y: auto; border: 1px solid #ddd; padding: 10px; margin-top: 10px;",
        lapply(sample_ids, function(id) {
          fluidRow(
            column(4, 
                   tags$label(class = "control-label", `for` = paste0("value_", id), id)
            ),
            column(8, 
                   textInput(
                     inputId = paste0("value_", id),
                     label = NULL,
                     value = "",
                     placeholder = paste("Value for", id)
                   )
            )
          )
        })
      )
    )
  })
  
  # Logic to create a new column
  observeEvent(input$createNewColumn, {
    req(values$seurat_list, input$newColumnName)
    
    # Validate the column name isn't empty
    if (trimws(input$newColumnName) == "") {
      showNotification("Column name cannot be empty", type = "error")
      return()
    }
    
    # Check if column already exists in any Seurat object
    column_exists <- FALSE
    for (seurat_obj in values$seurat_list) {
      if (input$newColumnName %in% colnames(seurat_obj@meta.data)) {
        column_exists <- TRUE
        break
      }
    }
    
    if (column_exists) {
      showModal(modalDialog(
        title = "Column already exists",
        "This column already exists in one or more Seurat objects. Do you want to overwrite it?",
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirmOverwrite", "Overwrite")
        )
      ))
      return()
    }
    
    # Create the column
    createColumn()
  })
  
  # Handle overwrite confirmation
  observeEvent(input$confirmOverwrite, {
    removeModal()
    createColumn()
  })
  
  # Function to create/update column in all Seurat objects
  createColumn <- function() {
    withProgress(message = 'Creating new column in all Seurat objects...', {
      # Get all sample IDs
      sample_ids <- getAllSampleIDs()
      
      # Get user-defined values for each sample
      values_map <- sapply(sample_ids, function(id) {
        input[[paste0("value_", id)]]
      })
      names(values_map) <- sample_ids
      
      # Update each Seurat object in the list
      for (i in seq_along(values$seurat_list)) {
        incProgress(1/length(values$seurat_list), 
                    detail = paste("Processing object", i, "of", length(values$seurat_list)))
        
        # Get sample IDs for this specific Seurat object
        obj_sample_ids <- levels(values$seurat_list[[i]]$orig.ident)
        
        # Filter values map for just this object's samples
        obj_values <- values_map[obj_sample_ids]
        
        # Create new metadata column for this object
        new_values <- plyr::mapvalues(values$seurat_list[[i]]$orig.ident, 
                                      from = obj_sample_ids, 
                                      to = obj_values)
        
        # Add or update the column in this object's metadata
        values$seurat_list[[i]][[input$newColumnName]] <- new_values
      }
      
      # Track added columns if not already doing so
      if(is.null(values$added_columns)) {
        values$added_columns <- c(input$newColumnName)
      } else {
        # Only add if it's not already in the list
        if(!(input$newColumnName %in% values$added_columns)) {
          values$added_columns <- c(values$added_columns, input$newColumnName)
        }
      }
      
      # Reset the column name input
      updateTextInput(session, "newColumnName", value = "")
      
      # Show success message
      showNotification(paste("Column", input$newColumnName, "created successfully in all Seurat objects!"), type = "message")
    })
  }
  
  # Display added columns
  output$originalIdentities <- renderUI({
    req(values$seurat_list)
    
    sample_ids <- getAllSampleIDs()
    
    tagList(
      tags$div(
        class = "sample-list",
        lapply(sample_ids, function(id) {
          tags$div(
            class = "sample-id-item",
            style = "padding: 8px; border-bottom: 1px solid #eee;",
            span(id)
          )
        })
      )
    )
  })
  
  # Handle edit button clicks for each column
  observe({
    req(values$added_columns, values$seurat_list)
    
    lapply(values$added_columns, function(col_name) {
      observeEvent(input[[paste0("editColumn_", col_name)]], {
        # Set the column name in the input
        updateTextInput(session, "newColumnName", value = col_name)
        
        # Get all unique sample IDs
        sample_ids <- getAllSampleIDs()
        
        # For each sample ID, find its value in the first Seurat object where it exists
        for(id in sample_ids) {
          # Find the first Seurat object that contains this sample ID
          for(seurat_obj in values$seurat_list) {
            if(id %in% levels(seurat_obj$orig.ident)) {
              # Get the current value for this sample
              current_value <- seurat_obj@meta.data[seurat_obj$orig.ident == id, col_name][1]
              
              # Update the text input
              updateTextInput(session, paste0("value_", id), value = current_value)
              
              # Move to the next sample ID
              break
            }
          }
        }
      })
    })
  })
  
  # Quality Control logic
  output$qcPlot <- renderPlot({
    req(values$seurat_list)
  
    # Combine meta.data across all Seurat objects with dataset label
    meta_combined <- do.call(rbind, lapply(names(values$seurat_list), function(name) {
      obj <- values$seurat_list[[name]]
      meta <- obj@meta.data
      meta$dataset <- name
      return(meta)
    }))
    
    # Plot using ggplot2
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    plots <- lapply(features, function(feat) {
      ggplot(meta_combined, aes_string(x = "dataset", y = feat, fill = "dataset")) +
        geom_violin(trim = FALSE) +
        theme_minimal() +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(feat)
    })
    
    # Display in a grid
    cowplot::plot_grid(plotlist = plots, ncol = 3)
  })
    
    
  observeEvent(input$runQC, {
    req(values$seurat_list)
    
    withProgress(message = 'Running QC filtering...', {
      # Filter each object individually
      filtered_list <- lapply(values$seurat_list, function(seurat_obj) {
        subset(seurat_obj,
               nFeature_RNA > input$minFeatures &
                 nFeature_RNA < input$maxFeatures &
                 nCount_RNA > input$minCount &
                 nCount_RNA < input$maxCount &
                 percent.mt < input$maxMito)
      })
      
      values$seurat_list <- filtered_list
      values$qc_done <- TRUE
    })
    
    # QC Summary
    output$qcSummary <- renderPrint({
      req(values$seurat_list, values$original_seurat_list)
      cat("QC filtering completed.\n")
      
      for (name in names(values$original_seurat_list)) {
        orig_cells <- ncol(values$original_seurat_list[[name]])
        filtered_cells <- ncol(values$seurat_list[[name]])
        cat(paste0("Dataset: ", name, "\n"))
        cat(paste("  Cells before filtering:", orig_cells, "\n"))
        cat(paste("  Cells after filtering:", filtered_cells, "\n"))
        cat(paste("  Cells removed:", orig_cells - filtered_cells, "\n\n"))
      }
    })
    
    updateTabItems(session, "tabs", "preprocess")
  })
 
  # Preprocessing logic triggered by commit button
  observeEvent(input$commitPreprocessing, {
    req(values$seurat_list)
    
    # Reset the log text
    log_text("")
    
    # Check which methods are selected
    selected_methods <- input$preprocessMethods
    
    # If no methods are selected, show a notification and return
    if (length(selected_methods) == 0) {
      showNotification("Please select at least one preprocessing method", 
                       type = "warning")
      return()
    }
    
    # Send initial message
    send_message("Starting preprocessing...")
    
    # Get the current Seurat list
    preprocessing_seurat_list <- values$seurat_list
    
    # Progress indication
    withProgress(message = 'Processing...', value = 0, {
      if ("doublet" %in% selected_methods) {
        incProgress(0.1, detail = "Running Doublet Removal...")
        
        send_message("Starting Doublet Removal, this may take some time...")
        
        preprocessing_seurat_list <- tryCatch({
          preprocess_DoubletRemoval(preprocessing_seurat_list)
        }, error = function(e) {
          send_message(paste("Error in Doublet Removal:", e$message))
          return(preprocessing_seurat_list)
        })
        
        incProgress(0.5)
      }
      
      if ("ambient" %in% selected_methods) {
        incProgress(0.1, detail = "Running Ambient RNA Removal...")
        
        send_message("Starting Ambient RNA Removal...")
        
        preprocessing_seurat_list <- tryCatch({
          preprocess_AmbientRNA(preprocessing_seurat_list)
        }, error = function(e) {
          send_message(paste("Error in Ambient RNA Removal:", e$message))
          return(preprocessing_seurat_list)
        })
        
        send_message("Filtering low RNA cells...")
        
        preprocessing_seurat_list <- tryCatch({
          remove_lowRNA(preprocessing_seurat_list)
        }, error = function(e) {
          send_message(paste("Error in Low RNA Filtering:", e$message))
          return(preprocessing_seurat_list)
        })
        
        incProgress(0.4)
      }
      
      # Update reactive values
      values$seurat_list <- preprocessing_seurat_list
      merged_seurat <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
      merged_seurat$orig.ident <- as.factor(merged_seurat$orig.ident)
      values$seurat <- merged_seurat
      values$preprocess_done <- TRUE
      send_message("Preprocessing and merging completed successfully!")
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "information")
  })

  
  
  observeEvent(input$commitLSI_1, {
    req(values$seurat)
    
    # Reset the log text
    log_text("")
    
    # Check which methods are selected
    selected_methods <- input$preprocessMethods
    
    # If no methods are selected, show a notification and return
    if (length(selected_methods) == 0) {
      showNotification("Running without ignoring mitotic cycle genes, X/Y chromosome genes and mitochondrial genes", 
                       type = "warning")
    }
    
    
    # Send initial message
    send_message("Starting log normalization...")
    
    
    LSI1_seurat <- values$seurat
    rawCounts <- GetAssayData(object=seurat, layer="counts")
    log2CP10k <- sparseLogX(rawCounts, logtype="log2", scale=TRUE, scaleFactor=10^4)
    LSI_seurat <- SetAssayData(object=LSI_seurat, layer="data", new.data=log2CP10k)
    
    umapNeighbors <- 50
    umapMinDist <- 0.5
    umapDistMetric <- "cosine"
    
  
    message("Running iterative LSI...")
    set.seed(1)
    
    #Initial Cluster allocation
    
    for(i in seq_along(resolution)){
      # If first round, compute variable genes on raw data first
      if(i == 1){
        message(sprintf("Identifying top %s variable genes among all cells...", nVarGenes))
        varGenes <- getVarGenes(log2CP10k, nvar=nVarGenes, blacklist=blacklist.genes)
      }else{
        # For remaining rounds, calculate variable genes using previous clusters
        clusterMat <- edgeR::cpm(groupSums(rawCounts, clusters, sparse=TRUE), log=TRUE, prior.count=3)
        message(sprintf("Identifying top %s variable genes from round %s LSI...", nVarGenes, i-1))
        varGenes <- getVarGenes(clusterMat, nvar=nVarGenes, blacklist=blacklist.genes)
      }
      # Now run LSI and find clusters
      LSIi <- runLSI(rawCounts[varGenes,], nComponents = max(nPCs), binarize = FALSE)
      
        message(sprintf("Harmonizing LSI SVD PCs for round %s...", i))
        harmonized_pcs <- HarmonyMatrix(
          data_mat  = LSIi$matSVD,
          meta_data = seurat_obj_merged@meta.data,
          vars_use  = covariates, # Covariates to 'harmonize'
          do_pca    = FALSE
        )
        LSIi$matSVD <- harmonized_pcs
      }
      
      reducName <- paste0("LSI_iter",i)
      seurat_obj_merged[[reducName]] <- CreateDimReducObject(embeddings = LSIi$matSVD, key = sprintf("LSI%s_", i), assay = "RNA")
      seurat_obj_merged <- FindNeighbors(object = seurat_obj_merged, reduction = reducName, dims = nPCs, force.recalc = TRUE)
      message(sprintf("Clustering with resolution %s...", resolution[i]))
      seurat_obj_merged <- FindClusters(object = seurat_obj_merged, resolution = resolution[i], algorithm = 4, method = "igraph")
      clusters <- Idents(seurat_obj_merged)
      #Store information
      lsiOut[[reducName]] <- list(
        lsiMat = LSIi$matSVD,
        svd = LSIi$svd,
        varFeatures = varGenes, 
        clusters = clusters
      )
    }
  )
  
  
  # Dimensionality reduction logic
  observeEvent(input$runPCA, {
    req(values$seurat, values$features_done)
    
    withProgress(message = 'Running PCA...', {
      values$seurat <- 
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