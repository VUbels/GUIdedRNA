# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)

message("Starting app.R...")
if (file.exists("global.R")) {
  message("Found global.R, loading it...")
  source("global.R")
} else {
  message("global.R not found in working directory:", getwd())
}

#library(GUIdedRNA)

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
      
      tags$li(class = "divider", style = "height: 1px; margin: 9px 0; overflow: hidden; background-color: #666;"),
      
      menuItem("Sample Information", tabName = "information", icon = icon("list")),
      menuItem("Quality Control", tabName = "qc", icon = icon("check-circle")),
      menuItem("Preprocessing", tabName = "preprocess", icon = icon("filter")),
      
      tags$li(class = "divider", style = "height: 1px; margin: 9px 0; overflow: hidden; background-color: #666;"),
      
      menuItem("LSI Round 1", tabName = "LSI_1", icon = icon("project-diagram")),
      menuItem("Initial Clustering", tabName = "annotate_broad", icon = icon("object-group")),
      
      tags$li(class = "divider", style = "height: 1px; margin: 9px 0; overflow: hidden; background-color: #666;"),
      
      menuItem("LSI Round 2", tabName = "LSI_2", icon = icon("project-diagram")),
      menuItem("Final Clustering", tabName = "annotate_specific", icon = icon("object-group")),
      
      tags$li(class = "divider", style = "height: 1px; margin: 9px 0; overflow: hidden; background-color: #666;"),
      
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
  };
  
  shinyjs.appendToLSI1_Console = function(message) {
    var consoleOutput = document.getElementById('lsi1_ConsoleOutput');
    if(consoleOutput) {
      var p = document.createElement('p');
      p.innerHTML = message;
      p.style.fontSize = '0.75em';
      p.style.margin = '0';
      p.style.padding = '2px 0';
      consoleOutput.appendChild(p);
      consoleOutput.scrollTop = consoleOutput.scrollHeight;
    }
  }",
      functions = c("appendToConsole", "appendToLSI1_Console")
    ),
    
    tabItems(
      # Upload Data tab
      tabItem(tabName = "upload",
              
              fluidRow(
                
                box(
                  title = "Set Output Directory",
                  width = 12,
                  tags$div(
                    style = "width: 20%; display: inline-block;",
                    shinyDirButton("folderBtn_output", "Select Directory", "Choose an output directory")),
                  helpText("Select directory that will contain all output"),
                  textOutput("outputFolder"),
                  actionButton("outputData_folder", "Set Output Folder",
                               class = "btn-success", 
                               style = "color: white; width: 20%;")
                ),
              ),
              
              fluidRow(
                box(
                  title = "Upload Folder or list of Folders for multiple 10X Genomic datasets",
                  width = 12,
                  shinyDirButton("folderBtn_single", "Select Directory", "Choose a directory with 10X data"),
                  helpText("Select directory containing files, or directory containing dataset folders"),
                  textOutput("selectedFolder"),
                  actionButton("loadData_folder", "Load Data",
                               class = "btn-success", 
                               style = "color: white; width: 20%;")
                ),
              ),
              
              fluidRow(
                box(
                  title = "Upload 10X Genomics files for Single Dataset",
                  width = 12,
                  fileInput("matrixFile", "Matrix (.mtx)"),
                  fileInput("featuresFile", "Features/Genes (.tsv/.txt)"),
                  fileInput("barcodesFile", "Barcodes (.tsv/.txt)"),
                  actionButton("loadData_file", "Load Data",
                               class = "btn-success", 
                               style = "color: white; width: 20%;")
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
                                       class = "btn-success", style = "color: white; width: 100%; margin-top: 15px;"),
                          
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
                  solidHeader = TRUE,
                  status = "primary",
                  height = "650px",
                  
                  # Side by side gene count inputs
                  fluidRow(
                    column(6, numericInput("minFeatures", "Min Genes", min = 0, value = 200, step = 100)),
                    column(6, numericInput("maxFeatures", "Max Genes", min = 0, value = 5000, step = 100))
                  ),
                  helpText(style = "margin-bottom: 15px; font-size: 0.9em;", 
                           "Determines minimum and maximum variety of genes in data"),
                  
                  # Side by side feature inputs
                  fluidRow(
                    column(6, numericInput("minCount", "Min Count", min = 0, value = 100, step = 100)),
                    column(6, numericInput("maxCount", "Max Count", min = 0, value = 3000, step = 100))
                  ),
                  helpText(style = "margin-bottom: 15px; font-size: 0.9em;", 
                           "Determines minimum and maximum total RNA count in data"),
                  
                  # Single mitochondrial percentage input
                  numericInput("maxMito", "Max Mitochondrial %", min = 0, max = 100, value = 20),
                  helpText(style = "margin-bottom: 15px; font-size: 0.9em;", 
                           "Determines minimum and maximum mitochondrial RNA in data"),
                  
                  # Add some padding above the button
                  div(style = "padding-top: 10px;",
                      actionButton("runQC", "Run QC",
                                   class = "btn-success", 
                                   style = "color: white; width: 100%;")
                  )
                ),
                box(
                  title = "QC Metrics",
                  width = 8,
                  height = "650px",
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
                  status = "primary",
                  height = "650px",
                  solidHeader = TRUE,
                  checkboxGroupInput("preprocessMethods", "Select Preprocessing Methods", 
                                     choices = c("Doublet Removal" = "doublet", 
                                                 "Ambient RNA Removal" = "ambient"),
                                     selected = c("doublet", "ambient")
                                     
                  ),
                  actionButton("commitPreprocessing", "Run Selected Methods", 
                               class = "btn-success", 
                               style = "color: white; width: 100%;")
                ),
                box(
                  title = "Preprocessing Results",
                  width = 8,
                  height = "650px",
                  # Add these style attributes for text alignment
                  div(id = "consoleOutput", 
                      style = "white-space: pre-wrap; height: 600px; overflow-y: auto; background-color: #f5f5f5; padding: 1px; font-family: monospace; text-align: left; vertical-align: top;")
                )
              )
            ),
      
      tabItem(tabName = "LSI_1",
              fluidRow(
                column(
                  width = 4,
                  height = "900px",
                  box(
                    title = "Blacklist Settings",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    checkboxGroupInput("blacklist_genes", "Ignore genes for LSI clustering (Not DE analysis)", 
                                       choices = c("Ignore Mitochondrial Genes" = "blacklist_mitogenes", 
                                                   "Ignore X/Y chomosomal genes" = "blacklist_sexgenes",
                                                   "Ignore Ribosomal genes" = "blacklist_rbgenes"),
                                       selected = c("blacklist_mitogenes", "blacklist_sexgenes", "blacklist_rbgenes")
                                       
                    )
                  ),
                  
                  box(
                    title = "LSI Parameters",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    
                    numericInput("nVarGenes", "Number of Variable Genes", 
                                 value = 2000, min = 100, max = 10000, step = 100),
                    
                    sliderInput("nPCs", "Number of Principal Components", 
                                value = 25, min = 10, max = 100, step = 5),
                    
                    textInput("resolution", "Clustering Resolutions (comma separated)", 
                              value = "0.2, 0.4, 0.6"),
                    
                    textInput("covariates", "Harmony Covariates", 
                              value = "orig.ident"),
                    
                    # UMAP parameters
                    numericInput("umapNeighbors", "UMAP Neighbors", 
                                 value = 50, min = 5, max = 200, step = 5),
                    
                    numericInput("umapMinDist", "UMAP Minimum Distance", 
                                 value = 0.5, min = 0.01, max = 0.99, step = 0.01),
                    
                    selectInput("umapDistMetric", "UMAP Distance Metric",
                                choices = c("cosine", "euclidean", "manhattan", "hamming"),
                                selected = "cosine"),
                    
                    actionButton("commitLSI_1", "Run LSI - Round 1", 
                                 class = "btn-success", 
                                 style = "color: white; width: 100%;")
                  )
                ),
                
                box(
                  title = "LSI Round 1 Results",
                  width = 8,
                  height = "930px",
                  div(id = "lsi1_ConsoleOutput", 
                      style = "white-space: pre-wrap; height: 850px; overflow-y: auto; background-color: #f5f5f5; padding: 1px; font-family: monospace; text-align: left; vertical-align: top;")
                )
              )
      ),
      
      
      # Cell Type Annotation tab
      tabItem(tabName = "annotate_broad",
              fluidRow(
                # Left column - DE genes and markers
                column(
                  width = 4,
                  box(
                    title = "Marker Gene Analysis",
                    width = NULL, # Full width of column
                    status = "primary",
                    solidHeader = TRUE,
                    
                    # Slider instead of numeric input for max markers
                    sliderInput("maxMarkers", "Max Markers per Cluster", 
                                min = 0, max = 30, value = 15, step = 1),
                    
                    # Action button with consistent styling
                    div(style = "margin-top: 15px; margin-bottom: 15px;",
                        actionButton("findMarkers", "Find Markers",
                                     class = "btn-success", 
                                     style = "color: white; width: 100%;")
                    )
                  ),
                  
                  box(
                    title = "Differentially Expressed Genes",
                    width = NULL, # Full width of column
                    height = "500px",
                    status = "primary",
                    solidHeader = TRUE,
                    
                    # Cluster selection dropdown
                    selectInput("clusterSelect", "Select Cluster", 
                                choices = NULL, # Will be populated dynamically
                                width = "100%"),
                    
                    # Marker gene table with scroll
                    div(style = "height: 400px; overflow-y: auto;",
                        DTOutput("markerTable") %>% withSpinner()
                    )
                  )
                ),
                
                # Right column - UMAP visualization
                column(
                  width = 8,
                  box(
                    title = "Cluster Visualization",
                    width = NULL, # Full width of column
                    height = "657px", # Match left column height
                    status = "primary",
                    solidHeader = TRUE,
                    
                    # Feature selection for UMAP coloring
                    fluidRow(
                      column(6, 
                             selectInput("umapDisplayType", "Display Type",
                                         choices = c("Clusters" = "clusters", 
                                                     "Gene Expression" = "gene"),
                                         selected = "clusters",
                                         width = "100%")
                      ),
                      column(6,
                             # Conditional UI for gene selection
                             conditionalPanel(
                               condition = "input.umapDisplayType == 'gene'",
                               selectInput("featureSelect", "Select Feature", 
                                           choices = NULL, # Will be populated dynamically
                                           width = "100%")
                             )
                      )
                    ),
                    
                    # UMAP plot with spinner
                    div(style = "height: 550px;",
                        plotOutput("clusterUMAP", height = "100%") %>% withSpinner()
                    )
                  )
                )
              ),
              
              # Bottom section - Cell type assignment
              fluidRow(
                box(
                  title = "Cell Type Assignment",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  
                  # Cell type assignment UI (will be rendered dynamically)
                  uiOutput("cellTypeUI"),
                  
                  # Action button for assigning cell types
                  div(style = "margin-top: 20px;",
                      actionButton("assignCellTypes", "Assign Cell Types",
                                   class = "btn-success", 
                                   style = "color: white; width: 20%;")
                  )
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
  
  options(shiny.maxRequestSize = 100 * 1024^3)
  
  
  # Initialize log text reactive values
  log_text <- reactiveVal("")
  lsi1_log_text <- reactiveVal("")
  
  display_ui_message <- function(msg) {
    # Function to display messages in the preprocessing UI console
    # Update the log text
    current <- isolate(log_text())
    log_text(paste0(current, msg, "\n"))
    
    # Send to JavaScript console
    js$appendToConsole(msg)
  }
  
  # 4. Create a new function specifically for LSI messages:
  display_lsi1_message <- function(msg) {
    # Function to display messages in the LSI1 UI console
    # Update the log text (if you want to share the log between consoles)
    current <- isolate(lsi1_log_text())
    lsi1_log_text(paste0(current, msg, "\n"))
    
    # Send to JavaScript console using the LSI1-specific function
    js$appendToLSI1_Console(msg)
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
  
  # For LSI1 tab
  assign("send_lsi1_message", function(msg) {
    display_lsi1_message(msg)
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
    LSI1_done = FALSE,
    broad_cluster_done = FALSE,
    LSI2_done = FALSE,
    specific_cluster_done = FALSE
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
    
    withProgress(message = 'Processing folder...', value = 0, {
      # Initialize a list to store Seurat objects
      seurat_list <- list()
      
      # 1. Process files in the main directory with pattern: name_matrix.mtx.gz, etc.
      main_files <- list.files(folder_path, pattern = "(_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz)$", 
                               full.names = TRUE, recursive = FALSE)
      
      # Find base names for files in the main directory
      main_base_names <- unique(sapply(strsplit(basename(main_files), 
                                                "_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz"), 
                                       `[`, 1))
      
      # Calculate the increment amount for each dataset
      total_datasets <- length(main_base_names)
      incr_amount <- 1 / total_datasets
      
      # Process datasets in the main directory
      for (i in seq_along(main_base_names)) {
        base <- main_base_names[i]
        
        # Update progress bar with dataset name and percentage
        incProgress(incr_amount, detail = paste0("Processing ", base, " (", i, "/", total_datasets, ")"))
        
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
        setProgress(detail = paste("Uploading folder", subdir_name))
        
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
            span(id),
            tags$span(
              class = "pull-right",
              style = "color: #888; font-size: 0.9em;"
            )
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
  output$addedColumnsUI <- renderUI({
    req(values$added_columns, values$seurat_list)
    
    if(length(values$added_columns) == 0) {
      return(tags$p("No columns added yet."))
    }
    
    tagList(
      tags$div(
        class = "added-columns-list",
        lapply(values$added_columns, function(col_name) {
          tags$div(
            class = "added-column-item",
            style = "padding: 8px; margin-bottom: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f9f9f9;",
            tags$span(style = "font-weight: bold;", col_name),
            tags$div(
              style = "margin-top: 5px;",
              actionButton(
                inputId = paste0("editColumn_", col_name),
                label = "Edit",
                class = "btn-xs btn-primary",
                style = "margin-right: 5px;"
              ),
              actionButton(
                inputId = paste0("deleteColumn_", col_name),
                label = "Delete",
                class = "btn-xs btn-danger"
              )
            )
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
  
  # Handle delete button clicks for each column
  observe({
    req(values$added_columns, values$seurat_list)
    
    lapply(values$added_columns, function(col_name) {
      observeEvent(input[[paste0("deleteColumn_", col_name)]], {
        # Confirm deletion
        showModal(modalDialog(
          title = "Confirm deletion",
          paste("Are you sure you want to delete the column", col_name, "from all Seurat objects?"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton(paste0("confirmDelete_", col_name), "Delete")
          )
        ))
      })
      
      # Handle confirmation
      observeEvent(input[[paste0("confirmDelete_", col_name)]], {
        removeModal()
        
        # Remove the column from all Seurat objects
        withProgress(message = paste('Deleting column', col_name, 'from all Seurat objects...'), {
          for (i in seq_along(values$seurat_list)) {
            incProgress(1/length(values$seurat_list), 
                        detail = paste("Processing object", i, "of", length(values$seurat_list)))
            
            if (col_name %in% colnames(values$seurat_list[[i]]@meta.data)) {
              values$seurat_list[[i]]@meta.data[[col_name]] <- NULL
            }
          }
          
          # Remove the column name from the list of added columns
          values$added_columns <- values$added_columns[values$added_columns != col_name]
          
          # Show success message
          showNotification(paste("Column", col_name, "deleted successfully!"), type = "message")
        })
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
      
      send_message("Merging datasets to single seurat object")
      merged_seurat <- merge(x=preprocessing_seurat_list[[1]], 
                             y=preprocessing_seurat_list[2:length(preprocessing_seurat_list)])
      merged_seurat$orig.ident <- as.factor(merged_seurat$orig.ident)
      
      send_message("Join count layers in seurat object")
      merged_seurat <- SeuratObject::JoinLayers(merged_seurat)
      values$seurat <- merged_seurat
      values$preprocess_done <- TRUE
      
      values$seurat_list <- NULL
      gc()
      
      send_message("Preprocessing and merging completed successfully!")
    })
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "LSI_1")
  })

  observeEvent(input$commitLSI_1, {
    req(values$seurat)
    
    # Reset the log text if needed
    log_text("")
    
    # Get the Seurat object
    LSI1_seurat <- values$seurat
    
    tryCatch({
      rawCounts <- Seurat::GetAssayData(object=LSI1_seurat, layer = "counts")
      
      # Verify rawCounts is valid
      if (is.null(rawCounts) || length(rawCounts) == 0) {
        send_lsi1_message("Error: Count matrix is empty or null. Attempting to fix...")
        
        # Try to recover counts from the Seurat object
        if ("RNA" %in% names(LSI1_seurat@assays)) {
          
          if ("counts" %in% names(LSI1_seurat@assays$RNA@layers)) {
            rawCounts <- LSI1_seurat@assays$RNA@layers$counts
            send_lsi1_message("Successfully recovered count matrix from RNA assay.")
          
            } else {
            # If counts layer is missing, try to recreate it from data
            send_lsi1_message("Counts layer not found. Creating count matrix from data layer...")
            data_matrix <- LSI1_seurat@assays$RNA@layers$data
            
            if (!is.null(data_matrix)) {
              # Create counts by un-log-transforming data (approximate)
              rawCounts <- 2^data_matrix - 1
              # Set small values to 0 and round
              rawCounts@x[rawCounts@x < 0.5] <- 0
              rawCounts <- round(rawCounts)
              send_lsi1_message("Created approximate count matrix from data layer.")
              
              # Store it back in the Seurat object
              LSI1_seurat@assays$RNA@layers$counts <- rawCounts
            
              } else {
              send_lsi1_message("ERROR: Cannot recover count data. Attempting normal Seurat workflow...")
              # Fall back to standard Seurat normalization
              LSI1_seurat <- Seurat::NormalizeData(LSI1_seurat)
              rawCounts <- Seurat::GetAssayData(object=LSI1_seurat, layer = "counts")
            }
          }
        } else {
          send_lsi1_message("ERROR: RNA assay not found. Cannot proceed.")
          return()
        }
      }
    })
    
    # Get selected parameters
    selected_blacklist <- input$blacklist_genes
    nVarGenes <- input$nVarGenes
    nPCs <- input$nPCs
    resolution <- as.numeric(unlist(strsplit(gsub(" ", "", input$resolution), ",")))
    covariates <- unlist(strsplit(gsub(" ", "", input$covariates), ","))
    umapNeighbors <- input$umapNeighbors
    umapMinDist <- input$umapMinDist
    umapDistMetric <- input$umapDistMetric

    # Log normalization on sparse matrix
    send_lsi1_message("Log normalization of sparse data matrix")
    rawCounts <- Seurat::GetAssayData(object=LSI1_seurat, layer = "counts")
    
    # Generate blacklist based on user selections
    blacklist.genes <- generate_GeneBlacklist(rawCounts, selected_blacklist, 
                                              function(msg) send_lsi1_message(msg))
    
    # If no options selected, show a notification
    if(length(selected_blacklist) == 0) {
      showNotification("Running without ignoring any genes", type = "warning")
    }
    
    send_lsi1_message("Starting LSI preprocessing...")
    
    # Show selected parameters
    send_lsi1_message(sprintf("Using LSI parameters:"))
    send_lsi1_message(sprintf("  - Number of variable genes: %d", nVarGenes))
    send_lsi1_message(sprintf("  - Number of PCs: %d", nPCs))
    send_lsi1_message(sprintf("  - Resolutions: %s", paste(resolution, collapse = ", ")))
    send_lsi1_message(sprintf("  - Covariates for harmony: %s", paste(covariates, collapse = ", ")))
    send_lsi1_message(sprintf("  - UMAP parameters: neighbors=%d, min_dist=%.2f, metric=%s", 
                              umapNeighbors, umapMinDist, umapDistMetric))
    
    # Run LSI
    lsi_results <- process_LSI(
      seurat_obj = LSI1_seurat,
      rawCounts = rawCounts,
      blacklist.genes = blacklist.genes,
      nVarGenes = nVarGenes,
      resolution = resolution,
      nPCs = nPCs,
      covariates = covariates,
      umapNeighbors = umapNeighbors,
      umapMinDist = umapMinDist,
      umapDistMetric = umapDistMetric,
      send_message = function(msg) send_lsi1_message(msg)
    )
    
    # Update values
    values$seurat <- lsi_results$seurat_obj
    
    send_lsi1_message("LSI round 1 processing complete!")
    
    values$LSI1_done <- TRUE
    
    # Navigate to next tab
    updateTabItems(session, "tabs", "BroadClustering")  
  })
  
  # Cell type annotation logic
  observeEvent(input$findMarkers, {
    req(values$seurat)
    
    withProgress(message = 'Finding marker genes...', {
      # Find markers for each cluster
      all_markers <- FindAllMarkers(values$seurat,
                                    test.use = "wilcox",
                                    only.pos = TRUE, 
                                    min.pct = 0.3, 
                                    logfc.threshold = 0.3)
      
      # Store markers
      values$markers <- all_markers %>%
        group_by(cluster) %>%
        top_n(n = input$maxMarkers, wt = avg_log2FC)
      
      values$markers_done <- TRUE
      
      # Update the cluster selection dropdown
      clusters <- sort(unique(as.character(values$markers$cluster)))
      updateSelectInput(session, "clusterSelect", 
                        choices = clusters, 
                        selected = clusters[1])
      
      # Update the feature selection dropdown with all marker genes
      gene_choices <- sort(unique(as.character(values$markers$gene)))
      updateSelectInput(session, "featureSelect", 
                        choices = gene_choices, 
                        selected = gene_choices[1])
      
      # Show notification
      showNotification("Marker genes found successfully", type = "message")
    })
    
    # Render the marker table for the first cluster initially
    output$markerTable <- renderDT({
      req(values$markers, input$clusterSelect)
      values$markers %>% 
        filter(cluster == input$clusterSelect) %>%
        arrange(desc(avg_log2FC)) %>%
        select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
        rename("Gene" = gene, 
               "Log2 FC" = avg_log2FC, 
               "% in Cluster" = pct.1, 
               "% in Other" = pct.2, 
               "Adj. p-value" = p_val_adj)
    }, options = list(pageLength = 10, scrollY = "350px", scroller = TRUE))
    
    # Create UI for cell type assignment
    output$cellTypeUI <- renderUI({
      req(values$seurat, values$markers_done)
      clusters <- sort(levels(values$seurat$seurat_clusters))
      
      tagList(
        fluidRow(
          column(12, 
                 h4("Assign cell types to clusters", 
                    style = "margin-bottom: 20px;")
          )
        ),
        fluidRow(
          lapply(seq_along(clusters), function(i) {
            cluster <- clusters[i]
            # Get top markers for this cluster for the label hint
            top_markers <- values$markers %>% 
              filter(cluster == !!cluster) %>% 
              top_n(3, avg_log2FC) %>% 
              pull(gene) %>% 
              paste(collapse = ", ")
            
            div(
              style = if(i %% 3 == 0) "clear: both;" else "",
              column(
                width = 4,
                div(
                  style = "border: 1px solid #ddd; border-radius: 5px; padding: 10px; margin-bottom: 15px; background-color: #f9f9f9;",
                  strong(paste("Cluster", cluster)),
                  p(style = "font-size: 0.8em; color: #666;", 
                    paste("Top markers:", top_markers)),
                  textInput(
                    inputId = paste0("label_cluster_", cluster),
                    label = NULL,
                    placeholder = "Enter cell type name",
                    value = ifelse(
                      is.null(values$seurat$cell_type) || 
                        !paste0("Cluster ", cluster) %in% unique(values$seurat$cell_type),
                      "", 
                      unique(values$seurat$cell_type[values$seurat$seurat_clusters == cluster])
                    )
                  )
                )
              )
            )
          })
        )
      )
    })
    
    # Update UMAP visualization
    output$clusterUMAP <- renderPlot({
      req(values$seurat)
      
      if(input$umapDisplayType == "clusters") {
        DimPlot(values$seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
          theme(legend.position = "right") +
          ggtitle("Clusters")
      } else {
        req(input$featureSelect)
        FeaturePlot(values$seurat, features = input$featureSelect, 
                    min.cutoff = "q10", max.cutoff = "q90") +
          scale_color_viridis_c() +
          ggtitle(paste("Expression:", input$featureSelect))
      }
    })
  })
  
  # 2. Updating marker table based on cluster selection
  observeEvent(input$clusterSelect, {
    req(values$markers, values$markers_done)
    
    output$markerTable <- renderDT({
      values$markers %>% 
        filter(cluster == input$clusterSelect) %>%
        arrange(desc(avg_log2FC)) %>%
        select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
        rename("Gene" = gene, 
               "Log2 FC" = avg_log2FC, 
               "% in Cluster" = pct.1, 
               "% in Other" = pct.2, 
               "Adj. p-value" = p_val_adj)
    }, options = list(pageLength = 10, scrollY = "350px", scroller = TRUE))
  })
  
  # 3. Cell type assignment
  observeEvent(input$assignCellTypes, {
    req(values$seurat, values$markers_done)
    
    withProgress(message = 'Assigning cell types...', {
      # Get cluster IDs and user-defined labels
      clusters <- levels(values$seurat$seurat_clusters)
      
      # Check if any labels are empty
      labels <- sapply(clusters, function(cluster) {
        label_value <- input[[paste0("label_cluster_", cluster)]]
        # Use default if empty
        if(is.null(label_value) || trimws(label_value) == "") {
          return(paste0("Cluster ", cluster))
        } else {
          return(label_value)
        }
      })
      
      # Create new identity column
      new_ids <- plyr::mapvalues(values$seurat$seurat_clusters, 
                                 from = clusters, 
                                 to = labels)
      values$seurat$cell_type <- new_ids
      
      # Add metadata about annotation
      values$seurat@misc$annotation_info <- list(
        date = Sys.time(),
        clusters = data.frame(
          cluster = clusters,
          label = labels,
          stringsAsFactors = FALSE
        )
      )
      
      values$annotation_done <- TRUE
      
      # Show success notification
      showNotification("Cell types assigned successfully", type = "message")
    })
    
    # Update UMAP plot to show cell types
    output$clusterUMAP <- renderPlot({
      req(values$seurat, values$annotation_done)
      
      # Add option to view by cell type in the display type dropdown
      updateSelectInput(session, "umapDisplayType",
                        choices = c("Clusters" = "clusters", 
                                    "Cell Types" = "celltypes",
                                    "Gene Expression" = "gene"),
                        selected = "celltypes")
      
      # Show the cell type UMAP
      DimPlot(values$seurat, reduction = "umap", group.by = "cell_type", label = TRUE) +
        theme(legend.position = "right") +
        ggtitle("Assigned Cell Types")
    })
    
    # Navigate to next tab (which will be LSI_2 now)
    updateTabItems(session, "tabs", "LSI_2")
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
shinyApp(ui = ui, server = server, options = list(
  width = 1200,
  height = 1200
))