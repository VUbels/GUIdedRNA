message("Starting GUIdedRNA app...")

# Check if running in package context
if (!requireNamespace("GUIdedRNA", quietly = TRUE)) {
  stop("GUIdedRNA package not found. Please install the package first.")
}

# Source global.R which is in the same directory
if (file.exists("global.R")) {
  message("Found global.R, loading it...")
  source("global.R")
} else {
  message("global.R not found in working directory:", getwd())
}


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
      
      menuItem("Integrate Results", tabName = "integrate_results", icon = icon("puzzle-piece")),
      
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
  };
      
  shinyjs.appendToLSI2_Console = function(message) {
    var consoleOutput = document.getElementById('lsi2_ConsoleOutput');
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
      functions = c("appendToConsole", "appendToLSI1_Console", "appendToLSI2_Console")
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
                  status = "primary",
                  fluidRow(
                    # Column 1: Display original identities from all objects
                    column(
                      width = 4,
                      status = "primary",
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
                      status = "primary",
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
                      style = "white-space: pre-wrap; height: 575px; overflow-y: auto; background-color: #f5f5f5; padding: 1px; font-family: monospace; text-align: left; vertical-align: top;")
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
                    checkboxGroupInput("blacklist_genes_1", "Ignore genes for LSI clustering (Not DE analysis)", 
                                       choices = c("Ignore Mitochondrial Genes" = "blacklist_mitogenes", 
                                                   "Ignore X/Y chomosomal genes" = "blacklist_sexgenes",
                                                   "Ignore Ribosomal genes" = "blacklist_rbgenes"),
                                       selected = c("blacklist_mitogenes", "blacklist_sexgenes", "blacklist_rbgenes")
                                       
                    )
                  ),
                  
                  box(
                    title = "LSI Round 1 Parameters",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    
                    numericInput("nVarGenes_1", "Number of Variable Genes", 
                                 value = 2000, min = 100, max = 5000, step = 100),
                    
                    sliderInput("nPCs_1", "Number of Principal Components", 
                                value = 25, min = 5, max = 50, step = 5),
                    
                    textInput("resolution_1", "Clustering Resolutions (comma separated)", 
                              value = "0.2, 0.4, 0.6"),
                    
                    textInput("covariates_1", "Harmony Covariates", 
                              value = "orig.ident"),
                    
                    # UMAP parameters
                    numericInput("umapNeighbors_1", "UMAP Neighbors", 
                                 value = 50, min = 5, max = 100, step = 5),
                    
                    numericInput("umapMinDist_1", "UMAP Minimum Distance", 
                                 value = 0.5, min = 0.01, max = 0.99, step = 0.01),
                    
                    selectInput("umapDistMetric_1", "UMAP Distance Metric",
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
                                min = 5, max = 30, value = 15, step = 5),
                    
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
                    height = "580px",
                    status = "primary",
                    solidHeader = TRUE,
                    
                    # Cluster selection dropdown
                    selectInput("clusterSelect", "Select Cluster", 
                                choices = NULL, # Will be populated dynamically
                                width = "100%"),
                    
                    # Marker gene table with scroll - simplified to show only genes and stats
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
                    height = "810px", # Match left column height
                    status = "primary",
                    solidHeader = TRUE,
                    
                    # Feature selection for UMAP coloring
                    fluidRow(
                      column(4, 
                             selectInput("umapDisplayType", "Display Type",
                                         choices = c("Clusters" = "clusters", 
                                                     "Gene Expression" = "gene"),
                                         selected = "clusters",
                                         width = "100%")
                      ),
                      column(4,
                             # Color By selection - always shown
                             selectInput("umapGroupBy", "Color By", 
                                         choices = c("seurat_clusters"), 
                                         selected = "seurat_clusters",
                                         width = "100%")
                      ),
                      column(4,
                             # Split By selection - always shown
                             selectInput("umapSplitBy", "Split By", 
                                         choices = c("None" = "none"), 
                                         selected = "none",
                                         width = "100%")
                      )
                    ),
                    
                    # Conditional UI for gene selection - shown in its own row when gene expression is selected
                    conditionalPanel(
                      condition = "input.umapDisplayType == 'gene'",
                      fluidRow(
                        column(12,
                               selectizeInput("featureSelect", "Select Gene", 
                                              choices = NULL,
                                              options = list(
                                                placeholder = 'Type to search genes',
                                                onInitialize = I('function() { this.setValue(""); }')
                                              ),
                                              width = "100%")
                        )
                      )
                    ),
                    
                    # UMAP plot with spinner - fixed aspect ratio and contained in its own div
                    div(style = "height: 550px; width: 100%; margin-top: 10px; overflow: hidden;",
                        plotOutput("clusterUMAP", height = "550px", width = "100%") %>% withSpinner()
                    ),
                    
                    # Download button positioned at bottom right
                    div(style = "position: absolute; bottom: 10px; right: 15px; z-index: 1000;",
                        actionButton("downloadUMAP_initial", "Download Figure",
                                     class = "btn-info btn-sm",
                                     style = "color: white; font-size: 12px;",
                                     icon = icon("download"))
                    )
                  )
                )
              ),
              # Cell type assignment UI
              fluidRow(
                column(12,
                       box(
                         title = "Cell Type Assignment",
                         width = NULL,
                         status = "primary",
                         solidHeader = TRUE,
                         uiOutput("cellTypeUI") %>% withSpinner()
                       )
                )
              )
      ),
      
      tabItem(tabName = "LSI_2",
              fluidRow(
                column(
                  width = 4,
                  height = "900px",
                  box(
                    title = "Blacklist Settings",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    checkboxGroupInput("blacklist_genes_2", "Ignore genes for LSI clustering (Not DE analysis)", 
                                       choices = c("Ignore Mitochondrial Genes" = "blacklist_mitogenes", 
                                                   "Ignore X/Y chomosomal genes" = "blacklist_sexgenes",
                                                   "Ignore Ribosomal genes" = "blacklist_rbgenes"),
                                       selected = c("blacklist_mitogenes", "blacklist_sexgenes", "blacklist_rbgenes")
                                       
                    )
                  ),
                  
                  box(
                    title = "LSI Round 2 Parameters",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    
                    numericInput("nVarGenes_2", "Number of Variable Genes", 
                                 value = 2000, min = 100, max = 5000, step = 100),
                    
                    sliderInput("nPCs_2", "Number of Principal Components", 
                                value = 15, min = 5, max = 50, step = 5),
                    
                    textInput("resolution_2", "Clustering Resolutions (comma separated)", 
                              value = "0.2, 0.3"),
                    
                    textInput("covariates_2", "Harmony Covariates", 
                              value = "orig.ident"),
                    
                    # UMAP parameters
                    numericInput("umapNeighbors_2", "UMAP Neighbors", 
                                 value = 35, min = 5, max = 500, step = 5),
                    
                    numericInput("umapMinDist_2", "UMAP Minimum Distance", 
                                 value = 0.4, min = 0.01, max = 0.99, step = 0.01),
                    
                    selectInput("umapDistMetric_2", "UMAP Distance Metric",
                                choices = c("cosine", "euclidean", "manhattan", "hamming"),
                                selected = "cosine"),
                    
                    actionButton("commitLSI_2", "Run LSI - Round 2", 
                                 class = "btn-success", 
                                 style = "color: white; width: 100%;")
                  )
                ),
                
                box(
                  title = "LSI Round 2 Results",
                  width = 8,
                  height = "930px",
                  div(id = "lsi2_ConsoleOutput", 
                      style = "white-space: pre-wrap; height: 850px; overflow-y: auto; background-color: #f5f5f5; padding: 1px; font-family: monospace; text-align: left; vertical-align: top;")
                )
              )
      ),
      
      # Cell Type Annotation tab
      tabItem(tabName = "annotate_specific",
              fluidRow(
                # Left column - Subset selector and marker analysis
                column(
                  width = 4,
                  # New subset selector box
                  box(
                    title = "Subset Selection",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    selectInput("subsetSelector", "Select Cell Type Subset", 
                                choices = NULL,
                                width = "100%"),
                    helpText("Switch between different cell type subsets processed in LSI Round 2")
                  ),
                  
                  box(
                    title = "Marker Gene Analysis",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    
                    sliderInput("maxMarkers_final", "Max Markers per Cluster", 
                                min = 5, max = 30, value = 15, step = 5),
                    
                    div(style = "margin-top: 15px; margin-bottom: 15px;",
                        actionButton("findAllMarkersAtOnce", "Calculate All Subset Markers",
                                     class = "btn-success", 
                                     style = "color: white; width: 100%;")
                    )
                  ),
                  
                  box(
                    title = "Differentially Expressed Genes",
                    width = NULL,
                    height = "520px",
                    status = "primary",
                    solidHeader = TRUE,
                    
                    selectInput("clusterSelect_final", "Select Cluster", 
                                choices = NULL,
                                width = "100%"),
                    
                    div(style = "height: 340px; overflow-y: auto;",
                        DTOutput("markerTable_final") %>% withSpinner()
                    )
                  )
                ),
                
                # Right column - UMAP visualization
                column(
                  width = 8,
                  box(
                    title = "Cluster Visualization",
                    width = NULL,
                    height = "810px",
                    status = "primary",
                    solidHeader = TRUE,
                    
                    # Feature selection for UMAP coloring
                    fluidRow(
                      column(4, 
                             selectInput("umapDisplayType_final", "Display Type",
                                         choices = c("Clusters" = "clusters", 
                                                     "Gene Expression" = "gene"),
                                         selected = "clusters",
                                         width = "100%")
                      ),
                      column(4,
                             selectInput("umapGroupBy_final", "Color By", 
                                         choices = c("seurat_clusters"), 
                                         selected = "seurat_clusters",
                                         width = "100%")
                      ),
                      column(4,
                             selectInput("umapSplitBy_final", "Split By", 
                                         choices = c("None" = "none"), 
                                         selected = "none",
                                         width = "100%")
                      )
                    ),
                    
                    # Conditional UI for gene selection
                    conditionalPanel(
                      condition = "input.umapDisplayType_final == 'gene'",
                      fluidRow(
                        column(12,
                               selectizeInput("featureSelect_final", "Select Gene", 
                                              choices = NULL,
                                              options = list(
                                                placeholder = 'Type to search genes',
                                                onInitialize = I('function() { this.setValue(""); }')
                                              ),
                                              width = "100%")
                        )
                      )
                    ),
                    
                    # UMAP plot
                    div(style = "height: 550px; width: 100%; margin-top: 10px; overflow: hidden;",
                        plotOutput("clusterUMAP_final", height = "550px", width = "100%") %>% withSpinner()
                    ),
                    
                    # Download button positioned at bottom right
                    div(style = "position: absolute; bottom: 10px; right: 15px; z-index: 1000;",
                        actionButton("downloadUMAP_final", "Download Figure",
                                     class = "btn-info btn-sm",
                                     style = "color: white; font-size: 12px;",
                                     icon = icon("download"))
                    )
                  )
                )
              ),
              # Cell type assignment UI for final clustering
              fluidRow(
                column(12,
                       box(
                         title = "Final Cell Type Assignment",
                         width = NULL,
                         status = "primary",
                         solidHeader = TRUE,
                         uiOutput("cellTypeUI_final") %>% withSpinner()
                       )
                )
              )
      ),
      
      # New Integration Results Tab
      tabItem(tabName = "integrate_results",
              fluidRow(
                box(
                  title = "Integration Controls",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  
                  fluidRow(
                    column(12,
                           h4("Integrate Final Cell Type Annotations"),
                           p("This will combine all final cell type annotations from the individual subsets 
                       back into the original Seurat object as a new metadata column."),
                           actionButton("integrateResults", "Integrate Final Annotations",
                                        class = "btn-success", 
                                        style = "color: white; width: 100%;")
                        )
                      )
                    )
              ),
              
              fluidRow(
                # Integration Summary
                column(6,
                       box(
                         title = "Integration Summary",
                         width = NULL,
                         status = "primary",
                         solidHeader = TRUE,
                         height = "600px",
                         verbatimTextOutput("integrationSummary")
                       )
                ),
                
                # Final UMAP visualization
                column(6,
                       box(
                         title = "Final Integrated Results",
                         width = NULL,
                         status = "primary", 
                         solidHeader = TRUE,
                         height = "600px",
                         plotOutput("finalIntegratedUMAP", height = "550px") %>% withSpinner()
                       )
                )
              )
      ),
      
      
      
      
      # Download Results tab
      tabItem(tabName = "download",
              fluidRow(
                box(
                  title = "Final Integrated Results",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  
                  checkboxGroupInput("finalExportItems", "Select Final Results to Export:",
                                     choices = c("Processed Annotated Integrated Seurat Object (RDS)" = "final_rds",
                                                 "Processed Count Matrix (MTX)" = "count_matrix",
                                                 "Cell Type Metadata (CSV)" = "final_meta",
                                                 "Subset-specific Marker Genes (CSV)" = "subset_markers",
                                                 "Dimensionality Reduction Coordinates (CSV)" = "dimred_coords",
                                                 "Integration Summary Report (HTML)" = "integration_report")),
                  
                  actionButton("downloadIntegratedData", "Save Results to Output Directory",
                               class = "btn-primary",
                               style = "color: white; width: 100%;")
                  
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
  lsi2_log_text <- reactiveVal("")
  
  #Preprocessing Message Logger
  display_ui_message <- function(msg) {
    # Function to display messages in the preprocessing UI console
    # Update the log text
    current <- isolate(log_text())
    log_text(paste0(current, msg, "\n"))
    
    # Send to JavaScript console
    js$appendToConsole(msg)
  }
  
  #LSI_1 Message Logger
  display_lsi1_message <- function(msg) {
    # Function to display messages in the LSI1 UI console
    # Update the log text (if you want to share the log between consoles)
    current <- isolate(lsi1_log_text())
    lsi1_log_text(paste0(current, msg, "\n"))
    
    # Send to JavaScript console using the LSI1-specific function
    js$appendToLSI1_Console(msg)
  }
  
  #LSI_2 Message Logger
  display_lsi2_message <- function(msg) {
    # Function to display messages in the LSI1 UI console
    # Update the log text (if you want to share the log between consoles)
    current <- isolate(lsi2_log_text())
    lsi2_log_text(paste0(current, msg, "\n"))
    
    # Send to JavaScript console using the LSI1-specific function
    js$appendToLSI2_Console(msg)
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
  
  assign("send_lsi2_message", function(msg) {
    display_lsi2_message(msg)
  }, envir = .GlobalEnv)
  
  # Render log text
  output$preprocessingLog <- renderText({
    log_text()
  })
  
  values <- reactiveValues(
    seurat = NULL,
    original_seurat = NULL,
    seurat_list = NULL,  
    original_seurat_list = NULL,  
    qc_done = FALSE,
    preprocess_done = FALSE,
    LSI1_done = FALSE,
    broad_cluster_done = FALSE,
    LSI2_done = FALSE,
    specific_cluster_done = FALSE,
    subset_objects = NULL,
    subset_markers = NULL,
    integration_summary = NULL,
    markers = NULL,  
    added_columns = NULL
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
    selected_blacklist <- input$blacklist_genes_1
    nVarGenes <- input$nVarGenes_1
    nPCs <- input$nPCs_1
    resolution <- as.numeric(unlist(strsplit(gsub(" ", "", input$resolution_1), ",")))
    covariates <- unlist(strsplit(gsub(" ", "", input$covariates_1), ","))
    umapNeighbors <- input$umapNeighbors_1
    umapMinDist <- input$umapMinDist_1
    umapDistMetric <- input$umapDistMetric_1

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
    updateTabItems(session, "tabs", "annotate_broad")  
  })
  
  
  # Render the marker table for differentially expressed genes - completely stripped down
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
        dplyr::top_n(n = input$maxMarkers, wt = avg_log2FC)
      
      values$markers_done <- TRUE
      
      # Update the cluster selection dropdown
      clusters <- sort(unique(as.character(values$markers$cluster)))
      updateSelectInput(session, "clusterSelect", 
                        choices = clusters, 
                        selected = clusters[1])
      
      # Initialize cell_type column if it doesn't exist
      if(!"cell_type" %in% colnames(values$seurat@meta.data)) {
        values$seurat$cell_type <- paste0("Cluster ", values$seurat$seurat_clusters)
      }
      
      # Show notification
      showNotification("Marker genes found successfully", type = "message")
    })
  })

  # Observer to initialize UI components when Seurat object is loaded
  observeEvent(values$seurat, {
    req(values$seurat)
    
    # Get all features for gene selection - INDEPENDENT of marker analysis
    all_features <- rownames(values$seurat)
    
    # Update feature selection with ALL genes from the dataset
    updateSelectizeInput(session, "featureSelect",
                         choices = all_features,
                         server = TRUE)
    
    # Get all metadata columns
    meta_columns <- colnames(values$seurat@meta.data)
    relevant_columns <- c("seurat_clusters", "orig.ident")
    
    # Add cell_type if it exists
    if("cell_type" %in% meta_columns) {
      relevant_columns <- c(relevant_columns, "cell_type")
    }
    
    # Add any custom columns from sample information tab
    if(!is.null(values$added_columns)) {
      relevant_columns <- c(relevant_columns, values$added_columns)
    }
    
    # Filter to only include columns that exist
    relevant_columns <- relevant_columns[relevant_columns %in% meta_columns]
    
    # Update the Color By dropdown
    updateSelectInput(session, "umapGroupBy", 
                      choices = relevant_columns, 
                      selected = "seurat_clusters")
    
    # Create split by options with "None" option first
    split_options <- c("None" = "none")
    
    # Add orig.ident if it exists
    if("orig.ident" %in% meta_columns) {
      split_options <- c(split_options, "orig.ident" = "orig.ident")
    }
    
    # Add cell_type if it exists
    if("cell_type" %in% meta_columns) {
      split_options <- c(split_options, "cell_type" = "cell_type")
    }
    
    # Add any custom columns for split.by
    if(!is.null(values$added_columns)) {
      for(col in values$added_columns) {
        if(col %in% meta_columns) {
          split_options[col] <- col
        }
      }
    }
    
    # Update the Split By dropdown
    updateSelectInput(session, "umapSplitBy", 
                      choices = split_options, 
                      selected = "none")
  })
  
  #
  output$markerTable <- renderDT({
    req(values$markers, input$clusterSelect)
    
    # Create a fresh dataframe with only the three needed columns
    marker_genes <- values$markers %>% 
      filter(cluster == input$clusterSelect) %>%
      arrange(desc(avg_log2FC))
    
    # Create a new dataframe with properly formatted columns
    display_table <- data.frame(
      Gene = marker_genes$gene,
      "Log2 FC" = round(marker_genes$avg_log2FC, 3),
      "Adj. p-value" = sprintf("%.2f", marker_genes$p_val_adj),
      stringsAsFactors = FALSE
    )
    
    # Use the DT package explicitly to have more control
    DT::datatable(
      display_table,
      options = list(
        dom = 't',  # table only, no search/pagination UI
        scrollY = "350px",
        scroller = TRUE,
        paging = FALSE,
        ordering = TRUE,
        columnDefs = list(
          list(targets = 0, className = "dt-left"),  # Gene column left aligned
          list(targets = 1:2, className = "dt-right") # Numeric columns right aligned
        )
      ),
      rownames = FALSE,  # No row names/indices
      selection = "none", # No row selection
      class = "compact stripe" # Styling
    )
  })
  
  observeEvent(input$applyLabels, {
    req(values$seurat)
    
    withProgress(message = 'Applying cell type labels...', {
      # Get current clusters
      clusters <- sort(levels(Idents(values$seurat)))
      
      # Initialize cell_type column if it doesn't exist
      if(!"cell_type" %in% colnames(values$seurat@meta.data)) {
        values$seurat$cell_type <- paste0("Cluster ", values$seurat$seurat_clusters)
      }
      
      # Apply new cell type labels
      for(cluster in clusters) {
        input_id <- paste0("label_cluster_", cluster)
        label_value <- input[[input_id]]
        
        # Only update if not empty
        if(!is.null(label_value) && nchar(trimws(label_value)) > 0) {
          values$seurat$cell_type[values$seurat$seurat_clusters == cluster] <- label_value
        }
      }
      
      # Update the group.by dropdown to include cell_type
      meta_columns <- colnames(values$seurat@meta.data)
      relevant_columns <- c("seurat_clusters", "orig.ident", "cell_type")
      if(!is.null(values$added_columns)) {
        relevant_columns <- c(relevant_columns, values$added_columns)
      }
      relevant_columns <- relevant_columns[relevant_columns %in% meta_columns]
      
      # Update both dropdowns
      updateSelectInput(session, "umapGroupBy", 
                        choices = relevant_columns, 
                        selected = "cell_type")
      
      # Update split by options to include cell_type
      split_options <- c("None" = "none")
      for(col in relevant_columns) {
        if(col != "seurat_clusters") {  # Don't include seurat_clusters in split options
          split_options[col] <- col
        }
      }
      
      updateSelectInput(session, "umapSplitBy", 
                        choices = split_options, 
                        selected = "none")
      
      # Show notification
      showNotification("Cell type labels applied successfully", type = "message")
    })
  }) 
  
  output$cellTypeUI <- renderUI({
    req(values$seurat)
    clusters <- sort(as.numeric(as.character(levels(Idents(values$seurat)))))
    
    tagList(
      fluidRow(
        column(12, 
               h4("Assign cell types to clusters", 
                  style = "margin-bottom: 20px;")
        )
      ),
      fluidRow(
        column(12,
               div(
                 style = "max-height: 400px; overflow-y: auto; border: 1px solid #ddd; border-radius: 5px; padding: 15px; background-color: #f9f9f9;",
                 lapply(clusters, function(cluster) {
                   # Get top markers for this cluster for the label hint
                   top_markers <- ""
                   if(!is.null(values$markers)) {
                     top_markers <- values$markers %>% 
                       filter(cluster == !!cluster) %>% 
                       dplyr::top_n(3, avg_log2FC) %>% 
                       pull(gene) %>% 
                       paste(collapse = ", ")
                   }
                   
                   div(
                     style = "margin-bottom: 8px; display: flex; align-items: center;",
                     div(
                       style = "flex: 0 0 150px; font-weight: bold;",
                       paste("Cluster", cluster)
                     ),
                     div(
                       style = "flex: 0 0 250px; font-size: 0.8em; color: #666;",
                       if(top_markers != "") paste("Top markers:", top_markers) else ""
                     ),
                     div(
                       style = "flex: 1;",
                       textInput(
                         inputId = paste0("label_cluster_", cluster),
                         label = NULL,
                         placeholder = "Enter cell type name",
                         value = ifelse(
                           is.null(values$seurat$cell_type) || 
                             !any(values$seurat$seurat_clusters == cluster),
                           "", 
                           unique(values$seurat$cell_type[values$seurat$seurat_clusters == cluster])
                         ),
                         width = "100%"
                       )
                     )
                   )
                 })
               )
        )
      ),
      # Add apply button to save cell type assignments
      fluidRow(
        column(12,
               div(
                 style = "margin-top: 20px;",
                 actionButton("applyLabels", "Apply Cell Type Labels",
                              class = "btn-success", 
                              style = "color: white; width: 100%;")
               )
        )
      )
    )
  })
  
  # Update UMAP visualization - this is the main fix
  output$clusterUMAP <- renderPlot({
    req(values$seurat)
    
    # Check if UMAP reduction exists
    if(!"umap" %in% names(values$seurat@reductions)) {
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "UMAP not computed yet. Please run LSI first.") +
          theme_void() +
          theme(
            plot.title = element_text(size = 16, hjust = 0.5),
            aspect.ratio = 1
          )
      )
    }
    
    # Define common theme elements for consistent plotting
    common_theme <- theme(
      legend.position = "right",
      aspect.ratio = 1,
      plot.margin = margin(10, 10, 10, 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
    
    if(input$umapDisplayType == "clusters") {
      # For cluster visualization
      if(is.null(input$umapSplitBy) || input$umapSplitBy == "none") {
        # No splitting
        p <- DimPlot(values$seurat, 
                     reduction = "umap", 
                     group.by = input$umapGroupBy, 
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste("Colored by:", input$umapGroupBy))
      } else {
        # With split.by
        p <- DimPlot(values$seurat, 
                     reduction = "umap", 
                     group.by = input$umapGroupBy, 
                     split.by = input$umapSplitBy,
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste("Colored by:", input$umapGroupBy, "| Split by:", input$umapSplitBy))
      }
    } else {
      # For gene expression
      req(input$featureSelect)
      if(is.null(input$featureSelect) || input$featureSelect == "" || 
         !input$featureSelect %in% rownames(values$seurat)) {
        # Return an error plot if feature not found
        p <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Selected gene not found in dataset") +
          theme_void() +
          common_theme
      } else if(is.null(input$umapSplitBy) || input$umapSplitBy == "none") {
        # Gene expression without splitting
        p <- FeaturePlot(values$seurat, 
                         features = input$featureSelect, 
                         min.cutoff = "q10", 
                         max.cutoff = "q90",
                         pt.size = 0.8) +
          scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                                name = "Expression") +
          common_theme +
          ggtitle(paste("Expression:", input$featureSelect))
      } else {
        # Gene expression with splitting
        tryCatch({
          # Get the split factor
          split_factor <- values$seurat@meta.data[[input$umapSplitBy]]
          unique_splits <- unique(split_factor)
          
          # Create individual plots for each split
          plot_list <- lapply(unique_splits, function(split_val) {
            # Get cells for this split
            cells_subset <- colnames(values$seurat)[split_factor == split_val]
            
            FeaturePlot(values$seurat, 
                        features = input$featureSelect,
                        cells = cells_subset,
                        min.cutoff = "q10", 
                        max.cutoff = "q90",
                        pt.size = 0.8) +
              scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                                   name = "Expression") +
              common_theme +
              ggtitle(paste(split_val))
          })
          
          # Combine the plots
          p <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
          
          # Add an overall title
          p <- cowplot::plot_grid(
            cowplot::ggdraw() + 
              cowplot::draw_label(paste("Expression:", input$featureSelect, "| Split by:", input$umapSplitBy), 
                                  fontface = 'bold', size = 14),
            p,
            ncol = 1,
            rel_heights = c(0.1, 1)
          )
        }, error = function(e) {
          p <- ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = paste("Error creating split plot:", e$message)) +
            theme_void() +
            common_theme
        })
      }
    }
    
    # Return the plot
    return(p)
  }, height = 550, width = 700) # Fixed dimensions to maintain aspect ratio
  
  # Add this observer to update dropdowns when markers are found
  observeEvent(values$markers, {
    req(values$markers)
    clusters <- sort(as.numeric(as.character(unique(values$markers$cluster))))
    updateSelectInput(session, "clusterSelect", 
                      choices = clusters, 
                      selected = clusters[1])
  })
  
  
  # LSI Round 2 Server Logic
  observeEvent(input$commitLSI_2, {
    req(values$seurat)
    
    # Check if cell type assignments exist
    if(!"cell_type" %in% colnames(values$seurat@meta.data)) {
      showNotification("Please complete initial clustering and cell type assignment first", type = "error")
      return()
    }
    
    # Reset the log text
    lsi2_log_text("")
    
    # Get unique cell types (excluding any that might be "Cluster X" format)
    cell_types <- unique(values$seurat$cell_type)
    cell_types <- cell_types[!grepl("^Cluster \\d+$", cell_types)]
    
    if(length(cell_types) == 0) {
      showNotification("No custom cell types found. Please assign cell types in Initial Clustering first.", type = "error")
      return()
    }
    
    send_lsi2_message(sprintf("Starting LSI Round 2 for %d cell types...", length(cell_types)))
    
    # Get selected parameters for LSI Round 2
    selected_blacklist_2 <- input$blacklist_genes_2
    nVarGenes_2 <- input$nVarGenes_2
    nPCs_2 <- input$nPCs_2
    resolution_2 <- as.numeric(unlist(strsplit(gsub(" ", "", input$resolution_2), ",")))
    covariates_2 <- unlist(strsplit(gsub(" ", "", input$covariates_2), ","))
    umapNeighbors_2 <- input$umapNeighbors_2
    umapMinDist_2 <- input$umapMinDist_2
    umapDistMetric_2 <- input$umapDistMetric_2
    
    # Initialize list to store subset objects
    subset_objects <- list()
    
    withProgress(message = 'Processing cell type subsets...', value = 0, {
      
      for(i in seq_along(cell_types)) {
        cell_type <- cell_types[i]
        
        incProgress(1/length(cell_types), detail = paste("Processing", cell_type))
        send_lsi2_message(sprintf("Processing subset %d/%d: %s", i, length(cell_types), cell_type))
        
        # Create subset for this cell type
        subset_cells <- colnames(values$seurat)[values$seurat$cell_type == cell_type]
        
        if(length(subset_cells) < 50) {
          send_lsi2_message(sprintf("  - Skipping %s: insufficient cells (%d)", cell_type, length(subset_cells)))
          next
        }
        
        send_lsi2_message(sprintf("  - Creating subset with %d cells", length(subset_cells)))
        seurat_subset <- subset(values$seurat, cells = subset_cells)
        
        # Get raw counts for this subset
        tryCatch({
          rawCounts_subset <- Seurat::GetAssayData(object = seurat_subset, layer = "counts")
          
          # Generate blacklist for this subset
          blacklist.genes_subset <- generate_GeneBlacklist(rawCounts_subset, selected_blacklist_2, 
                                                           function(msg) send_lsi2_message(paste("  -", msg)))
          
          send_lsi2_message(sprintf("  - Running LSI on %s subset...", cell_type))
          
          # Run LSI on subset
          lsi_results_subset <- process_LSI(
            seurat_obj = seurat_subset,
            rawCounts = rawCounts_subset,
            blacklist.genes = blacklist.genes_subset,
            nVarGenes = nVarGenes_2,
            resolution = resolution_2,
            nPCs = nPCs_2,
            covariates = covariates_2,
            umapNeighbors = umapNeighbors_2,
            umapMinDist = umapMinDist_2,
            umapDistMetric = umapDistMetric_2,
            send_message = function(msg) send_lsi2_message(paste("    -", msg))
          )
          
          # Store the processed subset
          subset_objects[[cell_type]] <- lsi_results_subset$seurat_obj
          
          send_lsi2_message(sprintf("  - Completed LSI for %s", cell_type))
          
        }, error = function(e) {
          send_lsi2_message(sprintf("  - Error processing %s: %s", cell_type, e$message))
        })
      }
    })
    
    # Store subset objects in reactive values
    values$subset_objects <- subset_objects
    values$LSI2_done <- TRUE
    
    send_lsi2_message(sprintf("LSI Round 2 complete! Processed %d subsets.", length(subset_objects)))
    
    # Navigate to final clustering tab
    updateTabItems(session, "tabs", "annotate_specific")
  })
  
  # Final Clustering Logic - Add subset selector observer
  observeEvent(values$subset_objects, {
    req(values$subset_objects)
    
    subset_names <- names(values$subset_objects)
    updateSelectInput(session, "subsetSelector", 
                      choices = subset_names, 
                      selected = subset_names[1])
  })
  
  # Reactive to get current subset
  current_subset <- reactive({
    req(values$subset_objects, input$subsetSelector)
    values$subset_objects[[input$subsetSelector]]
  })
  
  # Observer to update UI components when subset changes
  observeEvent(input$subsetSelector, {
    req(current_subset())
    
    # Get all features for gene selection from current subset
    all_features <- rownames(current_subset())
    
    # Update feature selection with genes from current subset
    updateSelectizeInput(session, "featureSelect_final",
                         choices = all_features,
                         server = TRUE)
    
    # Get metadata columns from current subset
    meta_columns <- colnames(current_subset()@meta.data)
    relevant_columns <- c("seurat_clusters", "orig.ident")
    
    # Add cell_type if it exists
    if("cell_type" %in% meta_columns) {
      relevant_columns <- c(relevant_columns, "cell_type")
    }
    
    # Add any custom columns
    if(!is.null(values$added_columns)) {
      relevant_columns <- c(relevant_columns, values$added_columns)
    }
    
    # Filter to only include columns that exist
    relevant_columns <- relevant_columns[relevant_columns %in% meta_columns]
    
    # Update dropdowns for final clustering
    updateSelectInput(session, "umapGroupBy_final", 
                      choices = relevant_columns, 
                      selected = "seurat_clusters")
    
    # Create split by options
    split_options <- c("None" = "none")
    for(col in relevant_columns) {
      if(col != "seurat_clusters") {
        split_options[col] <- col
      }
    }
    
    updateSelectInput(session, "umapSplitBy_final", 
                      choices = split_options, 
                      selected = "none")
    
    # Update cluster selection dropdown if markers exist for this subset
    if(!is.null(values$subset_markers) && !is.null(values$subset_markers[[input$subsetSelector]])) {
      clusters <- sort(unique(as.character(values$subset_markers[[input$subsetSelector]]$cluster)))
      updateSelectInput(session, "clusterSelect_final", 
                        choices = clusters, 
                        selected = clusters[1])
    } else {
      # Clear cluster selection if no markers calculated yet
      updateSelectInput(session, "clusterSelect_final", 
                        choices = character(0), 
                        selected = NULL)
    }
  })
  
  # Find markers for final clustering
  observeEvent(input$findAllMarkersAtOnce, {
    req(values$subset_objects)
    
    # Check if we have subset objects
    if(length(values$subset_objects) == 0) {
      showNotification("No subset objects found. Please complete LSI Round 2 first.", type = "error")
      return()
    }
    
    # Initialize the subset_markers list
    values$subset_markers <- list()
    
    withProgress(message = 'Calculating markers for all subsets...', value = 0, {
      
      subset_names <- names(values$subset_objects)
      
      for(i in seq_along(subset_names)) {
        subset_name <- subset_names[i]
        subset_obj <- values$subset_objects[[subset_name]]
        
        incProgress(1/length(subset_names), 
                    detail = paste("Finding markers for", subset_name, 
                                   paste0("(", i, "/", length(subset_names), ")")))
        
        # Check if subset has enough cells and clusters
        n_cells <- ncol(subset_obj)
        n_clusters <- length(unique(subset_obj$seurat_clusters))
        
        if(n_cells < 50) {
          showNotification(paste("Skipping", subset_name, ": insufficient cells (", n_cells, ")"), 
                           type = "warning")
          next
        }
        
        if(n_clusters < 2) {
          showNotification(paste("Skipping", subset_name, ": insufficient clusters (", n_clusters, ")"), 
                           type = "warning")
          next
        }
        
        tryCatch({
          # Find all markers for this subset
          all_markers <- FindAllMarkers(subset_obj,
                                        test.use = "wilcox",
                                        only.pos = TRUE, 
                                        min.pct = 0.3, 
                                        logfc.threshold = 0.3,
                                        verbose = FALSE)  # Suppress verbose output
          
          # Store top markers per cluster
          values$subset_markers[[subset_name]] <- all_markers %>%
            group_by(cluster) %>%
            dplyr::top_n(n = input$maxMarkers_final, wt = avg_log2FC) %>%
            arrange(cluster, desc(avg_log2FC))
          
          # Initialize cell_type_final column if it doesn't exist
          if(!"cell_type_final" %in% colnames(subset_obj@meta.data)) {
            values$subset_objects[[subset_name]]$cell_type_final <- 
              paste0("Cluster ", subset_obj$seurat_clusters)
          }
          
        }, error = function(e) {
          showNotification(paste("Error finding markers for", subset_name, ":", e$message), 
                           type = "warning")
        })
      }
      
      # Update the cluster selection dropdown for the current subset
      if(length(values$subset_markers) > 0 && !is.null(input$subsetSelector)) {
        current_subset_name <- input$subsetSelector
        if(!is.null(values$subset_markers[[current_subset_name]])) {
          clusters <- sort(unique(as.character(values$subset_markers[[current_subset_name]]$cluster)))
          updateSelectInput(session, "clusterSelect_final", 
                            choices = clusters, 
                            selected = clusters[1])
        }
      }
    })
    
    # Show completion message
    successful_subsets <- length(values$subset_markers)
    total_subsets <- length(values$subset_objects)
    
    showNotification(
      paste("Marker calculation complete!", successful_subsets, "of", total_subsets, "subsets processed"), 
      type = "message"
    )
  })
  
  # Render marker table for final clustering
  output$markerTable_final <- renderDT({
    req(values$subset_markers, input$subsetSelector, input$clusterSelect_final)
    
    if(is.null(values$subset_markers[[input$subsetSelector]])) {
      return(NULL)
    }
    
    marker_genes <- values$subset_markers[[input$subsetSelector]] %>% 
      filter(cluster == input$clusterSelect_final) %>%
      arrange(desc(avg_log2FC))
    
    display_table <- data.frame(
      Gene = marker_genes$gene,
      "Log2 FC" = round(marker_genes$avg_log2FC, 3),
      "Adj. p-value" = sprintf("%.2f", marker_genes$p_val_adj),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      display_table,
      options = list(
        dom = 't',
        scrollY = "350px",
        scroller = TRUE,
        paging = FALSE,
        ordering = TRUE,
        columnDefs = list(
          list(targets = 0, className = "dt-left"),
          list(targets = 1:2, className = "dt-right")
        )
      ),
      rownames = FALSE,
      selection = "none",
      class = "compact stripe"
    )
  })
  
  # UMAP visualization for final clustering
  output$clusterUMAP_final <- renderPlot({
    req(current_subset())
    
    if(!"umap" %in% names(current_subset()@reductions)) {
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "UMAP not computed yet. Please run LSI Round 2 first.") +
          theme_void() +
          theme(
            plot.title = element_text(size = 16, hjust = 0.5),
            aspect.ratio = 1
          )
      )
    }
    
    common_theme <- theme(
      legend.position = "right",
      aspect.ratio = 1,
      plot.margin = margin(10, 10, 10, 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
    
    if(input$umapDisplayType_final == "clusters") {
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        p <- DimPlot(current_subset(), 
                     reduction = "umap", 
                     group.by = input$umapGroupBy_final, 
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Colored by:", input$umapGroupBy_final))
      } else {
        p <- DimPlot(current_subset(), 
                     reduction = "umap", 
                     group.by = input$umapGroupBy_final, 
                     split.by = input$umapSplitBy_final,
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Colored by:", input$umapGroupBy_final, "| Split by:", input$umapSplitBy_final))
      }
    } else {
      req(input$featureSelect_final)
      if(is.null(input$featureSelect_final) || input$featureSelect_final == "" || 
         !input$featureSelect_final %in% rownames(current_subset())) {
        p <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Selected gene not found in dataset") +
          theme_void() +
          common_theme
      } else if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        p <- FeaturePlot(current_subset(), 
                         features = input$featureSelect_final, 
                         min.cutoff = "q10", 
                         max.cutoff = "q90",
                         pt.size = 0.8) +
          scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                               name = "Expression") +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Expression:", input$featureSelect_final))
      } else {
        tryCatch({
          split_factor <- current_subset()@meta.data[[input$umapSplitBy_final]]
          unique_splits <- unique(split_factor)
          
          plot_list <- lapply(unique_splits, function(split_val) {
            cells_subset <- colnames(current_subset())[split_factor == split_val]
            
            FeaturePlot(current_subset(), 
                        features = input$featureSelect_final,
                        cells = cells_subset,
                        min.cutoff = "q10", 
                        max.cutoff = "q90",
                        pt.size = 0.8) +
              scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                                   name = "Expression") +
              common_theme +
              ggtitle(paste(split_val))
          })
          
          p <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
          p <- cowplot::plot_grid(
            cowplot::ggdraw() + 
              cowplot::draw_label(paste(input$subsetSelector, "- Expression:", input$featureSelect_final, "| Split by:", input$umapSplitBy_final), 
                                  fontface = 'bold', size = 14),
            p,
            ncol = 1,
            rel_heights = c(0.1, 1)
          )
        }, error = function(e) {
          p <- ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = paste("Error creating split plot:", e$message)) +
            theme_void() +
            common_theme
        })
      }
    }
    
    return(p)
  }, height = 550, width = 700)
  
  # Cell type UI for final clustering
  output$cellTypeUI_final <- renderUI({
    req(current_subset())
    
    clusters <- sort(as.numeric(as.character(levels(Idents(current_subset())))))
    
    tagList(
      fluidRow(
        column(12, 
               h4(paste("Assign final cell types to", input$subsetSelector, "clusters"), 
                  style = "margin-bottom: 20px;")
        )
      ),
      fluidRow(
        column(12,
               div(
                 style = "max-height: 400px; overflow-y: auto; border: 1px solid #ddd; border-radius: 5px; padding: 15px; background-color: #f9f9f9;",
                 lapply(clusters, function(cluster) {
                   # Get top markers for this cluster
                   top_markers <- ""
                   if(!is.null(values$subset_markers) && !is.null(values$subset_markers[[input$subsetSelector]])) {
                     top_markers <- values$subset_markers[[input$subsetSelector]] %>% 
                       filter(cluster == !!cluster) %>% 
                       dplyr::top_n(3, avg_log2FC) %>% 
                       pull(gene) %>% 
                       paste(collapse = ", ")
                   }
                   
                   div(
                     style = "margin-bottom: 8px; display: flex; align-items: center;",
                     div(
                       style = "flex: 0 0 150px; font-weight: bold;",
                       paste("Cluster", cluster)
                     ),
                     div(
                       style = "flex: 0 0 250px; font-size: 0.8em; color: #666;",
                       if(top_markers != "") paste("Top markers:", top_markers) else ""
                     ),
                     div(
                       style = "flex: 1;",
                       textInput(
                         inputId = paste0("label_cluster_final_", input$subsetSelector, "_", cluster),
                         label = NULL,
                         placeholder = "Enter specific cell type name",
                         value = ifelse(
                           is.null(current_subset()$cell_type_final) || 
                             !any(current_subset()$seurat_clusters == cluster),
                           "", 
                           unique(current_subset()$cell_type_final[current_subset()$seurat_clusters == cluster])
                         ),
                         width = "100%"
                       )
                     )
                   )
                 })
               )
        )
      ),
      fluidRow(
        column(12,
               div(
                 style = "margin-top: 20px;",
                 actionButton("applyLabels_final", "Apply Final Cell Type Labels",
                              class = "btn-success", 
                              style = "color: white; width: 100%;")
               )
        )
      )
    )
  })
  
  
  
  
  
  # Apply final cell type labels
  observeEvent(input$applyLabels_final, {
    req(current_subset(), input$subsetSelector)
    
    withProgress(message = paste('Applying labels to', input$subsetSelector, '...'), {
      # Get current clusters for this subset
      clusters <- sort(levels(Idents(current_subset())))
      
      # Initialize cell_type_final column if it doesn't exist
      if(!"cell_type_final" %in% colnames(current_subset()@meta.data)) {
        values$subset_objects[[input$subsetSelector]]$cell_type_final <- 
          paste0("Cluster ", current_subset()$seurat_clusters)
      }
      
      # Apply new cell type labels
      for(cluster in clusters) {
        input_id <- paste0("label_cluster_final_", input$subsetSelector, "_", cluster)
        label_value <- input[[input_id]]
        
        # Only update if not empty
        if(!is.null(label_value) && nchar(trimws(label_value)) > 0) {
          values$subset_objects[[input$subsetSelector]]$cell_type_final[
            values$subset_objects[[input$subsetSelector]]$seurat_clusters == cluster] <- label_value
        }
      }
      
      showNotification(paste("Final cell type labels applied to", input$subsetSelector), type = "message")
    })
  })
  
  # Integration Results Logic
  observeEvent(input$integrateResults, {
    req(values$seurat, values$subset_objects)
    
    withProgress(message = 'Integrating final cell type annotations...', {
      
      # Initialize the final cell type column in the original object
      values$seurat$cell_type_final <- "Unassigned"
      
      # Counter for assigned cells
      total_assigned <- 0
      
      # Iterate through each subset and transfer labels
      for(subset_name in names(values$subset_objects)) {
        subset_obj <- values$subset_objects[[subset_name]]
        
        # Check if final cell types exist for this subset
        if("cell_type_final" %in% colnames(subset_obj@meta.data)) {
          
          # Get cell barcodes and their final annotations
          subset_cells <- colnames(subset_obj)
          subset_labels <- subset_obj$cell_type_final
          
          # Find matching cells in original object
          matching_cells <- intersect(subset_cells, colnames(values$seurat))
          
          if(length(matching_cells) > 0) {
            # Transfer labels to original object
            cell_indices <- match(matching_cells, colnames(values$seurat))
            label_indices <- match(matching_cells, subset_cells)
            
            values$seurat$cell_type_final[cell_indices] <- subset_labels[label_indices]
            total_assigned <- total_assigned + length(matching_cells)
          }
        }
      }
      
      # Summary statistics
      total_cells <- ncol(values$seurat)
      assigned_pct <- round(100 * total_assigned / total_cells, 2)
      
      # Show summary
      values$integration_summary <- list(
        total_cells = total_cells,
        assigned_cells = total_assigned,
        unassigned_cells = total_cells - total_assigned,
        assigned_percentage = assigned_pct,
        final_cell_types = table(values$seurat$cell_type_final)
      )
      
      showNotification(
        paste("Integration complete!", total_assigned, "cells assigned (", assigned_pct, "%)"), 
        type = "message"
      )
    })
  })
  
  # Display integration summary
  output$integrationSummary <- renderPrint({
    req(values$integration_summary)
    
    cat("Final Cell Type Integration Summary\n")
    cat("==================================\n\n")
    cat(paste("Total cells in dataset:", values$integration_summary$total_cells, "\n"))
    cat(paste("Cells with final annotations:", values$integration_summary$assigned_cells, "\n"))
    cat(paste("Unassigned cells:", values$integration_summary$unassigned_cells, "\n"))
    cat(paste("Assignment percentage:", values$integration_summary$assigned_percentage, "%\n\n"))
    
    cat("Final cell type distribution:\n")
    print(values$integration_summary$final_cell_types)
  })
  
  # Display final UMAP with integrated results
  output$finalIntegratedUMAP <- renderPlot({
    req(values$seurat, values$integration_summary)
    
    if(!"umap" %in% names(values$seurat@reductions)) {
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "UMAP not available in original object") +
          theme_void()
      )
    }
    
    DimPlot(values$seurat, 
            reduction = "umap", 
            group.by = "cell_type_final", 
            pt.size = 0.4,
            label = TRUE,
            label.size = 6,
            repel = TRUE) +
      theme(
        legend.position = "right",
        aspect.ratio = 1,
        plot.title = element_text(size = 14, hjust = 0.5)
      ) +
      ggtitle("Final Integrated Cell Type Annotations")
  }, height = 400, width = 600)
  
  
  # Download Logic
  observeEvent(input$downloadIntegratedData, {
    # Check if output folder is set
    if(!exists("output_folder") || is.null(output_folder) || output_folder == "") {
      showNotification("Please set an output directory in the Setup tab first", type = "error")
      return()
    }
    
    # Check if output folder exists
    if(!dir.exists(output_folder)) {
      showNotification("Output directory does not exist. Please select a valid directory in the Setup tab", type = "error")
      return()
    }
    
    # Check if any items are selected
    if(is.null(input$finalExportItems) || length(input$finalExportItems) == 0) {
      showNotification("Please select at least one item to export", type = "warning")
      return()
    }
    
    withProgress(message = 'Saving files to output directory...', value = 0, {
      
      # Create a timestamped subfolder for this export
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      export_folder <- file.path(output_folder, paste0("GUIdedRNA_Results_", timestamp))
      
      # Create the export directory
      dir.create(export_folder, recursive = TRUE, showWarnings = FALSE)
      
      total_items <- length(input$finalExportItems)
      progress_increment <- 1 / total_items
      
      files_saved <- c()
      
      if("final_rds" %in% input$finalExportItems && !is.null(values$seurat)) {
        incProgress(progress_increment, detail = "Saving Seurat object...")
        rds_file <- file.path(export_folder, "processed_annotated_integrated_seurat.rds")
        saveRDS(values$seurat, rds_file)
        files_saved <- c(files_saved, "Seurat Object (RDS)")
      }
      
      if("count_matrix" %in% input$finalExportItems && !is.null(values$seurat)) {
        incProgress(progress_increment, detail = "Saving count matrix...")
        matrix_file <- file.path(export_folder, "processed_count_matrix.mtx")
        Matrix::writeMM(GetAssayData(values$seurat, layer = "counts"), matrix_file)
        files_saved <- c(files_saved, "Count Matrix (MTX)")
      }
      
      if("final_meta" %in% input$finalExportItems && !is.null(values$seurat)) {
        incProgress(progress_increment, detail = "Saving metadata...")
        meta_file <- file.path(export_folder, "cell_type_metadata.csv")
        write.csv(values$seurat@meta.data, meta_file, row.names = TRUE)
        files_saved <- c(files_saved, "Cell Type Metadata (CSV)")
      }
      
      if("subset_markers" %in% input$finalExportItems && !is.null(values$subset_markers)) {
        incProgress(progress_increment, detail = "Saving marker genes...")
        marker_count <- 0
        for(subset_name in names(values$subset_markers)) {
          # Clean subset name for filename (remove special characters)
          clean_name <- gsub("[^A-Za-z0-9_-]", "_", subset_name)
          marker_file <- file.path(export_folder, paste0("markers_", clean_name, ".csv"))
          write.csv(values$subset_markers[[subset_name]], marker_file, row.names = FALSE)
          marker_count <- marker_count + 1
        }
        files_saved <- c(files_saved, paste0("Subset Marker Genes (", marker_count, " files)"))
      }
      
      if("dimred_coords" %in% input$finalExportItems && !is.null(values$seurat) && "umap" %in% names(values$seurat@reductions)) {
        incProgress(progress_increment, detail = "Saving UMAP coordinates...")
        coords <- as.data.frame(values$seurat@reductions$umap@cell.embeddings)
        coords_file <- file.path(export_folder, "dimensionality_reduction_coordinates.csv")
        write.csv(coords, coords_file, row.names = TRUE)
        files_saved <- c(files_saved, "UMAP Coordinates (CSV)")
      }
      
      if("integration_report" %in% input$finalExportItems && !is.null(values$integration_summary)) {
        incProgress(progress_increment, detail = "Saving integration report...")
        report_file <- file.path(export_folder, "integration_summary_report.txt")
        
        # Create the report content
        report_content <- paste(
          "GUIdedRNA Analysis Results",
          "==========================",
          "",
          paste("Analysis completed:", Sys.time()),
          paste("Export directory:", export_folder),
          "",
          "Final Cell Type Integration Summary",
          "==================================",
          "",
          paste("Total cells in dataset:", values$integration_summary$total_cells),
          paste("Cells with final annotations:", values$integration_summary$assigned_cells),
          paste("Unassigned cells:", values$integration_summary$unassigned_cells),
          paste("Assignment percentage:", values$integration_summary$assigned_percentage, "%"),
          "",
          "Final cell type distribution:",
          paste(capture.output(print(values$integration_summary$final_cell_types)), collapse = "\n"),
          "",
          "Files exported:",
          paste("-", files_saved, collapse = "\n"),
          sep = "\n"
        )
        
        writeLines(report_content, report_file)
        files_saved <- c(files_saved, "Integration Report (TXT)")
      }
      
      # Create a summary file with export information
      summary_file <- file.path(export_folder, "export_summary.txt")
      summary_content <- paste(
        "GUIdedRNA Export Summary",
        "========================",
        "",
        paste("Export timestamp:", timestamp),
        paste("Export directory:", export_folder),
        paste("Total files exported:", length(files_saved)),
        "",
        "Exported files:",
        paste("-", files_saved, collapse = "\n"),
        "",
        "File descriptions:",
        "- processed_annotated_integrated_seurat.rds: Complete Seurat object with all analyses",
        "- processed_count_matrix.mtx: Raw count matrix in Matrix Market format",
        "- cell_type_metadata.csv: Cell metadata including final annotations",
        "- markers_[celltype].csv: Differentially expressed genes for each cell type subset",
        "- dimensionality_reduction_coordinates.csv: UMAP coordinates for visualization",
        "- integration_summary_report.txt: Detailed analysis summary and statistics",
        sep = "\n"
      )
      
      writeLines(summary_content, summary_file)
    })
    
    # Show success notification with file count and location
    showNotification(
      HTML(paste0(
        "<strong>Export completed successfully!</strong><br>",
        "Files saved to: ", export_folder, "<br>",
        "Total files: ", length(files_saved) + 1  # +1 for the summary file
      )),
      type = "message",
      duration = 10  # Show for 10 seconds
    )
  })

# Download UMAP for Final Clustering
observeEvent(input$downloadUMAP_final, {
  # Check if output folder is set
  if(!exists("output_folder") || is.null(output_folder) || output_folder == "") {
    showNotification("Please set an output directory in the Setup tab first", type = "error")
    return()
  }
  
  # Check if output folder exists
  if(!dir.exists(output_folder)) {
    showNotification("Output directory does not exist. Please select a valid directory in the Setup tab", type = "error")
    return()
  }
  
  # Check if current subset and UMAP exist
  if(is.null(current_subset()) || !"umap" %in% names(current_subset()@reductions)) {
    showNotification("No UMAP available to download. Please run LSI Round 2 first.", type = "error")
    return()
  }
  
  withProgress(message = 'Saving UMAP figure...', {
    # Create timestamped filename
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    subset_name_clean <- gsub("[^A-Za-z0-9_-]", "_", input$subsetSelector)
    
    # Create descriptive filename based on current settings
    if(input$umapDisplayType_final == "clusters") {
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_UMAP_", input$umapGroupBy_final, "_", timestamp, ".png")
      } else {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_UMAP_", input$umapGroupBy_final, "_split_by_", input$umapSplitBy_final, "_", timestamp, ".png")
      }
    } else {
      feature_name <- gsub("[^A-Za-z0-9_-]", "_", input$featureSelect_final)
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_Expression_", feature_name, "_", timestamp, ".png")
      } else {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_Expression_", feature_name, "_split_by_", input$umapSplitBy_final, "_", timestamp, ".png")
      }
    }
    
    # Full file path
    file_path <- file.path(output_folder, filename)
    
    # Recreate the current plot
    common_theme <- theme(
      legend.position = "right",
      aspect.ratio = 1,
      plot.margin = margin(10, 10, 10, 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
    
    if(input$umapDisplayType_final == "clusters") {
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        p <- DimPlot(current_subset(), 
                     reduction = "umap", 
                     group.by = input$umapGroupBy_final, 
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Colored by:", input$umapGroupBy_final))
      } else {
        p <- DimPlot(current_subset(), 
                     reduction = "umap", 
                     group.by = input$umapGroupBy_final, 
                     split.by = input$umapSplitBy_final,
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Colored by:", input$umapGroupBy_final, "| Split by:", input$umapSplitBy_final))
      }
    } else {
      req(input$featureSelect_final)
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        p <- FeaturePlot(current_subset(), 
                         features = input$featureSelect_final, 
                         min.cutoff = "q10", 
                         max.cutoff = "q90",
                         pt.size = 0.8) +
          scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                               name = "Expression") +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Expression:", input$featureSelect_final))
      } else {
        split_factor <- current_subset()@meta.data[[input$umapSplitBy_final]]
        unique_splits <- unique(split_factor)
        
        plot_list <- lapply(unique_splits, function(split_val) {
          cells_subset <- colnames(current_subset())[split_factor == split_val]
          
          FeaturePlot(current_subset(), 
                      features = input$featureSelect_final,
                      cells = cells_subset,
                      min.cutoff = "q10", 
                      max.cutoff = "q90",
                      pt.size = 0.8) +
            scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                                 name = "Expression") +
            common_theme +
            ggtitle(paste(split_val))
        })
        
        p <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
        p <- cowplot::plot_grid(
          cowplot::ggdraw() + 
            cowplot::draw_label(paste(input$subsetSelector, "- Expression:", input$featureSelect_final, "| Split by:", input$umapSplitBy_final), 
                                fontface = 'bold', size = 14),
          p,
          ncol = 1,
          rel_heights = c(0.1, 1)
        )
      }
    }
    
    # Save the plot
    ggsave(file_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  })
  
  showNotification(
    HTML(paste0(
      "<strong>Figure saved successfully!</strong><br>",
      "File: ", basename(file_path), "<br>",
      "Location: ", output_folder
    )),
    type = "message",
    duration = 5
  )
})


# Download UMAP for Initial Clustering
observeEvent(input$downloadUMAP_initial, {
  # Check if output folder is set
  if(!exists("output_folder") || is.null(output_folder) || output_folder == "") {
    showNotification("Please set an output directory in the Setup tab first", type = "error")
    return()
  }
  
  # Check if output folder exists
  if(!dir.exists(output_folder)) {
    showNotification("Output directory does not exist. Please select a valid directory in the Setup tab", type = "error")
    return()
  }
  
  # Check if seurat object and UMAP exist
  if(is.null(values$seurat) || !"umap" %in% names(values$seurat@reductions)) {
    showNotification("No UMAP available to download. Please run LSI first.", type = "error")
    return()
  }
  
  withProgress(message = 'Saving UMAP figure...', {
    # Create timestamped filename
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    # Create descriptive filename based on current settings
    if(input$umapDisplayType == "clusters") {
      if(is.null(input$umapSplitBy) || input$umapSplitBy == "none") {
        filename <- paste0("Initial_Clustering_UMAP_", input$umapGroupBy, "_", timestamp, ".png")
      } else {
        filename <- paste0("Initial_Clustering_UMAP_", input$umapGroupBy, "_split_by_", input$umapSplitBy, "_", timestamp, ".png")
      }
    } else {
      feature_name <- gsub("[^A-Za-z0-9_-]", "_", input$featureSelect)
      if(is.null(input$umapSplitBy) || input$umapSplitBy == "none") {
        filename <- paste0("Initial_Clustering_Expression_", feature_name, "_", timestamp, ".png")
      } else {
        filename <- paste0("Initial_Clustering_Expression_", feature_name, "_split_by_", input$umapSplitBy, "_", timestamp, ".png")
      }
    }
    
    # Full file path
    file_path <- file.path(output_folder, filename)
    
    # Recreate the current plot
    common_theme <- theme(
      legend.position = "right",
      aspect.ratio = 1,
      plot.margin = margin(10, 10, 10, 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
    
    if(input$umapDisplayType == "clusters") {
      if(is.null(input$umapSplitBy) || input$umapSplitBy == "none") {
        p <- DimPlot(values$seurat, 
                     reduction = "umap", 
                     group.by = input$umapGroupBy, 
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste("Colored by:", input$umapGroupBy))
      } else {
        p <- DimPlot(values$seurat, 
                     reduction = "umap", 
                     group.by = input$umapGroupBy, 
                     split.by = input$umapSplitBy,
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste("Colored by:", input$umapGroupBy, "| Split by:", input$umapSplitBy))
      }
    } else {
      req(input$featureSelect)
      if(is.null(input$umapSplitBy) || input$umapSplitBy == "none") {
        p <- FeaturePlot(values$seurat, 
                         features = input$featureSelect, 
                         min.cutoff = "q10", 
                         max.cutoff = "q90",
                         pt.size = 0.8) +
          scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                               name = "Expression") +
          common_theme +
          ggtitle(paste("Expression:", input$featureSelect))
      } else {
        split_factor <- values$seurat@meta.data[[input$umapSplitBy]]
        unique_splits <- unique(split_factor)
        
        plot_list <- lapply(unique_splits, function(split_val) {
          cells_subset <- colnames(values$seurat)[split_factor == split_val]
          
          FeaturePlot(values$seurat, 
                      features = input$featureSelect,
                      cells = cells_subset,
                      min.cutoff = "q10", 
                      max.cutoff = "q90",
                      pt.size = 0.8) +
            scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                                 name = "Expression") +
            common_theme +
            ggtitle(paste(split_val))
        })
        
        p <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
        p <- cowplot::plot_grid(
          cowplot::ggdraw() + 
            cowplot::draw_label(paste("Expression:", input$featureSelect, "| Split by:", input$umapSplitBy), 
                                fontface = 'bold', size = 14),
          p,
          ncol = 1,
          rel_heights = c(0.1, 1)
        )
      }
    }
    
    # Save the plot
    ggsave(file_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  })
  
  showNotification(
    HTML(paste0(
      "<strong>Figure saved successfully!</strong><br>",
      "File: ", basename(file_path), "<br>",
      "Location: ", output_folder
    )),
    type = "message",
    duration = 5
  )
})

# Download UMAP for Final Clustering
observeEvent(input$downloadUMAP_final, {
  # Check if output folder is set
  if(!exists("output_folder") || is.null(output_folder) || output_folder == "") {
    showNotification("Please set an output directory in the Setup tab first", type = "error")
    return()
  }
  
  # Check if output folder exists
  if(!dir.exists(output_folder)) {
    showNotification("Output directory does not exist. Please select a valid directory in the Setup tab", type = "error")
    return()
  }
  
  # Check if current subset and UMAP exist
  if(is.null(current_subset()) || !"umap" %in% names(current_subset()@reductions)) {
    showNotification("No UMAP available to download. Please run LSI Round 2 first.", type = "error")
    return()
  }
  
  withProgress(message = 'Saving UMAP figure...', {
    # Create timestamped filename
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    subset_name_clean <- gsub("[^A-Za-z0-9_-]", "_", input$subsetSelector)
    
    # Create descriptive filename based on current settings
    if(input$umapDisplayType_final == "clusters") {
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_UMAP_", input$umapGroupBy_final, "_", timestamp, ".png")
      } else {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_UMAP_", input$umapGroupBy_final, "_split_by_", input$umapSplitBy_final, "_", timestamp, ".png")
      }
    } else {
      feature_name <- gsub("[^A-Za-z0-9_-]", "_", input$featureSelect_final)
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_Expression_", feature_name, "_", timestamp, ".png")
      } else {
        filename <- paste0("Final_Clustering_", subset_name_clean, "_Expression_", feature_name, "_split_by_", input$umapSplitBy_final, "_", timestamp, ".png")
      }
    }
    
    # Full file path
    file_path <- file.path(output_folder, filename)
    
    # Recreate the current plot
    common_theme <- theme(
      legend.position = "right",
      aspect.ratio = 1,
      plot.margin = margin(10, 10, 10, 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
    
    if(input$umapDisplayType_final == "clusters") {
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        p <- DimPlot(current_subset(), 
                     reduction = "umap", 
                     group.by = input$umapGroupBy_final, 
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Colored by:", input$umapGroupBy_final))
      } else {
        p <- DimPlot(current_subset(), 
                     reduction = "umap", 
                     group.by = input$umapGroupBy_final, 
                     split.by = input$umapSplitBy_final,
                     pt.size = 0.8,
                     label = TRUE,
                     label.size = 5) +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Colored by:", input$umapGroupBy_final, "| Split by:", input$umapSplitBy_final))
      }
    } else {
      req(input$featureSelect_final)
      if(is.null(input$umapSplitBy_final) || input$umapSplitBy_final == "none") {
        p <- FeaturePlot(current_subset(), 
                         features = input$featureSelect_final, 
                         min.cutoff = "q10", 
                         max.cutoff = "q90",
                         pt.size = 0.8) +
          scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                               name = "Expression") +
          common_theme +
          ggtitle(paste(input$subsetSelector, "- Expression:", input$featureSelect_final))
      } else {
        split_factor <- current_subset()@meta.data[[input$umapSplitBy_final]]
        unique_splits <- unique(split_factor)
        
        plot_list <- lapply(unique_splits, function(split_val) {
          cells_subset <- colnames(current_subset())[split_factor == split_val]
          
          FeaturePlot(current_subset(), 
                      features = input$featureSelect_final,
                      cells = cells_subset,
                      min.cutoff = "q10", 
                      max.cutoff = "q90",
                      pt.size = 0.8) +
            scale_color_gradient(low = "#00bfc4", high = "#f8766d", na.value = "grey50",
                                 name = "Expression") +
            common_theme +
            ggtitle(paste(split_val))
        })
        
        p <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
        p <- cowplot::plot_grid(
          cowplot::ggdraw() + 
            cowplot::draw_label(paste(input$subsetSelector, "- Expression:", input$featureSelect_final, "| Split by:", input$umapSplitBy_final), 
                                fontface = 'bold', size = 14),
          p,
          ncol = 1,
          rel_heights = c(0.1, 1)
        )
      }
    }
    
    # Save the plot
    ggsave(file_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  })
  
  showNotification(
    HTML(paste0(
      "<strong>Figure saved successfully!</strong><br>",
      "File: ", basename(file_path), "<br>",
      "Location: ", output_folder
    )),
    type = "message",
    duration = 5
  )
})
 

}
 
# Run the application
shinyApp(ui = ui, server = server, options = list(
  width = 1200,
  height = 1200
))