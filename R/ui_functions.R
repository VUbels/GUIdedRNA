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