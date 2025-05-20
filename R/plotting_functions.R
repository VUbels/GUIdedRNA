#' Function to run UMAP on LSI results.
#'
#' @param df Seurat object.
#' @param rawCounts count matrix from Seurat object.
#' @param cmap colour map to pass to function.
#' @param coverLabel cover label.
#' @param point.size individual point size per cell.
#' @param namedColors vector colour map.
#' @param plotTitle title to plot.
#' @param colorLims colour limits for minimal and maximal percentage based on expression data.
#' @param na.value colour for cells not expression gene.
#' @param useRaster raster plot.

plotUMAP <- function(df, dataType = "qualitative", cmap = NULL, covarLabel = "", point_size=0.5, 
                     namedColors=FALSE, plotTitle=NULL, colorLims=NULL, na.value="grey35", useRaster=TRUE){
  # Given a df containing the UMAP x and y coords and a third column, 
  # plot the UMAP
  
  if(useRaster){
    p <- (
      ggplot(df, aes(x = df[,1], y = df[,2], color = df[,3]))
      + geom_point_rast(size = point_size)
      + theme_BOR()
      + theme(
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        aspect.ratio=1
      )
      + xlab("UMAP1")
      + ylab("UMAP2")
    )
  }else{
    p <- (
      ggplot(df, aes(x = df[,1], y = df[,2], color = df[,3]))
      + geom_point(size = point_size)
      + theme_BOR()
      + theme(
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        aspect.ratio=1
      )
      + xlab("UMAP1")
      + ylab("UMAP2")
    )
  }
  
  # Set plot title
  if(!is.null(plotTitle)){
    p <- p + ggtitle(plotTitle)
  }else{
    p <- p + ggtitle(sprintf("n = %s", nrow(df)))
  }
  # If colormap provided, update colors
  if(!is.null(cmap)){
    if(namedColors){
      # Named colormap corresponding to discrete values in third column
      p <- p + scale_color_manual(values=cmap, limits=names(cmap), name=covarLabel, na.value=na.value)
      p <- p + guides(fill = guide_legend(title=covarLabel), 
                      colour = guide_legend(override.aes = list(size=5)))
    }else{
      # Remove names
      names(cmap) <- NULL
      if(dataType == "qualitative"){
        # check to make sure you have enough colors for qualitative mapping
        nvals <- length(unique(df[,3]))
        cmap <- getColorMap(cmap, n=nvals)
        p <- p + scale_color_manual(values=cmap, name=covarLabel, na.value=na.value)
        p <- p + guides(fill = guide_legend(title=covarLabel), 
                        colour = guide_legend(override.aes = list(size=5)))
      }else{
        if(!is.null(colorLims)){
          p <- p + scale_color_gradientn(colors=cmap, name=covarLabel, limits=colorLims, na.value=na.value)
        }else{
          p <- p + scale_color_gradientn(colors=cmap, name=covarLabel, na.value=na.value)
        }
      }
    }
  }
  p
}