

#' Remove low quality cells based on the median absolute deviation.
#'
#' @param object Seurat object
#' @param nmads Median absolute deviation that should be counted as outliner
#' @param type If upper, lower , both boundaries should be used for filtering
#' @param mttype If upper, lower , both boundaries should be used for mitochondrial percent filtering
#'
#' @return object Seurat object
#' @export
#'
#' @examples
#' sobj <- mad_filtering(sobj)

mad_filtering <- function(object = objec, nmads = 3, type = "both", mttype = "lower"){

  nCount_ol <- scater::isOutlier(object@meta.data$nCount_RNA,
                                 nmads = nmads,
                                 type = type,
                                 log = TRUE)

  nFeature_ol <- scater::isOutlier(object@meta.data$nFeature_RNA,
                                   nmads = nmads,
                                   type = type,
                                   log = TRUE)

  pMito_ol <- scater::isOutlier(object@meta.data$nFeature_RNA,
                               nmads = nmads,
                               type = mttype,
                               log = FALSE)



  object@meta.data$mad_filtered <- (nCount_ol | nFeature_ol | pMito_ol)



  ## filter visualization

  DF <- data.frame(nCount_RNA = object@meta.data$nCount_RNA, nFeature_RNA = object@meta.data$nFeature_RNA , mito_percent = object@meta.data$mito_percent, Filtered = object@meta.data$mad_filtered)
  p <- ggplot(DF, aes(x = nCount_RNA, y = nFeature_RNA , color = Filtered)) +
    geom_point() + theme_bw() +
    scale_color_manual(values = c("darkgreen", "darkred")) +
    ggtitle("MAD filtered cells") +
    labs(caption = paste0("#filtered cells: ", sum(DF$Filtered =="TRUE"), " from ", length(DF$Filtered), "cells", "\n",
                          "nCount_th: ", round(attr(nCount_ol, "thresholds")["lower"], digits = 1), ", " ,round(attr(nCount_ol, "thresholds")["higher"], digits = 1), "\n",
                          "nFeature_th: ", round(attr(nFeature_ol, "thresholds")["lower"], digits = 1), ", " ,round(attr(nFeature_ol, "thresholds")["higher"], digits = 1), "\n",
                          "pMito_th: ", round(attr(pMito_ol, "thresholds")["lower"], digits = 1), ", " ,round(attr(pMito_ol, "thresholds")["higher"], digits = 1))
         ) +
    geom_vline(xintercept = attr(nCount_ol, "thresholds"), color = 'red', linetype = "dashed") +
    geom_hline(yintercept = attr(nFeature_ol, "thresholds"), color = 'red', linetype = "dashed") +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10')
  print(p)


  ## remove filtered cells

  object <- subset(object, subset = mad_filtered == "FALSE")
  return(object)
}







#' Remove Heterotypic doublets from sequencing data
#'
#' @param object A Seurat object
#' @param samples meta.data information with  the sequencing samples
#'
#' @return object A Seurat object with removed doublets
#' @export
#'
#' @examples
#' remove_doublets(object = object, samples = "batch")

remove_doublets <- function(object = object, samples = "samples"){

  #### remove dublicates with scDblFinder
  sce <- as.SingleCellExperiment(object)
  sce <- scDblFinder::scDblFinder(sce, clusters = NULL, samples = samples)
  a <- Seurat::as.Seurat(sce)


  stat <-  table(object@meta.data[["scDblFinder.class"]])

  object_visualization <- object %>%  preprocess_data() %>% Seurat::RunPCA(npcs = 100 ) %>% Seurat::RunUMAP(dims = 1:30)

  p <- DimPlot(object_visualization, reduction = "umap", group.by = "scDblFinder.class") +
    labs(title = "Detected doublets",
         caption = paste0("#singlets: ", stat[["singlet"]], " ; ",
                          "#doublets: ", stat[["doublet"]] ))
  print(p)

  object <- subset(object, subset = scDblFinder.class == "singlet")


  ## remove unnecessary information from Seurat object
  object@meta.data[["scDblFinder.cluster"]] <- NULL
  object@meta.data[["scDblFinder.class"]] <- NULL
  object@meta.data[["scDblFinder.score"]] <- NULL
  object@meta.data[["scDblFinder.weighted"]] <- NULL
  object@meta.data[["scDblFinder.difficulty"]] <- NULL
  object@meta.data[["scDblFinder.cxds_score"]] <- NULL
  object@meta.data[["scDblFinder.mostLikelyOrigin"]] <- NULL
  object@meta.data[["scDblFinder.originAmbiguous"]] <- NULL

  return(object)
}


