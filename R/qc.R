

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







remove_doublets <- function(object){



}


