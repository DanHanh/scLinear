

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
#' \dontrun{
#' sobj <- mad_filtering(sobj)
#' }

mad_filtering <- function(object = objec, samples = NULL, nmads = 3, type = "both", mttype = "higher", remove_cells = TRUE, ...){

  ##
  if(is.null(samples)){
    batch <- NULL
  }else{
    batch <- object@meta.data %>% dplyr::pull(samples)
  }

  nCount_ol <- scater::isOutlier(object@meta.data$nCount_RNA,
                                 nmads = nmads,
                                 type = type,
                                 log = TRUE,
                                 batch = batch
                                 )

  nFeature_ol <- scater::isOutlier(object@meta.data$nFeature_RNA,
                                   nmads = nmads,
                                   type = type,
                                   log = TRUE,
                                   batch = batch
                                   )

  pMito_ol <- scater::isOutlier(object@meta.data$mito_percent,
                               nmads = nmads,
                               type = mttype,
                               log = FALSE,
                               batch = batch
                               )


  object@meta.data$mad_filtered <- (nCount_ol | nFeature_ol | pMito_ol)

  ## filter visualization
  metadata <- object@meta.data %>% dplyr::select(tidyselect::any_of(c("nCount_RNA", "nFeature_RNA", "mito_percent", "mad_filtered" ))) %>% dplyr::rename(Filtered = "mad_filtered")
  p <- ggplot2::ggplot(metadata, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA , color = Filtered)) +
    ggplot2::geom_point() + ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("darkgreen", "darkred")) +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::scale_y_continuous(trans='log10')
  base::print(p)

  if(is.null(samples)){
    p <- ggplot2::ggplot(metadata, ggplot2::aes(x = "", fill = Filtered, label = ggplot2::after_stat(count))) + ggplot2::theme_bw() +
      ggplot2::geom_bar(position = "identity", stat = "count") + ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
      ggplot2::geom_text(stat = "count", vjust = -1) + ggplot2::labs(title = "Number of filtered cells by sample", fill = "Filtered") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + ggplot2::xlab("Sample") + ggplot2::ylab("# cells")
    base::print(p)

  }else{
    metadata <- cbind(metadata, data.frame(samples = object@meta.data %>% dplyr::pull(samples)))


  p <- ggplot2::ggplot(metadata, ggplot2::aes(x = samples, fill = Filtered, label = ggplot2::after_stat(count))) + ggplot2::theme_bw() +
    ggplot2::geom_bar(position = "identity", stat = "count") + ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
    ggplot2::geom_text(stat = "count", vjust = -1) + ggplot2::labs(title = "Number of quality filtered cells by sample", fill = "Filtered") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + ggplot2::xlab("Sample") + ggplot2::ylab("# cells")
  base::print(p)
  }

  if(remove_cells){
    ## remove filtered cells
    object <- base::subset(object, subset = mad_filtered == "FALSE")
    object@meta.data[["mad_filtered"]] <- NULL
  }

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
#' \dontrun{
#' remove_doublets(object = object, samples = "batch")
#' }

remove_doublets <- function(object = object, samples = NULL, remove_cells = TRUE ,seed = 42, ...){
  set.seed(seed = seed)

  #### remove doublets with scDblFinder
  sce <- Seurat::as.SingleCellExperiment(object)
  sce <- scDblFinder::scDblFinder(sce, samples = samples, ...)

  if(is.null(samples)){
    metadata <- sce@colData@listData %>% as.data.frame() %>%
      dplyr::select("scDblFinder.class")

    p <- ggplot2::ggplot(metadata, ggplot2::aes(x = "", fill = scDblFinder.class, label = ggplot2::after_stat(count))) + ggplot2::theme_bw() +
      ggplot2::geom_bar(position = "identity", stat = "count") + ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
      ggplot2::geom_text(stat = "count", vjust = -1) + ggplot2::labs(title = "Number of doublets by sample", fill = "Type") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + ggplot2::xlab("Sample") + ggplot2::ylab("# cells")
    base::print(p)

  }else{
    metadata <- sce@colData@listData %>% as.data.frame() %>%
                  dplyr::select("scDblFinder.sample", "scDblFinder.class")

    p <- ggplot2::ggplot(metadata, ggplot2::aes(x = scDblFinder.sample, fill = scDblFinder.class, label = ggplot2::after_stat(count))) + ggplot2::theme_bw() +
      ggplot2::geom_bar(position = "identity", stat = "count") + ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
      ggplot2::geom_text(stat = "count", vjust = -1) + ggplot2::labs(title = "Number of doublets by sample", fill = "Type") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + ggplot2::xlab("Sample") + ggplot2::ylab("# cells")
    base::print(p)
}
  ## add singlet/doublet information to initial Seurat object
  object[["scDblFinder.class"]] <- sce@colData@listData[["scDblFinder.class"]]

  if(remove_cells){
    ## remove doublets
    object <- base::subset(object, subset = scDblFinder.class == "singlet")

    ## remove unnecessary information from Seurat object
    object@meta.data[["scDblFinder.class"]] <- NULL
  }

  return(object)
}


#' Remove empty droplets
#'
#' @param object A Seurat object
#' @param lower A
#' @param FDR A
#'
#' @return object A Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' empty_drops(object = object, lower = 100, FDR = 0.01)
#' }
empty_drops <- function(object, lower = 100, FDR = 0.01, samples = NULL){

  object <- pbmc10k
  length(colnames(object))
  object@meta.data[["samples"]] <- c(rep("s1", times = 8420/2), rep("s2", times = 8420/2))


  object_split <- Seurat::SplitObject(object, split.by = samples)

  for (i in names(object_split)){
    print(i)

    e.out <- DropletUtils::emptyDrops(object_split[[i]]@assays$RNA@counts, lower = lower)
    is.cell <- e.out$FDR <= FDR
    object_split[[i]] <- object_split[[i]][, which(is.cell)]

  }

  object <- merge(object_split[[1]], y = object_split[2:length(object_split)])
  return(object)
}
