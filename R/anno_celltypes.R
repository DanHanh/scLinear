

#' Annotate each cell with the corresponding cell type
#'
#' @param object Seurat object
#'
#' @return A Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- anno_celltypes(object = sobj, anno_level = 2, species = "Hs")
#' }

anno_celltypes <- function(object, anno_level = 2, selfClusters = NULL ,species = "Hs", ...){
  Seurat::DefaultAssay(object) <- "RNA"

  ## load Human Panglao database
  base::load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA"))


  anno_res <- scMRMA::scMRMA(input = object,
                    species = species,
                    db = "panglaodb",
                    selfClusters = selfClusters,
                    ...
                     )
  if(length(anno_level) == 1){
    object@meta.data$cell_type <-
            anno_res$multiR$annotationResult[[paste0("Level",anno_level)]]
  }else{
    for(i in anno_level){
      object@meta.data[[paste0("cell_type_", i)]] <-
        anno_res$multiR$annotationResult[[paste0("Level",i)]]

    }
  }


  return(object)

}



#' Create A UMAP from the annotated Seurat object
#'
#' @param object A Seurat object
#' @param group.by Paraneter by witch to group visualization by
#' @param ndims Number of dimensions for dimensional reduction
#'
#' @return p A Dim plot
#' @export
#'
#' @examples
#' \dontrun{
#' p <- visualize_data(sobj, group.by = "cell_type", ndims = 30)
#' }

visualize_data <- function(object, group.by = "cell_type", ndims = NULL, ...){



  object <- object %>% Seurat::RunPCA(npcs = max(ndims, 100) )

  ## estimate dimensionality of data if no ndims is supplied
  if(is.null(ndims)){
  ndims <- ceiling(intrinsicDimension::maxLikGlobalDimEst(object@reductions[[paste0("pca")]]@cell.embeddings, k = 20)[["dim.est"]])
  }

  object <- object %>% Seurat::RunUMAP(dims = 1:ndims)

  p <- Seurat::DimPlot(object, group.by = group.by, ...)
  base::print(p)

  return(p)
}




cluster_with_integration <- function(object, method = "rpca", samples = "samples", resolution  = 0.8){

  if(typeof(object) == "list"){
    object_list <- object
  }else{
    Seurat::SplitObject(sobj, split.by = "batch")
    object_list <- Seurat::SplitObject(object, split.by = samples)
  }

  ## normalize data and find variable features for each sample
  object_list <- lapply(object_list, function(ob){
    ob <- ob %>% Seurat::NormalizeData() %>% Seurat::FindVariableFeatures();
    ob
  })

  ## select features that are repeatedly variable across samples
  features <- Seurat::SelectIntegrationFeatures(object.list = object_list)

  ## scale and calculate PCs for each sample
  object_list <- lapply(object_list, function(ob) {
    ob <- ob %>% Seurat::ScaleData(verbose = FALSE, features = features) %>%
      Seurat::RunPCA(verbose = FALSE, npcs = 100, features = features);
    ob
  })


  ## find integration anchors
  integration.anchors <- Seurat::FindIntegrationAnchors(object.list = object_list, anchor.features = features, reduction = "rpca")

  object_integrated <- Seurat::IntegrateData(anchorset = integration.anchors)

  Seurat::DefaultAssay(object_integrated) <- "integrated"

  object_integrated <- object_integrated %>%  Seurat::ScaleData() %>%
                          Seurat::RunPCA(assay = "integrated", npcs = 100)

  ## estimate the optimal number of PCs for clustering
  ndims <- ceiling(intrinsicDimension::maxLikGlobalDimEst(object_integrated@reductions[[paste0("pca")]]@cell.embeddings, k = 20)[["dim.est"]])

  object_integrated <- object_integrated %>% Seurat::FindNeighbors(dims = 1:ndims, reduction= "pca")

  object_integrated <- object_integrated %>% Seurat::FindClusters(resolution = resolution)

  p <- visualize_data(object_integrated, group.by = "seurat_clusters")

  return(object_integrated)

}

