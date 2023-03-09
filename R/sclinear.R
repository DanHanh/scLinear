

#' Prepare data for modality prediction
#'
#' @param object A Seurat object
#'
#' @return object A Seurat object containing with annotated cell types
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, cluster_with_integrated_data = FALSE, resolution = 0.8)
#' }

prepare_data <- function(object = object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, cluster_with_integrated_data = FALSE, resolution = 0.8, seed = 42){
  set.seed(seed)

  Seurat::DefaultAssay(object) <- "RNA"

  object <- object %>% preprocess_data()

  if(remove_doublets){
    object <- object %>% remove_doublets(samples = samples)
  }

  if(low_qc_cell_removal){
    object <- object %>% mad_filtering(samples = samples)
  }

  if(cluster_with_integrated_data){

    object <- cluster_with_integration(object, samples = samples, resolution = resolution)
    Seurat::Idents(object) <- object@meta.data[["seurat_clusters"]]
    object <- object %>% anno_celltypes(anno_level = 2, selfClusters = Seurat::Idents(.))

  }else{
    object <- object %>% anno_celltypes(anno_level = 2)
    object <- object %>% Seurat::RunPCA(npcs = 100 )
    ## estimate the optimal number of PCs for clustering
    ndims <- ceiling(intrinsicDimension::maxLikGlobalDimEst(object@reductions[[paste0("pca")]]@cell.embeddings, k = 20)[["dim.est"]])
    object <- object %>% Seurat::RunUMAP(dims = 1:ndims)

  }



  p1 <- Seurat::DimPlot(object, group.by = "cell_type", label = TRUE, repel = TRUE) + ggplot2::theme(legend.position = "null")
  base::print(p1)



  return(object)

}


#' Predict modalities based on gene expression data
#'
#' @param object
#'
#' @return object A Seurat object containing additional single cell modalities
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj)
#' }
scLinear <- function(object = object, cell_type){
  reticulate::source_python("./inst/python/preprocessing.py")
  reticulate::source_python("./inst/python/evaluate.py")
  reticulate::source_python("./inst/python/prediction.py")

  object <- object %>% base::subset(subset = cell_type == "T")

  gexp_matrix <- t(as.matrix(object@assays$RNA@counts))
  adt_matrix <- t(as.matrix(object@assays$ADT@counts))



  gexp_matrix_py <- reticulate::r_to_py(gexp_matrix)
  adt_matrix_py <- reticulate::r_to_py(adt_matrix)

  pipe <- ADTPredictor(do_log1p = FALSE)
  pipe$fit(gexp_matrix_py, adt_matrix)
  adt_pred <- pipe$predict(gexp_matrix_py)
  evaluate(adt_pred, adt_matrix_py)

}


#' Create a Gene expression to ADT assay predictor
#'
#' @param do_log1p
#'
#' @return pipe an adt predictor object
#' @export
#'
#' @examples
#' \dontrun{
#' create_adt_predictor(do_log1p = FALSE)
#' }
create_adt_predictor <- function(do_log1p = FALSE){
    reticulate::source_python("./inst/python/prediction.py")
    reticulate::source_python("./inst/python/preprocessing.py")
    pipe <- ADTPredictor(do_log1p = do_log1p)
  return(pipe)
}

#' Train a predictor object
#'
#' @param pipe A predictor object
#' @param gex_train Gene expression assay
#' @param adt_train ADT assay
#'
#' @return pipe A trained predictor object
#' @export
#'
#' @examples
#' \dontrun{
#' fit_predictor(pipe = pipe, = object@assays$RNA , adt_train = object@assays$ADT)
#' }
fit_predictor <- function(pipe, gex_train , adt_train){
  gexp_matrix <- as.matrix(gex_train@counts)
  adt_matrix <- as.matrix(adt_train@counts)

  gexp_matrix_py <- reticulate::r_to_py(t(gexp_matrix))
  adt_matrix_py <- reticulate::r_to_py(t(adt_matrix))


  pipe$fit(gexp_matrix_py, adt_matrix_py)

  return(pipe)
}


#' Predict ADT values from gene expression
#'
#' @param gexp
#'
#' @return A A
#' @export
#'
#' @examples
#' \dontrun{
#' adt_predict(gextp)
#' }
adt_predict <- function(pipe, gexp){
  reticulate::source_python("./inst/python/prediction.py")

  gexp_matrix <- as.matrix(gexp@counts)
  gexp_matrix_py <- reticulate::r_to_py(t(gexp_matrix))


  predicted_adt <- pipe$predict(gexp_matrix_py)

  predicted_adt <- t(predicted_adt)

  colnames(predicted_adt) <- colnames(gexp_matrix)



  Seurat::CreateAssayObject(counts)

}

#' Evaluate the adt predictor
#'
#' @param pipe A
#' @param gex_test A
#' @param adt_test A
#'
#' @return A
#' @export
#'
#' @examples
#' \dontrun{
#' evaluate_predictor(pipe, gex_test, adt_test)
#' }
evaluate_predictor <- function(pipe, gex_test, adt_test){
  reticulate::source_python("./inst/python/prediction.py")
  reticulate::source_python("./inst/python/evaluate.py")

  gexp_matrix <- as.matrix(gex_test@counts)
  adt_matrix <- as.matrix(adt_test@counts)

  gexp_matrix_py <- reticulate::r_to_py(t(gexp_matrix))
  adt_matrix_py <- reticulate::r_to_py(t(adt_matrix))

  predicted_adt <- pipe$predict(gexp_matrix_py)

  ev_res <- evaluate(predicted_adt, adt_matrix_py)

  return(ev_res)
}
