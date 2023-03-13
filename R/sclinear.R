

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

prepare_data <- function(object = object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, annotation_selfCluster = FALSE, resolution = 0.8, seed = 42){
  set.seed(seed)

  Seurat::DefaultAssay(object) <- "RNA"

  object <- object %>% Seurat::NormalizeData() %>%
                        Seurat::FindVariableFeatures() %>%
                        Seurat::ScaleData()

  if(remove_doublets){
    object <- object %>% remove_doublets(samples = samples)
  }

  if(low_qc_cell_removal){
    object <- object %>% mad_filtering(samples = samples)
  }

  if(integrate_data){
    object <- integrate_samples(object, samples = samples, resolution = resolution)
  }

  object <- cluster_data(object, resolution = resolution)
  Seurat::Idents(object) <- object@meta.data[["seurat_clusters"]]

  if(annotation_selfCluster){
    object <- object %>% anno_celltypes(anno_level = anno_level, selfClusters = Seurat::Idents(.))
  }else{
    object <- object %>% anno_celltypes(anno_level = anno_level, selfClusters = NULL)
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
#' @return adt_assay retuns an adt assay object
#' @export
#'
#' @examples
#' \dontrun{
#' adt_predict(gextp)
#' }
adt_predict <- function(pipe, gexp){

  gexp_matrix <- t(as.matrix(gexp@counts))
  gexp_matrix_py <- reticulate::r_to_py(gexp_matrix)

  predicted_adt <- pipe$predict(gexp_matrix_py, gex_names = colnames(gexp_matrix))

  ## adt matrix
  adt <- as.matrix(predicted_adt[[1]])
  ## names of predicted proteins
  adt_names <- predicted_adt[[2]]$to_list()

  ## add
  colnames(adt) <- adt_names
  ## add initial cell names
  rownames(adt) <- rownames(gexp_matrix)
  ## transpose back for assay
  adt <- t(adt)

  adt_assay <- Seurat::CreateAssayObject(adt)

  return(adt_assay)
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

  gexp_matrix <- t(as.matrix(gexp@counts))
  gexp_matrix_py <- reticulate::r_to_py(gexp_matrix)

  ## predict adt from test set
  predicted_adt <- pipe$predict(gexp_matrix_py, gex_names = colnames(gexp_matrix))
  ## adt matrix
  adt_res <- as.matrix(predicted_adt[[1]])
  ## names of predicted proteins
  adt_res_names <- predicted_adt[[2]]$to_list()
  ## add
  colnames(adt_res) <- adt_res_names
  ## add initial cell names
  rownames(adt_res) <- rownames(gexp_matrix)

  ## preprocess adt test matrix
  adt_test_matrix <- t(as.matrix(adt_test@counts))

  ## subset pipe and test marix to the intersection proteins predicted and in test matrix
  adt_res <- adt_res[,colnames(adt_res) %in% colnames(adt_matrix)]
  adt_test_matrix <- adt_test_matrix[,colnames(adt_test_matrix) %in% colnames(adt_res)]

  ## reorder adt text matrix to the same order as predicted adt
  adt_test_matrix <- adt_test_matrix[,match(colnames(adt_res), colnames(adt_test_matrix))]

  ## lognormalize values
  adt_res <- log1p(adt_res/rowSums(adt_res) * 10000)
  adt_test_matrix <- log1p(adt_test_matrix/rowSums(adt_test_matrix) * 10000)


  adt_res_py <- reticulate::r_to_py(adt_res)
  adt_test_matrix_py<- reticulate::r_to_py(adt_test_matrix)


  ev_res <- evaluate$evaluate(adt_res_py, adt_test_matrix_py)

  return(ev_res)
}
