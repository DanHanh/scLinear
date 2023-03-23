

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

prepare_data <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE,remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = FALSE, resolution = 0.8, seed = 42){
  set.seed(seed)

  Seurat::DefaultAssay(object) <- "RNA"

  if(remove_empty_droplets){

    object <- empty_drops(object = object, lower = lower, FDR = FDR)
  }

  if(!("mito_percent" %in% names(object@meta.data))){
    object$mito_percent <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
  }


  object <- object %>% Seurat::NormalizeData() %>%
                        Seurat::FindVariableFeatures() %>%
                        Seurat::ScaleData()

  if(remove_doublets){
    print("Start remove doublets")
    object <- object %>% remove_doublets(samples = samples)
  }

  if(low_qc_cell_removal){
    print("Start low quality cell removal")
    object <- object %>% mad_filtering(samples = samples)
  }

  if(integrate_data){
    print("Start integrate data")
    object <- integrate_samples(object, samples = samples, resolution = resolution)
  }

  print("Start clustering data")
  object <- cluster_data(object, resolution = resolution)
  Seurat::Idents(object) <- object@meta.data[["seurat_clusters"]]

  print("Start cell type annotation")
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
    pipe <- prediction$ADTPredictor(do_log1p = do_log1p)
  return(pipe)
}

#' Train a predictor object
#'
#' @param pipe A predictor object
#' @param gex_train Gene expression assay
#' @param adt_train ADT assay
#' @param normalize Normalize GEX and ATD assay before fitting.
#'
#' @return pipe A trained predictor object
#' @export
#'
#' @examples
#' \dontrun{
#' fit_predictor(pipe = pipe, gex_train = object@assays$RNA , adt_train = object@assays$ADT)
#' }
fit_predictor <- function(pipe, gexp_train , adt_train, normalize = TRUE){

  gexp_matrix <- as.matrix(gexp_train@counts)
  adt_matrix <- as.matrix(adt_train@counts)

  if(normalize){
    ## normalize data GEX
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gexp_matrix))
    clusters <- scran::quickCluster(sce)
    sce <- scran::computeSumFactors(sce, clusters=clusters)
    sce <- scuttle::logNormCounts(sce, pseudo.count = 1, center.size.factors = FALSE, log = FALSE)
    gexp_matrix <- sce@assays@data@listData[["normcounts"]]
    gexp_matrix <- base::log1p(gexp_matrix)
  }

  if(normalize){
    adt_matrix <- Seurat::NormalizeData(adt_matrix, normalization.method = "CLR", margin = 2)
  }

  gexp_matrix_py <- reticulate::r_to_py(t(gexp_matrix))
  adt_matrix_py <- reticulate::r_to_py(t(adt_matrix))

  pipe$fit(gexp_matrix_py, adt_matrix_py, gex_names = rownames(gexp_matrix), adt_names = rownames(adt_matrix))

  return(pipe)
}


#' Predict ADT values from gene expression
#'
#' @param gexp A
#' @param pipe A
#' @param do_log1p A
#'
#' @return adt_assay retuns an adt assay object
#' @export
#'
#' @examples
#' \dontrun{
#' adt_predict(gextp)
#' }
adt_predict <- function(pipe, gexp, normalize = TRUE){

  gexp_matrix <- gexp@counts

  if(normalize){
    ## normalize data GEX
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gexp_matrix))
    clusters <- scran::quickCluster(sce)
    sce <- scran::computeSumFactors(sce, clusters=clusters)
    sce <- scuttle::logNormCounts(sce, pseudo.count = 1, center.size.factors = FALSE, log = FALSE)
    gexp_matrix <- sce@assays@data@listData[["normcounts"]]
    gexp_matrix <- base::log1p(gexp_matrix)
  }


  gexp_matrix <- t(as.matrix(gexp_matrix))

  gexp_matrix_py <- reticulate::r_to_py(gexp_matrix)

  predicted_adt <- pipe$predict(gexp_matrix_py, gex_names = colnames(gexp_matrix))

  ## adt matrix
  adt <- as.matrix(predicted_adt[[1]])
  ## names of predicted proteins
  adt_names <- predicted_adt[[2]] #$to_list()

  ## add
  colnames(adt) <- adt_names
  ## add initial cell names
  rownames(adt) <- rownames(gexp_matrix)
  ## transpose back for assay
  adt <- t(adt)

  adt_assay <- Seurat::CreateAssayObject(data = adt)

  return(adt_assay)
}

#' Evaluate the adt predictor
#'
#' @param pipe A
#' @param gexp_test A
#' @param adt_test A
#' @param do_log1p A
#'
#' @return A
#' @export
#'
#' @examples
#' \dontrun{
#' evaluate_predictor(pipe, gex_test, adt_test)
#' }
evaluate_predictor <- function(pipe, gexp_test, adt_test, normalize = TRUE){

  predicted_adt <- adt_predict(pipe, gexp_test, normalize = normalize)

  ## subset to only comone features
  p_adt <- subset(predicted_adt,features = which(rownames(predicted_adt) %in% rownames(adt_test)) )
  t_adt <- subset(adt_test,features = which(rownames(adt_test) %in% rownames(predicted_adt)) )


  ### CLR transform data test data
  if(normalize){
    t_adt <- Seurat::NormalizeData(t_adt, normalization.method = "CLR", margin = 2)
  }
  ## transpose to fit anndata format
  p_adt_matrix <- t(as.matrix(p_adt@data))
  t_adt_matrix <- t(as.matrix(t_adt@data))

  ## reorder adt text matrix to the same order as predicted adt
  t_adt_matrix <- t_adt_matrix[,match(colnames(p_adt_matrix), colnames(t_adt_matrix))]

  p_adt_matrix_py <- reticulate::r_to_py(p_adt_matrix)
  t_adt_matrix_py<- reticulate::r_to_py(t_adt_matrix)

  ev_res <- evaluate$evaluate(p_adt_matrix_py, t_adt_matrix_py)

  return(ev_res)
}


#' Load a pretrained model
#'
#' @param pipe A pipe
#' @param model Choose pretrained model: all, bcells, tcells, nkcells
#'
#' @return pipe
#' @export
#'
#' @examples
#' \dontrun{
#' load_pretrained_model(pipe, model = "all")
#' }
load_pretrained_model <- function(pipe, model = "all"){

  load_path <-  base::system.file("python",package = "scLinear")


  m <- switch(model,
           "all" = "ADTPredictor_neuripstrain_alltypes.joblib",
           "bcell" = "ADTPredictor_neuripstrain_Bcells.joblib",
           "nkcell" = "ADTPredictor_neuripstrain_NKcells.joblib",
           "tcell" = "ADTPredictor_neuripstrain_Tcells.joblib")


  pipe$load(paste0(load_path,"/",m))


  return(pipe)

}
