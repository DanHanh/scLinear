

#' Prepare data for modality prediction
#'
#' @param object A Seurat object
#'
#' @return object A pre-processed Seurat object  with annotated cell types
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, resolution = 0.8)
#' }

prepare_data <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE,remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = FALSE, resolution = 0.8, seed = 42, return_plots = FALSE, print_plots = TRUE){
  set.seed(seed)

  plot_list <- list()

  Seurat::DefaultAssay(object) <- "RNA"

  if(remove_empty_droplets){
    object <- empty_drops(object = object, lower = lower, FDR = FDR, samples = samples)
    plot_list[["empty_dropts"]] <- object[[2]]
    object <- object[[1]]
  }

  if(!("mito_percent" %in% names(object@meta.data))){
    object$mito_percent <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
  }

  if(remove_doublets){
    print("Start remove doublets")
    object <- object %>% remove_doublets(samples = samples, print_plots = print_plots)
    plot_list[["doublets"]] <- object[[2]]
    object <- object[[1]]
  }

  if(low_qc_cell_removal){
    print("Start low quality cell removal")
    object <- object %>% mad_filtering(samples = samples, print_plots = print_plots)
    plot_list[["low_qc_cells"]] <- object[[2]]
    object <- object[[1]]
  }


  # object <- object %>% Seurat::NormalizeData() %>%
  #   Seurat::FindVariableFeatures() %>%
  #   Seurat::ScaleData()

  if(integrate_data){
    print("Start integrate data")
    object <- integrate_samples(object, samples = samples)
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
  if(print_plots){base::print(p1)}

  if(return_plots){
    return_object <- list(object = object, plots = plot_list)
  }else{
    return_object <- object
  }


  return(return_object)

}


#' Predict modalities based on gene expression data
#'
#' @param object A Seurat object
#'
#' @return object A Seurat object containing additional single cell modalities
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj)
#' }
scLinear <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = FALSE, resolution = 0.8, seed = 42, return_plots = FALSE, model = "all", assay_name = "RNA", print_plots = FALSE){
  set.seed(seed)
  object <- prepare_data(object,
                         remove_doublets = remove_doublets,
                         low_qc_cell_removal = low_qc_cell_removal,
                         anno_level = anno_level,
                         samples = samples,
                         integrate_data = integrate_data,
                         remove_empty_droplets = remove_empty_droplets,
                         lower = lower,
                         FDR = FDR,
                         annotation_selfCluster = annotation_selfCluster,
                         resolution = resolution,
                         seed = seed,
                         return_plots = FALSE,
                         print_plots = print_plots)

  pipe <- create_adt_predictor()
  pipe <- load_pretrained_model(pipe, model = model)

  object@assays["predicted_ADT"] <-  adt_predict(pipe = pipe,
                                                  gexp = Seurat::GetAssay(object, assay = assay_name),
                                                  normalize = TRUE)

  return(object)

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
#' create_adt_predictor()
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
fit_predictor <- function(pipe, gexp_train , adt_train, slot_gex = "counts", slot_adt = "counts", normalize_gex = TRUE, normalize_adt = TRUE){


  gexp_matrix <- Seurat::GetAssayData(gexp_train, slot = slot_gex)
  adt_matrix <- Seurat::GetAssayData(adt_train, slot = slot_adt)



  if(normalize_gex){
    gexp_matrix <- gexp_normalize(gexp_matrix)
  }
  if(normalize_adt){
    adt_matrix <- Seurat::NormalizeData(adt_matrix, normalization.method = "CLR", margin = 2)
  }

  gexp_matrix_py <- reticulate::r_to_py(Matrix::t(gexp_matrix))
  adt_matrix_py <- reticulate::r_to_py(Matrix::t(adt_matrix))

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
adt_predict <- function(pipe, gexp, slot = "counts", normalize = TRUE){

  gexp_matrix <- Seurat::GetAssayData(gexp, slot = slot)

  if(normalize){
    gexp_matrix <- gexp_normalize(gexp_matrix)
  }
  gexp_matrix <- Matrix::t(gexp_matrix)

  gexp_matrix_py <- reticulate::r_to_py(as.matrix(gexp_matrix))

  predicted_adt <- pipe$predict(gexp_matrix_py, gex_names = colnames(gexp_matrix))

  ## adt matrix
  adt <- predicted_adt[[1]]
  ## names of predicted proteins
  if(typeof(predicted_adt[[2]]) == "environment"){
    adt_names <- predicted_adt[[2]]$to_list()
  }else{
    adt_names <- predicted_adt[[2]]
  }
  ## add
  colnames(adt) <- adt_names
  ## add initial cell names
  rownames(adt) <- rownames(gexp_matrix)
  ## transpose back for assay
  adt <- Matrix::t(adt)

  adt_assay <- Seurat::CreateAssayObject(data = adt)

  #Seurat::Key(adt_assay) <- "predictedadt_"

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
evaluate_predictor <- function(pipe, gexp_test, adt_test, slot = "counts", normalize_gex = TRUE, normalize_adt = TRUE){

  ### CLR transform test data
  if(normalize_adt){
    adt_test <- Seurat::NormalizeData(adt_test, normalization.method = "CLR", margin = 2)
  }


  predicted_adt <- adt_predict(pipe, gexp_test, slot = slot,  normalize = normalize_gex)

  ## subset to features found in predicted and in test matrix
  p_adt <- subset(predicted_adt,features = which(rownames(predicted_adt) %in% rownames(adt_test)) )
  t_adt <- subset(adt_test,features = which(rownames(adt_test) %in% rownames(predicted_adt)) )

  ## transpose to fit anndata format
  p_adt_matrix <- Matrix::t(p_adt@data)
  t_adt_matrix <- Matrix::t(t_adt@data)

  ## reorder adt text matrix to the same order as predicted adt
  t_adt_matrix <- t_adt_matrix[,match(colnames(p_adt_matrix), colnames(t_adt_matrix))]

  p_adt_matrix_py <- reticulate::r_to_py(p_adt_matrix)
  t_adt_matrix_py<- reticulate::r_to_py(t_adt_matrix)

  ev_res <- evaluate$evaluate(p_adt_matrix_py, t_adt_matrix_py)
  return_df <- data.frame(RMSE = ev_res[[1]], Pearson = ev_res[[2]], Spearman = ev_res[[3]])

  return(return_df)
}


#' Load a pre-trained model
#'
#' @param pipe A pipe
#' @param model Choose a pre-trained model to load. The pretraines models were
#' trained on the NeurIPS data. In most cases we would recommend to use the model ("all") trained on all
#' cell types. models available: all, bcells, tcells, nkcells.
#'
#' @return pipe Returns the pipe with a loaded pre-trained model.
#' @export
#'
#' @examples
#' \dontrun{
#' load_pretrained_model(pipe, model = "all")
#' }
load_pretrained_model <- function(pipe, model = "all"){

  load_path <-  base::system.file("data",package = "scLinear")

  m <- switch(model,
           "all" = "ADTPredictor_neuripstrain_alltypes.joblib",
           "bcell" = "ADTPredictor_neuripstrain_Bcells.joblib",
           "nkcell" = "ADTPredictor_neuripstrain_NKcells.joblib",
           "tcell" = "ADTPredictor_neuripstrain_Tcells.joblib")

  pipe$load(paste0(load_path,"/",m))

  return(pipe)

}

