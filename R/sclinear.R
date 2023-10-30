

#' Prepare data for modality prediction
#'
#' @param object A Seurat object
#' @param remove_doublets Should doublets be removed? (TRUE,FALSE)
#' @param low_qc_cell_removal Should low quality cells be removed using median absolute deviation?
#' @param anno_level Level of annotation used for cell type annotation. Either a single number or a vector with multiple numbers c(1,2,3,4).
#' @param samples Variable that contains the sample names. If set, will use rpca to integrate the samples.
#' @param integrate_data Should the sample be integrated through the rpca method?
#' @param remove_empty_droplets If the raw unfiltered matrix was used to create the Seurat object, empty droplets should be filtered out. (TRUE, FALSE)
#' @param lower Lower boundary to define empty droplets. All cells with a lower amount of nFeatures are assumed to be empty.
#' @param FDR FDR threshold to define non empty cells.
#' @param annotation_selfCluster Should the clusters determined by Seurat be used for cell type annotation.
#' @param resolution Resolution for louvain clustering.
#' @param seed Used seed.
#' @param return_plots Should plots be returned from function
#' @param print_plots Print plots. (TRUE, FALSE)
#' @param species Species. Relevant for cell type annotation. ("Hs", "Mm")
#'
#' @return object A pre-processed Seurat object  with annotated cell types
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, resolution = 0.8)
#' }

prepare_data <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE,remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = FALSE, resolution = 0.8, seed = 42, return_plots = FALSE, print_plots = TRUE, species = "Hs"){
  set.seed(seed)

  plot_list <- list()

  Seurat::DefaultAssay(object) <- "RNA"

  if(remove_empty_droplets){
    object <- empty_drops(object = object, lower = lower, FDR = FDR, samples = samples, seed = seed)
    #plot_list[["empty_dropts"]] <- object[[2]]
    #object <- object[[1]]
  }

  if(!("mito_percent" %in% names(object@meta.data))){
    object$mito_percent <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
  }

  if(remove_doublets){
    print("Start remove doublets")
    object <- object %>% remove_doublets(samples = samples, print_plots = print_plots, seed = seed)
    plot_list[["doublets"]] <- object[[2]]
    object <- object[[1]]
  }

  if(low_qc_cell_removal){
    print("Start low quality cell removal")
    object <- object %>% mad_filtering(samples = samples, print_plots = print_plots, seed = seed)
    plot_list[["low_qc_cells"]] <- object[[2]]
    object <- object[[1]]
  }


  # object <- object %>% Seurat::NormalizeData() %>%
  #   Seurat::FindVariableFeatures() %>%
  #   Seurat::ScaleData()

  if(integrate_data){
    print("Start integrate data")
    object <- integrate_samples(object, samples = samples, seed = seed)
  }

  print("Start clustering data")
  object <- cluster_data(object, resolution = resolution, seed = seed)
  Seurat::Idents(object) <- object@meta.data[["seurat_clusters"]]

  print("Start cell type annotation")
  if(annotation_selfCluster){
    object <- object %>% anno_celltypes(anno_level = anno_level, selfClusters = Seurat::Idents(.), species = species, seed = seed)
  }else{
    object <- object %>% anno_celltypes(anno_level = anno_level, selfClusters = NULL, species = species, seed = seed)
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
scLinear <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = FALSE, resolution = 0.8, seed = 42, return_plots = FALSE, model = "all", assay_name = "RNA", print_plots = FALSE, species = "Hs"){
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
                         print_plots = print_plots,
                         species = species)

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
#' @param gexp Matrix with gene expression data
#' @param pipe Trained ADT predictor
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
#' @param pipe Trained ADT predictor
#' @param gexp_test Matrix with gene expression data.
#' @param adt_test Matrix with ADT count data.
#' @param do_log1p A
#'
#' @return Returns a data frame containing RSME, Pearson correlation and Spearman correlation of the tested data.
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
#' @param pipe A ADT predictor object
#' @param model Choose a pre-trained model to load. The pre-traines models were
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



#' Normalize gene expression matrix with scran and scuttle
#'
#' @param gexp_matrix A gene expression matrix
#' @param center.size.factors A
#' @param log A
#' @param ... For the method, additional arguments passed to logNormCounts.
#'
#' @return Normalized expression matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Normalize expression matirx
#' normalized_matrix <- gexp_normalize(sobj\@assays[["RNA"]]\@counts)
#' # Add normalized matrix back to RNA assay in Seurat.
#' sobj\@assays[["RNA"]]\@data <- normalized_matrix
#' }
gexp_normalize <- function(gexp_matrix, center.size.factors = FALSE, log = FALSE, ...){
  ## normalize data GEX
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gexp_matrix))
  clusters <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters=clusters)
  sce <- scuttle::logNormCounts(sce, center.size.factors = center.size.factors, log = log, ...)
  gexp_matrix <- sce@assays@data@listData[["normcounts"]]
  gexp_matrix <- base::log1p(gexp_matrix)
  return(gexp_matrix)
}

