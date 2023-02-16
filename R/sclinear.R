

#' Prepare data for modality prediction
#'
#' @param object A Seurat object
#'
#' @return object A Seurat object containing with annotated cell types
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj)
#' }

prepare_data <- function(object = object, remove_doublets = FALSE, low_qc_cell_removal = FALSE, anno_level = 2, samples = samples){

  Seurat::DefaultAssay(object) <- "RNA"

  object <- object %>% preprocess_data()

  if(remove_doublets){
    object <- object %>% remove_doublets()
  }

  if(low_qc_cell_removal){
    object <- object %>% mad_filtering()
  }

  object <- object %>% anno_celltypes(anno_level = 2)

  object <- object %>% Seurat::RunPCA(npcs = 100 ) %>% Seurat::RunUMAP(dims = 1:30)


  p1 <- Seurat::DimPlot(object, group.by = "cell_type", label = TRUE, repel = TRUE) + ggplot2::theme(legend.position = "null")
  p2 <- Seurat::DimPlot(object, group.by = "original_cell_type", label = TRUE, repel = TRUE) + ggplot2::theme(legend.position = "null")

  print(p1 + p2)



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

  object <- object %>% subset(subset = cell_type == "T")

  gexp_matrix <- Matrix::as.matrix(object@assays$RNA@counts)

  gexp_matrix_py <- reticulate::r_to_py(gexp_matrix)

  pipe <- ADTPredictor(do_log1p = FALSE)
  pipe.fit((gexp_matrix_py))

}
