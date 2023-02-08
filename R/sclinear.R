

#' Predict modalities based on gene expression data
#'
#' @param object A Seurat object
#'
#' @return A Seurat object containing additional single cell modalities
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj)
#' }

scLinear <- function(object = object, remove_doublets = FALSE, low_qc_cell_removal = FALSE, anno_level = 2){

  reticulate::source_python("./inst/python/example.py")


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


  p1 <- DimPlot(sobj, group.by = "cell_type", label = TRUE, repel = TRUE) + theme(legend.position = "null")
  p2 <- DimPlot(sobj, group.by = "original_cell_type", label = TRUE, repel = TRUE) + theme(legend.position = "null")

  print(p1 + p2 )

  return(object)

}
