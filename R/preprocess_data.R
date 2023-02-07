

#' Preprocess single cell sequencing data
#'
#' @param object A seurat object
#'
#' @return A seurat object with preprocessed data
#' @export
#'
#' @examples
#' preprocess_data(object = sobj)


preprocess_data <- function(object = object){
  Seurat::DefaultAssay(object) <- "RNA"

  object <- object %>% Seurat::NormalizeData() %>%
                Seurat::FindVariableFeatures() %>%
                Seurat::ScaleData()

  return(object)
}
