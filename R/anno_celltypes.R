

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

anno_celltypes <- function(object, anno_level = 2, species = "Hs", ...){
  ## load Human Panglao database
  base::load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA"))


  anno_res <- scMRMA::scMRMA(input = object,
                    species = species,
                    db = "panglaodb",
                    ...
                     )

  object@meta.data$cell_type <-
          anno_res$multiR$annotationResult[[paste0("Level",anno_level)]]


  return(object)

}
