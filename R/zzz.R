numpy <- NULL
anndata <- NULL
scanpy <- NULL
joblib <- NULL
os <- NULL
torch <- NULL
pytorch_lightning <- NULL
sklearn <- NULL
warnings <- NULL
scipy <- NULL
typing <- NULL


preprocessing <- NULL
evaluate <- NULL
prediction <- NULL


#' onLoad function
#'
#' @param libname A
#' @param pkgname A
#'
#' @return NULL
.onLoad <- function(libname, pkgname){

  if(reticulate::py_module_available("numpy")){numpy <<- reticulate::import("numpy", delay_load = TRUE)}
  if(reticulate::py_module_available("joblib")){joblib <<- reticulate::import("joblib", delay_load = TRUE)}
  if(reticulate::py_module_available("torch")){torch <<- reticulate::import("torch", delay_load = TRUE)}
  if(reticulate::py_module_available("pytorch_lightning")){pytorch_lightning <<- reticulate::import("pytorch_lightning", delay_load = TRUE)}
  if(reticulate::py_module_available("sklearn")){sklearn <<- reticulate::import("sklearn", delay_load = TRUE)}
  if(reticulate::py_module_available("anndata")){anndata <<- reticulate::import("anndata", delay_load = TRUE)}
  if(reticulate::py_module_available("scanpy")){scanpy <<- reticulate::import("scanpy", delay_load = TRUE)}
  if(reticulate::py_module_available("os")){os <<- reticulate::import("os", delay_load = TRUE)}
  if(reticulate::py_module_available("warnings")){warnings <<- reticulate::import("warnings", delay_load = TRUE)}
  if(reticulate::py_module_available("scipy")){scipy <<- reticulate::import("scipy", delay_load = TRUE)}
  if(reticulate::py_module_available("typing")){typing <<- reticulate::import("typing", delay_load = TRUE)}


  # load package specific python modules
  module_path <-  base::system.file("python",package = utils::packageName())
  if(module_path == ""){module_path <-  base::system.file("inst/python",package = utils::packageName())}

  if(all(
    reticulate::py_module_available("numpy"),
    reticulate::py_module_available("sklearn"),
    reticulate::py_module_available("anndata"),
    reticulate::py_module_available("scanpy"),
    reticulate::py_module_available("warnings"),
    reticulate::py_module_available("typing")
  )){preprocessing <<- reticulate::import_from_path("preprocessing",module_path,delay_load = TRUE)}else{
    print("Loading preprocessing module failed. Check if all necessary python dependencies are installed: https://github.com/DanHanh/scLinear")
  }

  if(all(
    reticulate::py_module_available("sklearn"),
    reticulate::py_module_available("scipy")
  )){evaluate <<- reticulate::import_from_path("evaluate",module_path,delay_load = TRUE)}else{
    print("Loading evaluate module failed. Check if all necessary python dependencies are installed: https://github.com/DanHanh/scLinear")
  }

  if(all(
    reticulate::py_module_available("numpy"),
    reticulate::py_module_available("joblib"),
    reticulate::py_module_available("pytorch_lightning"),
    reticulate::py_module_available("torch"),
    reticulate::py_module_available("sklearn"),
    reticulate::py_module_available("anndata"),
    reticulate::py_module_available("os"),
    reticulate::py_module_available("warnings"),
    reticulate::py_module_available("scipy"),
    reticulate::py_module_available("typing")
  )){prediction <<- reticulate::import_from_path("prediction",module_path,delay_load = TRUE)}else{
    print("Loading prediction module failed. Check if all necessary python dependencies are installed: https://github.com/DanHanh/scLinear")
  }
}
