





preprocessing <- NULL
prediction <- NULL
evaluate <- NULL


.onLoad <- function(libname, pkgname){
  module_path <-  base::system.file("python",package = utils::packageName())
  preprocessing <<- reticulate::import_from_path("preprocessing",module_path,delay_load = TRUE)
  prediction <<- reticulate::import_from_path("prediction",module_path,delay_load = TRUE)
  evaluate <<- reticulate::import_from_path("evaluate",module_path,delay_load = TRUE)

  # utils::assignInMyNamespace("preprocessing",preprocessing)
  # utils::assignInMyNamespace("prediction",prediction)
  # utils::assignInMyNamespace("evaluate",evaluate)

}




# preprocessing <- reticulate::import_from_path("preprocessing",file.path("inst", "python"))
# prediction <- reticulate::import_from_path("prediction",file.path("inst", "python"))
# evaluate <- reticulate::import_from_path("evaluate",file.path("inst", "python"))
#
# for (obj in names(preprocessing)) {assign(obj, NULL)}
# for (obj in names(prediction)) {assign(obj, NULL)}
# for (obj in names(evaluate)) {assign(obj, NULL)}
#
#
# .onLoad <- function(libname, pkgname) {
#   preprocessing <- reticulate::import_from_path("preprocessing",
#                 system.file("python", package = packageName()),delay_load = TRUE)
#   prediction <- reticulate::import_from_path("prediction",
#                 system.file("python", package = packageName()),delay_load = TRUE)
#   evaluate <- reticulate::import_from_path("evaluate",
#                 system.file("python", package = packageName()),delay_load = TRUE)
#
#
#
#   for (obj in names(preprocessing)){assignInMyNamespace(obj, preprocessing[[obj]])}
#   for (obj in names(prediction)){assignInMyNamespace(obj, prediction[[obj]])}
#   for (obj in names(evaluate)){assignInMyNamespace(obj, evaluate[[obj]])}
#
#
#
# }



# numpy <- NULL
# scanpy <- NULL
# joblib <- NULL
# evaluate <- NULL
# prediction <- NULL
# preprocessing <- NULL
# torch <- NULL
# pytorch_lightning <- NULL
# sklearn <- NULL
# anndata <- NULL
#
# .onLoad <- function(libname, pkgname){
#     #reticulate::configure_environment(pkgname)
#     module_path <-  base::system.file("python",package = "scLinear")
#
#     #tryCatch({
#       numpy <<- reticulate::import("numpy", delay_load = TRUE)
#       joblib <<- reticulate::import("joblib", delay_load = TRUE)
#       torch <<- reticulate::import(module = "torch", delay_load = TRUE)
#       pytorch_lightning <<- reticulate::import(module = "pytorch_lightning", delay_load = TRUE)
#       sklearn <<- reticulate::import(module = "sklearn", delay_load = TRUE)
#       anndata <<- reticulate::import("anndata", delay_load = TRUE)
#       scanpy <<- reticulate::import("scanpy", delay_load = TRUE)
#
#       preprocessing <<- reticulate::import_from_path(module = "preprocessing", path = module_path, delay_load = TRUE)
#       prediction <<- reticulate::import_from_path(module = "prediction", path = module_path, delay_load = TRUE)
#       evaluate <<- reticulate::import_from_path(module = "evaluate", path = module_path, delay_load = TRUE)
#     #}, error = function(e){
#      # packageStartupMessage("Some python packages could not be loaded. Try install_pyton_dependencies to install missing dependencies!")
#     #  return(NULL)
#     #  }
#    # )
#
# }







#' Install all python dependencies for scLinear
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' install_pyton_dependencies()
#' }
install_pyton_dependencies <- function(){
  # Install python dependencies
  if(any(
    !all(sapply(c("numpy", "sklearn", "joblib", "anndata", "torch", "scanpy"), function(x){reticulate::py_module_available(x)}))
  )){
    install.dependencies <- readline('Do you want to install the necessary python packages (yes/no)?')
    if(install.dependencies == "y" || install.dependencies == "yes"){
      if(!reticulate::py_module_available("numpy")) suppressWarnings(suppressMessages(reticulate::conda_install(packages = "numpy")))
      if(!reticulate::py_module_available("sklearn")) suppressWarnings(suppressMessages(reticulate::conda_install(packages = "scikit-learn")))
      if(!reticulate::py_module_available("anndata")) suppressWarnings(suppressMessages(reticulate::conda_install(packages = "anndata")))
      if(!reticulate::py_module_available("joblib")) suppressWarnings(suppressMessages(reticulate::conda_install(packages = "joblib")))
      if(!reticulate::py_module_available("torch")) suppressWarnings(suppressMessages(reticulate::conda_install(packages = "pytorch-lightning")))
      if(!reticulate::py_module_available("scanpy")) suppressWarnings(suppressMessages(reticulate::conda_install(packages = "scanpy")))
    }
  }else{
    print("All python dependencies are already available!")
    }

}
