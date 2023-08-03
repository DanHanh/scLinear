
numpy <- NULL
scanpy <- NULL
joblib <- NULL
evaluate <- NULL
prediction <- NULL
preprocessing <- NULL
torch <- NULL
pytorch_ligthning <- NULL

.onLoad <- function(libname, pkgname){
    reticulate::configure_environment(pkgname)
    module_path <-  base::system.file("python",package = "scLinear")
    scanpy <<- reticulate::import("scanpy", delay_load = TRUE)
    numpy <<- reticulate::import("numpy", delay_load = TRUE)
    joblib <<- reticulate::import("joblib", delay_load = TRUE)
    preprocessing <<- reticulate::import_from_path(module = "preprocessing", path = module_path, delay_load = TRUE)
    prediction <<- reticulate::import_from_path(module = "prediction", path = module_path, delay_load = TRUE)
    evaluate <<- reticulate::import_from_path(module = "evaluate", path = module_path, delay_load = TRUE)
    torch <<- reticulate::import_from_path(module = "torch", path = module_path, delay_load = TRUE)
    pytorch_ligthning <<- reticulate::import_from_path(module = "pytorch-lightning", path = module_path, delay_load = TRUE)
    sklearn <<- reticulate::import_from_path(module = "scikit-learn", path = module_path, delay_load = TRUE)

    }
