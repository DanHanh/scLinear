
numpy <- NULL
scanpy <- NULL
joblib <- NULL
evaluate <- NULL
prediction <- NULL
preprocessing <- NULL

.onLoad <- function(libname, pkgname){
    module_path <-  base::system.file("python",package = "scLinearDev")
    numpy <<- reticulate::import("numpy", delay_load = TRUE)
    scanpy <<- reticulate::import("scanpy", delay_load = TRUE)
    joblib <<- reticulate::import("joblib", delay_load = TRUE)
    preprocessing <<- reticulate::import_from_path(module = "preprocessing", path = module_path, delay_load = TRUE)
    prediction <<- reticulate::import_from_path(module = "prediction", path = module_path, delay_load = TRUE)
    evaluate <<- reticulate::import_from_path(module = "evaluate", path = module_path, delay_load = TRUE)
    # reticulate::import_from_path(module = "preprocessing", path = "inst/python", delay_load = TRUE)
    # reticulate::import_from_path(module = "prediction", path = "inst/python", delay_load = TRUE)
    # reticulate::import_from_path(module = "evaluate", path = "inst/python", delay_load = TRUE)

    # reticulate::source_python("./inst/python/preprocessing.py")
    # reticulate::source_python("./inst/python/prediction.py")
    # reticulate::source_python("./inst/python/evaluate.py")
}
