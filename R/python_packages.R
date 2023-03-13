
numpy <- NULL
scanpy <- NULL
joblib <- NULL
evaluate <- NULL
prediction <- NULL
preprocessing <- NULL
.onLoad <- function(libname, pkgname){
    numpy <<- reticulate::import("numpy", delay_load = TRUE)
    scanpy <<- reticulate::import("scanpy", delay_load = TRUE)
    joblib <<- reticulate::import("joblib", delay_load = TRUE)
    evaluate <<- reticulate::import_from_path(module = "evaluate", path = "./inst/python/", delay_load = TRUE)
    prediction <<- reticulate::import_from_path(module = "prediction", path = "./inst/python/", delay_load = TRUE)
    preprocessing <<- reticulate::import_from_path(module = "preprocessing", path = "./inst/python/", delay_load = TRUE)
}
