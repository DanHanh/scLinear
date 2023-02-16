
numpy <- NULL
scanpy <- NULL
.onLoad <- function(libname, pkgname){
    numpy <<- reticulate::import("numpy", delay_load = TRUE)
    scanpy <<- reticulate::import("scanpy", delay_load = TRUE)
}
