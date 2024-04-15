


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



  ## python environment setup scaffold from github TomKellyGenetics/leiden
  if(!reticulate::py_available()){
    ## Try to install miniconda if conda is not availavle
    tryCatch({reticulate::conda_list()},
             error = function(e){
               packageStartupMessage(e)
               packageStartupMessage("Conda is not available")
               install.conda <- readline("install miniconda (yes/no)?")
               if(install.conda == "yes" || install.conda == "y"){
                 reticulate::install_miniconda()
                 reticulate::conda_update()}
             }
    )
    tryCatch({
      is.reticulate.env <- any(grepl("r-reticulate", reticulate::conda_list()$python))
      # create conda env if no base image found
      if(!(is.reticulate.env)){
        if(interactive()){
          install.deps <- readline("create conda environment (yes/no)?")
          packageStartupMessage(install.deps)
        } else {
          packageStartupMessage("create conda environment (yes/no)?")
          install.deps <- "no (use interactive mode)"
          packageStartupMessage("no (use interactive mode)")
        }
        if(install.deps == "yes" || install.deps == "y"){
          reticulate::miniconda_update()
          reticulate::conda_create(envname = "r-reticulate")
          reticulate::conda_install(envname = "r-reticulate", packages = "conda")
        }
      }
      # use "r-reticulate" or "base" image (which ever is used by reticulate if installed already)
      reticulate.env <- reticulate::conda_list()$name[grep("r-reticulate", reticulate::conda_list()$python)][1]
      packageStartupMessage(paste(c("using environment:",  reticulate.env), collapse = " "))
      suppressWarnings(suppressMessages(reticulate::use_python(reticulate::conda_python())))
      suppressWarnings(suppressMessages(reticulate::use_condaenv(reticulate.env)))
    }, error = function(e){
      packageStartupMessage("Unable to set up conda environment r-reticulate")
      packageStartupMessage("run in terminal:")
      packageStartupMessage("conda init")
      packageStartupMessage("conda create -n r-reticulate")
    },
    finally = packageStartupMessage("conda environment r-reticulate"))
  }


  tryCatch({
    is.reticulate.env <- any(grepl("r-reticulate", reticulate::conda_list()$python))
    if(reticulate::py_available() || is.reticulate.env ){

      ## test if any python dependency is not available
      if(!all(
        reticulate::py_module_available("numpy"),
        reticulate::py_module_available("joblib"),
        reticulate::py_module_available("pytorch_lightning"),
        reticulate::py_module_available("torch"),
        reticulate::py_module_available("sklearn"),
        reticulate::py_module_available("anndata"),
        reticulate::py_module_available("scanpy"),
        reticulate::py_module_available("os"),
        reticulate::py_module_available("warnings"),
        reticulate::py_module_available("scipy"),
        reticulate::py_module_available("typing")
      )){


        if(interactive()){
          install.deps <- readline("install dependencies (yes/no)?")
          packageStartupMessage(install.deps)
        } else {
          packageStartupMessage("install dependencies (yes/no)?")
          install.deps <- "no (use interactive mode)"
          packageStartupMessage("no (use interactive mode)")
        }
        if(install.deps == "yes" || install.deps == "y"){
          reticulate.env <- reticulate::conda_list()$name[grep("r-reticulate", reticulate::conda_list()$python)][1]
          packageStartupMessage(paste(c("using environment:",  reticulate.env), collapse = " "))
          install_python_modules <- function(method = "auto", conda = "auto") {
            if(!is.null(reticulate::conda_binary())){
              reticulate::use_python(reticulate::conda_python())
              if(!(is.reticulate.env)){
                reticulate::conda_create(envname = reticulate.env)
                if(!reticulate::py_module_available("conda")) reticulate::conda_install(envname = reticulate.env, packages = "conda")
              }
              suppressWarnings(suppressMessages(reticulate::use_condaenv(reticulate.env)))
              if(.Platform$OS.type == "windows"){
                utils::install.packages("devtools",  quiet = TRUE)
                devtools::install_github("rstudio/reticulate", ref = "86ebb56",  quiet = TRUE)
                if(!reticulate::py_module_available("numpy")) ((reticulate::conda_install(envname = reticulate.env, packages = "numpy")))
                if(!reticulate::py_module_available("joblib")) ((reticulate::conda_install(envname = reticulate.env, packages = "joblib")))
                if(!reticulate::py_module_available("pytorch_lightning")) ((reticulate::conda_install(envname = reticulate.env, packages = "pytorch-lightning")))
                if(!reticulate::py_module_available("torch")) ((reticulate::conda_install(envname = reticulate.env, packages = "pytorch-lightning")))
                if(!reticulate::py_module_available("sklearn")) ((reticulate::conda_install(envname = reticulate.env, packages = "scikit-learn")))
                if(!reticulate::py_module_available("anndata")) ((reticulate::conda_install(envname = reticulate.env, packages = "anndata")))
                if(!reticulate::py_module_available("scanpy")) ((reticulate::conda_install(envname = reticulate.env, packages = "scanpy")))
                if(!reticulate::py_module_available("os")) ((reticulate::conda_install(envname = reticulate.env, packages = "os")))
                if(!reticulate::py_module_available("warnings")) ((reticulate::conda_install(envname = reticulate.env, packages = "warnings")))
                if(!reticulate::py_module_available("scipy")) ((reticulate::conda_install(envname = reticulate.env, packages = "scipy")))
                if(!reticulate::py_module_available("typing")) ((reticulate::conda_install(envname = reticulate.env, packages = "typing")))
                utils::install.packages("reticulate",  quiet = TRUE)
              } else {
                if(!reticulate::py_module_available("numpy")) ((reticulate::conda_install(reticulate.env, "numpy")))
                if(!reticulate::py_module_available("joblib")) ((reticulate::conda_install(reticulate.env, "joblib")))
                if(!reticulate::py_module_available("pytorch_lightning")) ((reticulate::conda_install(reticulate.env, "pytorch-lightning")))
                if(!reticulate::py_module_available("torch")) ((reticulate::conda_install(reticulate.env, "pytorch-lightning")))
                if(!reticulate::py_module_available("sklearn")) ((reticulate::conda_install(reticulate.env, "scikit-learn")))
                if(!reticulate::py_module_available("anndata")) ((reticulate::conda_install(reticulate.env, "anndata")))
                if(!reticulate::py_module_available("scanpy")) ((reticulate::conda_install(reticulate.env, "scanpy")))
                if(!reticulate::py_module_available("os")) ((reticulate::conda_install(reticulate.env, "os")))
                if(!reticulate::py_module_available("warnings")) ((reticulate::conda_install(reticulate.env, "warnings")))
                if(!reticulate::py_module_available("scipy")) ((reticulate::conda_install(reticulate.env, "scipy")))
                if(!reticulate::py_module_available("typing")) ((reticulate::conda_install(reticulate.env, "typing")))
                Sys.setenv(RETICULATE_PYTHON = reticulate::conda_python())
              }
            } else {

              if(!reticulate::py_module_available("numpy")) ((reticulate::py_install(reticulate.env, "numpy")))
              if(!reticulate::py_module_available("joblib")) ((reticulate::py_install(reticulate.env, "joblib")))
              if(!reticulate::py_module_available("pytorch_lightning")) ((reticulate::py_install(reticulate.env, "pytorch-lightning")))
              if(!reticulate::py_module_available("torch")) ((reticulate::py_install(reticulate.env, "pytorch-lightning")))
              if(!reticulate::py_module_available("sklearn")) ((reticulate::py_install(reticulate.env, "scikit-learn")))
              if(!reticulate::py_module_available("anndata")) ((reticulate::py_install(reticulate.env, "anndata")))
              if(!reticulate::py_module_available("scanpy")) ((reticulate::py_install(reticulate.env, "scanpy")))
              if(!reticulate::py_module_available("os")) ((reticulate::py_install(reticulate.env, "os")))
              if(!reticulate::py_module_available("warnings")) ((reticulate::py_install(reticulate.env, "warnings")))
              if(!reticulate::py_module_available("scipy")) ((reticulate::py_install(reticulate.env, "scipy")))
              if(!reticulate::py_module_available("typing")) ((reticulate::py_install(reticulate.env, "typing")))
              Sys.setenv(RETICULATE_PYTHON = reticulate::py_config()$python)
            }
          }
        }
        quiet <- function(expr, all = TRUE) {
          if (Sys.info()['sysname'] == "Windows") {
            file <- "NUL"
          } else {
            file <- "/dev/null"
          }

          if (all) {
            suppressWarnings(suppressMessages(suppressPackageStartupMessages(
              capture.output(expr, file = file)
            )))
          } else {
            capture.output(expr, file = file)
          }
        }
        install_python_modules()
      }
    }
  }, error = function(e){
    packageStartupMessage("Unable to install python modules")
    packageStartupMessage("run in terminal:")
    packageStartupMessage("conda install -n r-reticulate -c conda-forge <pkg_name>")
  },
  finally = packageStartupMessage("all python modules installed"))


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
