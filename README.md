
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scLinear

<!-- badges: start -->
<!-- badges: end -->

The goal of scLinear is to predict antibody derived tags (ADT) data from
gene expression data in scRNA-seq data. It includes all the necessary pre-processing steps, comes equiped with pre-trained models and also allows the training of new models.  
<p align="center"><img src="man/figures/schematic.v5.2.png" width="75%" height="75%" align="middle"/><p>  
  
- [Installation](#Installation)
- [Example](#Example)
- [Other Functions](#Other-functions)
- [Citation](#Citation)
- [Contact us](#Contact-us)
## Installation

You can install the development version of scLinear using devtools.  

``` r
# install.packages("devtools")
devtools::install_github("DanHanh/scLinear")
```
## Example
After installing scLinear you may run the following example.

### Get data
The PBMC data can be downloaded from the 10X Genomics website (https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3). Then the following code can be used to generate a Seurat object.

``` r
set.seed(42)

library(scLinear)

## Example data from the Seurat PBMC10K example data set for the tutorial: https://satijalab.org/seurat/articles/multimodal_vignette.html.
# Download link: "https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3".
# File: "Feature / cell matrix (filtered)"

pbmc10k.data <- Seurat::Read10X(data.dir = "./local/raw_feature_bc_matrix/pbmc_10k_protein_v3_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
rownames(x = pbmc10k.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = pbmc10k.data[["Antibody Capture"]]))
pbmc10k <- Seurat::CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 1, min.features = 1)
pbmc10k[["ADT"]] <- Seurat::CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])
Seurat::DefaultAssay(pbmc10k) <- "RNA"
```

### Running scLinear
You may run scLinear directly providing the Seurat object as input. List of optional parameters:  
  
* `remove_doublets` Removal of doublets. TRUE (default) or FALSE.  
* `low_qc_cell_removal` Removal of low quality cells. TRUE (default) or FALSE.  
* `anno_level` Level of annotation. 1, 2, 3 or 4. See https://github.com/JiaLiVUMC/scMRMA for more details.  
* `samples` NULL (default).  
* `integrate_data` TRUE or FALSE (default).  
* `remove_empty_droplets` Removal of empty droplets. TRUE or FALSE (default).  
* `lower` = 100 (default).  
* `FDR` = 0.01 (default).  
* `annotation_selfCluster` TRUE or FALSE (default).  
* `resolution` = 0.8 (default).  
* `seed` = 42 (default).  
* `return_plots` TRUE or FALSE (default).  
* `model` Available models "all" (default), "bcell", "tcell" and "nkcell".  
* `assay_name` = "RNA" (default).  
* `print_plots` TRUE or FALSE (default).  

The scLinear function uses the counts slot from the RNA assay to predict the ADT assay. The functions performs the default preprocessing steps and returns a Seurat object with the added "predicted_ADT" assay.
  
``` r
pbmc10k_adt_predicted <- scLinear(pbmc10k)
#> [1] "Start remove doublets"
#> [1] "Start low quality cell removal"
#> [1] "Start clustering data"
#> [1] "Number of used dimensions for clustering: 26"
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 6703
#> Number of edges: 285702
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.8810
#> Number of communities: 14
#> Elapsed time: 0 seconds
#> [1] "Start cell type annotation"
#> Pre-defined cell type database panglaodb will be used.
#> Multi Resolution Annotation Started. 
#> Level 1 annotation started. 
#> Level 2 annotation started. 
#> Level 3 annotation started. 
#> Level 4 annotation started. 
#> Uniform Resolution Annotation Started.
```
## Other functions
scLinear calls different sub-workflows which can also be called independentaly.
### Prepare data
`prepare_data()` performes all the necessay pre-processing steps. Parameters are the same as [above](##Running-scLinear") 

``` r
pbmc10k <- prepare_data(pbmc10k,
        integrate_data = FALSE,
        annotation_selfCluster = TRUE, 
        remove_empty_droplets = FALSE)
#> [1] "Start remove doublets"
```

<p align="center"><img src="man/figures/README-unnamed-chunk-4-1.png" width="75%" height="75%"/><p>

    #> [1] "Start low quality cell removal"

<p align="center"><img src="man/figures/README-unnamed-chunk-4-2.png" width="75%" height="75%"/><img src="man/figures/README-unnamed-chunk-4-3.png" width="75%" height="75%"/><p>

    #> [1] "Start clustering data"
    #> [1] "Number of used dimensions for clustering: 26"
    #> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    #> 
    #> Number of nodes: 6703
    #> Number of edges: 285702
    #> 
    #> Running Louvain algorithm...
    #> Maximum modularity in 10 random starts: 0.8810
    #> Number of communities: 14
    #> Elapsed time: 0 seconds
    #> [1] "Start cell type annotation"
    #> Pre-defined cell type database panglaodb will be used.
    #> Multi Resolution Annotation Started. 
    #> Level 1 annotation started. 
    #> Level 2 annotation started. 
    #> Level 3 annotation started. 
    #> Level 4 annotation started. 
    #> Uniform Resolution Annotation Started.

<p align="center"><img src="man/figures/README-unnamed-chunk-4-4.png" width="75%" height="75%"/><p>

### Use a pre-trained model
User may manually load pre-trained models (available models: all, bcell, tcell, nkcell). If a pretrained model is used it is advided to use the raw data slot from the RNA assay, and normalization = TRUE, to ensure that the input data is normalized the same way as for the training data. `adt_predict()` is then used to predict the ADT values with parameters:  
* `pipe` Pretrained model.
* `gexp` Gene expression matrix.
* `slot` Seurat slot to use. "counts" (default)
* `normalize` TRUE (default) or FALSE.
  An example can be found below:

``` r
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, model = "all")

pbmc10k@assays["predicted_ADT"] <-  adt_predict(pipe = pipe,
                        gexp = pbmc10k@assays[["RNA"]],
                        normalize = TRUE)
```

## Train a new model

To train a new model the following commands need to be used.   
`create_adt_predictor()` to initialize predictor.  
`fit_predictor()` with parameters:
* `pipe` predictor initialized above.  
* `gexp_train` gene expression matrix of training set (i.e. RNA assay from Seurat object).  
* `adt_train` ADT matrix of training set.   
* `normalize_gex` Gene expression normalization. TRUE (default) or FALSE.
* `normalize_adt` ADT normalization. TRUE (default) or FALSE.
  
Subsequently, the `evaluate_predictor` command can be used (same parameters with `fit_predictor`) to return the RMSE, Pearson and Spearman of the training process. An example of this process can be found below.  

``` r
## Create a training and a test set
indx <- sample(1:length(colnames(pbmc10k)), size = length(colnames(pbmc10k)), replace = FALSE)
pbmc10k_train <- pbmc10k[,indx[1:5000]]
pbmc10k_test <- pbmc10k[,indx[5001:length(colnames(pbmc10k))]]

## create predictor
pipe <- create_adt_predictor()

## train predictor
pipe <- fit_predictor(pipe = pipe,
 gexp_train = pbmc10k_train@assays[["RNA"]],
              adt_train = pbmc10k_train@assays[["ADT"]],
              normalize_gex = TRUE,
              normalize_adt = TRUE)

## evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                  gexp_test = pbmc10k_test@assays[["RNA"]],
                  adt_test = pbmc10k_test@assays[["ADT"]],
                  normalize_gex = TRUE,
                  normalize_adt = TRUE)

print(eval_res)
#>        RMSE   Pearson  Spearman
#> 1 0.3473135 0.9436867 0.8726635

## add the predicted adt assay
pbmc10k_test@assays["predicted_ADT"] <-  adt_predict(pipe = pipe,
                        gexp = pbmc10k_test@assays[["RNA"]],
                        normalize = TRUE)
```
## Citation
Daniel Hanhart et al., "scLinear; against the deep current to predict protein abundance at single-cell resolution"

## Contact us
For any request or question you may contact:
* Daniel Hanhart daniel.hanhart@unibe.ch
* Panagiotis Chouvardas panagiotis.chouvardas@unibe.ch
