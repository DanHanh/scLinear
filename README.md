
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scLinear

<!-- badges: start -->
<!-- badges: end -->

The goal of scLinear is to …

## Installation

You can install the development version of scLinear from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DanHanh/scLinear")
```

## Get example data

This is a basic example which shows you how to solve a common problem:

``` r
library(scLinear)

## get some example data from the Seurat PBMC10K example data set for the tutorial: https://satijalab.org/seurat/articles/multimodal_vignette.html.
# The download link: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3.
# The File: "Feature / cell matrix (raw)"

pbmc10k.data <- Seurat::Read10X(data.dir = "./local/raw_feature_bc_matrix")
#> 10X data contains more than one type and is being returned as a list containing matrices of each type.
rownames(x = pbmc10k.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = pbmc10k.data[["Antibody Capture"]]))
pbmc10k <- Seurat::CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 1, min.features = 1)
pbmc10k[["ADT"]] <- Seurat::CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])
Seurat::DefaultAssay(pbmc10k) <- "RNA"
```

## Prepare data

``` r
pbmc10k <- prepare_data(pbmc10k, integrate_data = FALSE, annotation_selfCluster = TRUE, remove_empty_droplets = TRUE)
#> Centering and scaling data matrix
#> [1] "Start remove doublets"
#> Clustering cells...
#> 11 clusters
#> Creating ~5054 artifical doublets...
#> Dimensional reduction
#> Evaluating kNN...
#> Training model...
#> iter=0, 497 cells excluded from training.
#> Threshold found:0.237
#> 702 (8.3%) doublets called
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

    #> [1] "Start low quality cell removal"

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

    #> [1] "Start clustering data"
    #> Centering and scaling data matrix
    #> Warning in PrepDR(object = object, features = features, verbose = verbose): The
    #> following 50 features requested have zero variance (running reduction without
    #> them): PBK, SDC1, IGHV3-23, IGLV1-44, MCM10, HJURP, E2F8, IGLV2-8, SPACA3,
    #> IGHV4-4, IGKV1D-39, IGHV6-1, IGKV1-8, IGKV1D-12, CXCL13, CCNA1, IGLVI-70, CA1,
    #> LCN6, IGKV2D-28, SLC6A9, IGHGP, DPEP1, SPATA3, FAM19A5, TLX2, SAMD7, SHISA3,
    #> PF4V1, MAMDC2-AS1, HASPIN, CDC25C, AC106897.1, TEX49, IGHV3-49, AC012073.1,
    #> TM4SF1, AC112722.1, TCIM, AL353795.2, AMHR2, IGHV3-66, IGHV1-69D, LRRC49, PPY,
    #> SLC24A3, GYPA, TPSAB1, SORBS2, CRYAB
    #> PC_ 1 
    #> Positive:  LYZ, FCN1, CST3, MNDA, S100A9, S100A8, SERPINA1, VCAN, CD14, CSTA 
    #>     TYMP, SPI1, GRN, CTSS, CFD, S100A12, MS4A6A, CYBB, AIF1, LST1 
    #>     CD68, CLEC7A, MPEG1, FGL2, FCER1G, KLF4, CLEC12A, CEBPD, FTL, FTH1 
    #> Negative:  IL32, TRBC2, TRAC, LTB, IL7R, LDHB, NPM1, HSPA8, ISG20, TRBC1 
    #>     CTSW, CD27, ITM2A, SPOCK2, LEF1, CCL5, GZMA, RORA, CST7, HMGN1 
    #>     C12orf75, NKG7, PRF1, KLRB1, AQP3, PTMA, CDC25B, TRAT1, MT-CYB, KLRD1 
    #> PC_ 2 
    #> Positive:  CD79A, MS4A1, BANK1, IGHM, LINC00926, IGHD, CD79B, HLA-DQA1, HLA-DQB1, TCL1A 
    #>     FCER2, TNFRSF13C, BCL11A, CD22, RALGPS2, SPIB, VPREB3, TSPAN13, ADAM28, P2RX5 
    #>     FCRLA, FCRL1, FAM129C, HLA-DOB, IGKC, IGLC2, HLA-DRA, CD19, EAF2, MEF2C 
    #> Negative:  NKG7, S100A4, GZMA, CTSW, CST7, CCL5, IL32, PRF1, GNLY, KLRD1 
    #>     PFN1, GZMH, FGFBP2, KLRF1, FCGR3A, ACTG1, GZMB, SRGN, HOPX, SPON2 
    #>     S100A6, CMC1, ACTB, MATK, TRDC, CCL4, CALM1, MYOM2, C12orf75, PRSS23 
    #> PC_ 3 
    #> Positive:  NKG7, GZMB, GZMH, FGFBP2, GNLY, KLRD1, CST7, PRF1, FCGR3A, KLRF1 
    #>     GZMA, CCL5, SPON2, CTSW, HOPX, CCL4, HLA-DPB1, TRDC, HLA-DPA1, PRSS23 
    #>     CMC1, CD74, CLIC3, C12orf75, MATK, MYOM2, RHOC, HLA-DRB1, ADGRG1, HLA-DQA2 
    #> Negative:  IL7R, LEF1, LDHB, TRAC, CD27, LTB, TRAT1, RGCC, NPM1, CD40LG 
    #>     AQP3, NELL2, CD5, RPLP0, LRRN3, SPOCK2, ITM2A, ANP32B, TRBC2, CD28 
    #>     SOCS3, PASK, MYC, ICOS, H1FX, SLC2A3, TSHZ2, SIRPG, CORO1B, CRIP2 
    #> PC_ 4 
    #> Positive:  FGFBP2, GZMH, KLRD1, KLRF1, GNLY, CYP1B1, SPON2, VNN2, S100A12, TRDC 
    #>     PRSS23, PRF1, RBP7, MYOM2, ADGRG1, CCL4, G0S2, PADI4, TKTL1, S100A8 
    #>     NKG7, CST7, FCRL6, MARC1, CXCL8, CDA, TYROBP, GZMB, ANPEP, MT-CO3 
    #> Negative:  PLD4, GAS6, SERPINF1, IL3RA, ITM2C, YBX1, PPP1R14B, LILRA4, RPLP0, FCER1A 
    #>     CCDC50, LDHB, CLEC4C, AL096865.1, CLEC10A, GPR183, DERL3, SCT, CORO1B, DNASE1L3 
    #>     AQP3, UGCG, NPM1, LRRC26, LDHA, SMPD3, TRAC, SLC25A5, CDKN1C, LMNA 
    #> PC_ 5 
    #> Positive:  CDKN1C, HES4, TCF7L2, C1QA, SIGLEC10, CTSL, IFITM3, AC064805.1, CKB, CSF1R 
    #>     MS4A7, BATF3, NEURL1, HMOX1, FCGR3A, MTSS1, ZNF703, RRAS, RHOC, CAMK1 
    #>     ICAM4, ABI3, LYPD2, MARCKS, LRRC25, IFI30, SMIM25, GPBAR1, MS4A4A, FMNL2 
    #> Negative:  LILRA4, SERPINF1, DERL3, CLEC4C, SCT, DNASE1L3, LRRC26, PTCRA, IL3RA, ITM2C 
    #>     SMPD3, MZB1, GAS6, RHEX, SCAMP5, JCHAIN, SMIM5, TPM2, PLD4, ZFAT 
    #>     AL096865.1, LAMP5, UGCG, CCDC50, PTGDS, TNFRSF21, APP, TCF4, IRF4, MAP1A
    #> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    #> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    #> This message will be shown once per session
    #> 13:40:51 UMAP embedding parameters a = 0.9922 b = 1.112
    #> 13:40:51 Read 6792 rows and found 30 numeric columns
    #> 13:40:51 Using Annoy for neighbor search, n_neighbors = 30
    #> 13:40:51 Building Annoy index with metric = cosine, n_trees = 50
    #> 0%   10   20   30   40   50   60   70   80   90   100%
    #> [----|----|----|----|----|----|----|----|----|----|
    #> **************************************************|
    #> 13:40:52 Writing NN index file to temp file /tmp/Rtmpm4iDmA/file1f98228f8aa43
    #> 13:40:52 Searching Annoy index using 1 thread, search_k = 3000
    #> 13:40:53 Annoy recall = 100%
    #> 13:40:54 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    #> 13:40:56 Initializing from normalized Laplacian + noise (using irlba)
    #> 13:40:57 Commencing optimization for 500 epochs, with 296664 positive edges
    #> 13:41:03 Optimization finished
    #> Computing nearest neighbor graph
    #> Computing SNN
    #> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    #> 
    #> Number of nodes: 6792
    #> Number of edges: 302688
    #> 
    #> Running Louvain algorithm...
    #> Maximum modularity in 10 random starts: 0.8736
    #> Number of communities: 15
    #> Elapsed time: 0 seconds
    #> Warning in PrepDR(object = object, features = features, verbose = verbose): The
    #> following 50 features requested have zero variance (running reduction without
    #> them): PBK, SDC1, IGHV3-23, IGLV1-44, MCM10, HJURP, E2F8, IGLV2-8, SPACA3,
    #> IGHV4-4, IGKV1D-39, IGHV6-1, IGKV1-8, IGKV1D-12, CXCL13, CCNA1, IGLVI-70, CA1,
    #> LCN6, IGKV2D-28, SLC6A9, IGHGP, DPEP1, SPATA3, FAM19A5, TLX2, SAMD7, SHISA3,
    #> PF4V1, MAMDC2-AS1, HASPIN, CDC25C, AC106897.1, TEX49, IGHV3-49, AC012073.1,
    #> TM4SF1, AC112722.1, TCIM, AL353795.2, AMHR2, IGHV3-66, IGHV1-69D, LRRC49, PPY,
    #> SLC24A3, GYPA, TPSAB1, SORBS2, CRYAB
    #> PC_ 1 
    #> Positive:  LYZ, FCN1, CST3, MNDA, S100A9, S100A8, SERPINA1, VCAN, CD14, CSTA 
    #>     TYMP, SPI1, GRN, CTSS, CFD, S100A12, MS4A6A, CYBB, AIF1, LST1 
    #>     CD68, CLEC7A, MPEG1, FGL2, FCER1G, KLF4, CLEC12A, CEBPD, FTL, FTH1 
    #> Negative:  IL32, TRBC2, TRAC, LTB, IL7R, LDHB, NPM1, HSPA8, ISG20, TRBC1 
    #>     CTSW, CD27, ITM2A, SPOCK2, LEF1, CCL5, GZMA, RORA, CST7, HMGN1 
    #>     C12orf75, NKG7, PRF1, KLRB1, AQP3, PTMA, CDC25B, TRAT1, MT-CYB, KLRD1 
    #> PC_ 2 
    #> Positive:  CD79A, MS4A1, BANK1, IGHM, LINC00926, IGHD, CD79B, HLA-DQA1, HLA-DQB1, TCL1A 
    #>     FCER2, TNFRSF13C, BCL11A, CD22, RALGPS2, SPIB, VPREB3, TSPAN13, ADAM28, P2RX5 
    #>     FCRLA, FCRL1, FAM129C, HLA-DOB, IGKC, IGLC2, HLA-DRA, CD19, EAF2, MEF2C 
    #> Negative:  NKG7, S100A4, GZMA, CTSW, CST7, CCL5, IL32, PRF1, GNLY, KLRD1 
    #>     PFN1, GZMH, FGFBP2, KLRF1, FCGR3A, ACTG1, GZMB, SRGN, HOPX, SPON2 
    #>     S100A6, CMC1, ACTB, MATK, TRDC, CCL4, CALM1, MYOM2, C12orf75, PRSS23 
    #> PC_ 3 
    #> Positive:  NKG7, GZMB, GZMH, FGFBP2, GNLY, KLRD1, CST7, PRF1, FCGR3A, KLRF1 
    #>     GZMA, CCL5, SPON2, CTSW, HOPX, CCL4, HLA-DPB1, TRDC, HLA-DPA1, PRSS23 
    #>     CMC1, CD74, CLIC3, C12orf75, MATK, MYOM2, RHOC, HLA-DRB1, ADGRG1, HLA-DQA2 
    #> Negative:  IL7R, LEF1, LDHB, TRAC, CD27, LTB, TRAT1, RGCC, NPM1, CD40LG 
    #>     AQP3, NELL2, CD5, RPLP0, LRRN3, SPOCK2, ITM2A, ANP32B, TRBC2, CD28 
    #>     SOCS3, PASK, MYC, ICOS, H1FX, SLC2A3, TSHZ2, SIRPG, CORO1B, CRIP2 
    #> PC_ 4 
    #> Positive:  FGFBP2, GZMH, KLRD1, KLRF1, GNLY, CYP1B1, SPON2, VNN2, S100A12, TRDC 
    #>     PRSS23, PRF1, RBP7, MYOM2, ADGRG1, CCL4, G0S2, PADI4, TKTL1, S100A8 
    #>     NKG7, CST7, FCRL6, MARC1, CXCL8, CDA, TYROBP, GZMB, ANPEP, MT-CO3 
    #> Negative:  PLD4, GAS6, SERPINF1, IL3RA, ITM2C, YBX1, PPP1R14B, LILRA4, RPLP0, FCER1A 
    #>     CCDC50, LDHB, CLEC4C, AL096865.1, CLEC10A, GPR183, DERL3, SCT, CORO1B, DNASE1L3 
    #>     AQP3, UGCG, NPM1, LRRC26, LDHA, SMPD3, TRAC, SLC25A5, CDKN1C, LMNA 
    #> PC_ 5 
    #> Positive:  CDKN1C, HES4, TCF7L2, C1QA, SIGLEC10, CTSL, IFITM3, AC064805.1, CKB, CSF1R 
    #>     MS4A7, BATF3, NEURL1, HMOX1, FCGR3A, MTSS1, ZNF703, RRAS, RHOC, CAMK1 
    #>     ICAM4, ABI3, LYPD2, MARCKS, LRRC25, IFI30, SMIM25, GPBAR1, MS4A4A, FMNL2 
    #> Negative:  LILRA4, SERPINF1, DERL3, CLEC4C, SCT, DNASE1L3, LRRC26, PTCRA, IL3RA, ITM2C 
    #>     SMPD3, MZB1, GAS6, RHEX, SCAMP5, JCHAIN, SMIM5, TPM2, PLD4, ZFAT 
    #>     AL096865.1, LAMP5, UGCG, CCDC50, PTGDS, TNFRSF21, APP, TCF4, IRF4, MAP1A 
    #> 13:41:35 UMAP embedding parameters a = 0.9922 b = 1.112
    #> 13:41:35 Read 6792 rows and found 30 numeric columns
    #> 13:41:35 Using Annoy for neighbor search, n_neighbors = 30
    #> 13:41:35 Building Annoy index with metric = cosine, n_trees = 50
    #> 0%   10   20   30   40   50   60   70   80   90   100%
    #> [----|----|----|----|----|----|----|----|----|----|
    #> **************************************************|
    #> 13:41:35 Writing NN index file to temp file /tmp/Rtmpm4iDmA/file1f9826a8afe87
    #> 13:41:35 Searching Annoy index using 1 thread, search_k = 3000
    #> 13:41:37 Annoy recall = 100%
    #> 13:41:38 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    #> 13:41:40 Initializing from normalized Laplacian + noise (using irlba)
    #> 13:41:40 Commencing optimization for 500 epochs, with 296664 positive edges
    #> 13:41:47 Optimization finished

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" />

    #> [1] "Start cell type annotation"
    #> Warning: replacing previous import 'Matrix::cov2cor' by 'stats::cov2cor' when
    #> loading 'scMRMA'
    #> Warning: replacing previous import 'Matrix::toeplitz' by 'stats::toeplitz' when
    #> loading 'scMRMA'
    #> Warning: replacing previous import 'Matrix::update' by 'stats::update' when
    #> loading 'scMRMA'
    #> Warning: replacing previous import 'Matrix::expand' by 'tidyr::expand' when
    #> loading 'scMRMA'
    #> Warning: replacing previous import 'Matrix::pack' by 'tidyr::pack' when loading
    #> 'scMRMA'
    #> Warning: replacing previous import 'Matrix::unpack' by 'tidyr::unpack' when
    #> loading 'scMRMA'
    #> Loading required package: Seurat
    #> Attaching SeuratObject

<img src="man/figures/README-unnamed-chunk-3-4.png" width="100%" />

    #> Pre-defined cell type database panglaodb will be used.
    #> Multi Resolution Annotation Started. 
    #> Level 1 annotation started. 
    #> Level 2 annotation started. 
    #> Level 3 annotation started. 
    #> Level 4 annotation started. 
    #> Uniform Resolution Annotation Started.

<img src="man/figures/README-unnamed-chunk-3-5.png" width="100%" />

## Train new model

``` r
## create a training and test set
set.seed(42)
indx <- sample(1:length(colnames(pbmc10k)), size = length(colnames(pbmc10k)), replace = FALSE)
pbmc10k_train <- pbmc10k[,indx[1:5000]]
pbmc10k_test <- pbmc10k[,indx[5001:length(colnames(pbmc10k))]]

## create predictor
pipe <- create_adt_predictor()

## train predictor
pipe <- fit_predictor(pipe = pipe, gexp_train = pbmc10k_train@assays[["RNA"]],
              adt_train = pbmc10k_train@assays[["ADT"]],
              normalize = TRUE)
#> Normalizing across cells

## evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                  gexp_test = pbmc10k_test@assays[["RNA"]],
                  adt_test = pbmc10k_test@assays[["ADT"]],
                  normalize = TRUE)
#> Normalizing across cells

## add predicted adt assay
pbmc10k_test@assays[["predicted_ADT"]] <-  adt_predict(pipe = pipe, gexp = pbmc10k_test@assays[["RNA"]], normalize = TRUE)
```

## Use pretrained model

``` r
library(scLinear)

## subset pbmc data to only T-cells
b_cells <- pbmc10k %>% base::subset(subset = cell_type == "T")

pipe <- create_adt_predictor()
pipe$gex_preprocessor$do_log1p <- FALSE
## load pre-trained model (available models: all, bcell, tcell, nkcell)
pipe <- load_pretrained_model(pipe, model = "all")

pipe$gex_preprocessor$do_log1p <- FALSE
evaluate_predictor(pipe, b_cells@assays$RNA, b_cells@assays$ADT, normalize = TRUE)
#> Normalizing across cells
#> [[1]]
#> [1] 0.6037294
#> 
#> [[2]]
#> [1] 0.8506006
#> 
#> [[3]]
#> [1] 0.7471668
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
