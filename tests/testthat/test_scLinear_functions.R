context("scLinear functions")


small_sobj_1 <- readRDS("./../testdata/test_object_1.rds")
small_sobj_2 <- readRDS("./../testdata/test_object_2.rds")

small_sobj_1@meta.data$sample <- "S1"
small_sobj_2@meta.data$sample <- "S2"

combined_small_sobj <- merge(small_sobj_1, y = small_sobj_2)

Seurat::DefaultAssay(small_sobj_1) <- "RNA"
Seurat::DefaultAssay(small_sobj_2) <- "RNA"
Seurat::DefaultAssay(combined_small_sobj) <- "RNA"


# Test if the remove doublets function works
test_that("remove_doublets function", {
  # Test with one sample object
  object_dedoublified <- remove_doublets(object = small_sobj_1)
  expect_equal(Seurat::DefaultAssay(small_sobj_1), "RNA")
  expect_type(object_dedoublified, "list") # returns list with plot and seurat object
  expect_s3_class(object_dedoublified$plots, "ggplot")
  expect_s4_class(object_dedoublified$object, "Seurat")
  expect_lt(ncol(object_dedoublified$object[["RNA"]]), ncol(small_sobj_1[["RNA"]])) # Test if doublets were removed

  # Test with multiple sample object
  object_combined_dedoublified <- remove_doublets(object = combined_small_sobj, samples = "sample")
  expect_equal(Seurat::DefaultAssay(combined_small_sobj), "RNA")
  expect_type(object_combined_dedoublified, "list") # returns list with plot and seurat object
  expect_s3_class(object_combined_dedoublified$plots, "ggplot")
  expect_s4_class(object_combined_dedoublified$object, "Seurat")
  expect_lt(ncol(object_combined_dedoublified$object[["RNA"]]), ncol(combined_small_sobj[["RNA"]])) # Test if doublets were removed
})




# Test if the remove doublets function works
test_that("mad_filtering function", {
  # Test with one sample object
  object_filtered <- mad_filtering(small_sobj_1)
  expect_equal(Seurat::DefaultAssay(small_sobj_1), "RNA")
  expect_equal(Seurat::DefaultAssay(object_filtered$object), "RNA")
  expect_type(object_filtered, "list") # returns list with plot and seurat object
  expect_s3_class(object_filtered$plots, "ggplot")
  expect_s4_class(object_filtered$object, "Seurat")
  expect_lt(ncol(object_filtered$object[["RNA"]]), ncol(small_sobj_1[["RNA"]])) # Test if doublets were removed

  # Test with multiple sample object
  object_combined_filtered <- mad_filtering(combined_small_sobj, samples = "sample")
  expect_equal(Seurat::DefaultAssay(combined_small_sobj), "RNA")
  expect_equal(Seurat::DefaultAssay(object_combined_filtered$ob), "RNA")
  expect_type(object_combined_filtered, "list") # returns list with plot and seurat object
  expect_s3_class(object_combined_filtered$plots, "ggplot")
  expect_s4_class(object_combined_filtered$object, "Seurat")
  expect_lt(ncol(object_combined_filtered$object[["RNA"]]), ncol(combined_small_sobj[["RNA"]])) # Test if doublets were removed
})


# Test if the remove doublets function works
test_that("integrate_samples function", {
  combined_small_sobj_integrated <- integrate_samples(combined_small_sobj, samples = "sample", npcs = 30, k.weight = 20)
  expect_true("integrated" %in% Seurat::Assays(combined_small_sobj_integrated)) # Test if integrated assay was created
  expect_equal(Seurat::DefaultAssay(combined_small_sobj_integrated), "integrated") # Test if integrated assay was set as default assay
})



# Test if the cluster_data function works
test_that("cluster_data function", {
  ## Test for single sample
  small_sobj_clustered <- cluster_data(small_sobj_1, npcs_calculate = 30)
  expect_true("seurat_clusters" %in% colnames(small_sobj_clustered@meta.data)) # Test if seurat_clusters were calculated
  ## Test for multiple samples
  combined_small_sobj_integrated <- integrate_samples(combined_small_sobj, samples = "sample", npcs = 30, k.weight = 20)
  expect_equal(Seurat::DefaultAssay(combined_small_sobj_integrated), "integrated") # Test if integrated assay was set as default assay
  combined_small_sobj_integrated_clustered <- cluster_data(combined_small_sobj_integrated, npcs_calculate = 30)
  expect_true("seurat_clusters" %in% colnames(combined_small_sobj_integrated_clustered@meta.data)) # Test if seurat_clusters were calculated
  expect_equal(combined_small_sobj_integrated_clustered@commands[["FindClusters"]]@assay.used, "integrated") # Test if integrated assay was used for clustering the data
})


# Test if the cluster_data function works
test_that("anno_celltypes function", {
  # Test for single sample
  small_sobj_clustered <- cluster_data(small_sobj_1, npcs_calculate = 30)
  small_sobj_annotated <- anno_celltypes(small_sobj_clustered, anno_level = c(1, 2, 3, 4))
  expect_true(all(c("cell_type", "cell_type_1", "cell_type_2", "cell_type_3", "cell_type_4") %in% colnames(small_sobj_annotated@meta.data))) # Test if seurat_clusters were calculated

  # Test if it works for self suplied clusters
  small_sobj_annotated <- anno_celltypes(small_sobj_clustered, anno_level = c(1, 2, 3, 4), selfClusters = Seurat::Idents(small_sobj_clustered))
  expect_true(all(c("cell_type", "cell_type_1", "cell_type_2", "cell_type_3", "cell_type_4") %in% colnames(small_sobj_annotated@meta.data))) # Test if seurat_clusters were calculated

  # Test for multiple samples
  combined_small_sobj_integrated <- integrate_samples(combined_small_sobj, samples = "sample", npcs = 30, k.weight = 20)
  combined_small_sobj_integrated_clustered <- cluster_data(combined_small_sobj_integrated, npcs_calculate = 30)
  combined_small_sobj_integrated_clustered_annotated <- anno_celltypes(combined_small_sobj_integrated_clustered, anno_level = c(1, 2, 3, 4))
  expect_true(all(c("cell_type", "cell_type_1", "cell_type_2", "cell_type_3", "cell_type_4") %in% colnames(combined_small_sobj_integrated_clustered_annotated@meta.data))) # Test if seurat_clusters were calculated
})



# Test if the adt predictor
test_that("adt predictor functions", {
  pipe <- create_adt_predictor()
  expect_type(pipe, "environment") # Test if pipe is an environment object and not NULL

  # Loading the pretrained models
  pipe_all <- load_pretrained_model(pipe = pipe, model = "all")
  expect_type(pipe_all, "environment") # Test if pipe_all is an environment object and not NULL

  pipe_bcell <- load_pretrained_model(pipe = pipe, model = "bcell")
  expect_type(pipe_bcell, "environment") # Test if pipe_bcell is an environment object and not NULL

  nkcell <- load_pretrained_model(pipe = pipe, model = "nkcell")
  expect_type(nkcell, "environment") # Test if nkcell is an environment object and not NULL

  pipe_tcell <- load_pretrained_model(pipe = pipe, model = "tcell")
  expect_type(pipe_tcell, "environment") # Test if pipe_tcell is an environment object and not NULL

  # Test if the predict function works with an seurat object
  predicted_small_sobj <- small_sobj_1
  adt_assay_predicted_all <- adt_predict(pipe_all, small_sobj_1)
  predicted_small_sobj[["predicted_ADT"]] <- adt_assay_predicted_all
  expect_s4_class(adt_assay_predicted_all, "Assay")

  # Test if prediction works for Seurat object, Assay object and matrix
  seurat_object <- small_sobj_1
  assay_object <- Seurat::GetAssay(small_sobj_1, assay = "RNA")
  dgmatrix_object <- Seurat::GetAssayData(small_sobj_1, layer = "counts")
  matrix_object <- as.matrix(Seurat::GetAssayData(small_sobj_1, layer = "counts"))

  adt_assay_predicted_all_seurat_object <- adt_predict(pipe_all, seurat_object)
  adt_assay_predicted_all_assay_object <- adt_predict(pipe_all, assay_object)
  adt_assay_predicted_all_dgmatrix_object <- adt_predict(pipe_all, dgmatrix_object)
  adt_assay_predicted_all_matrix_object <- adt_predict(pipe_all, matrix_object)

  expect_s4_class(adt_assay_predicted_all_seurat_object, "Assay")
  expect_s4_class(adt_assay_predicted_all_assay_object, "Assay")
  expect_s4_class(adt_assay_predicted_all_dgmatrix_object, "Assay")
  expect_s4_class(adt_assay_predicted_all_matrix_object, "Assay")
})
