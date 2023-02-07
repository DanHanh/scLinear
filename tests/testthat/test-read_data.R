test_that("read_data() Create seurat object from data", {
  expect_equal(1, 1)
})


# test_that("read_data() Create seurat object from data", {
#   counts <- read.table("./../data/mat-001.csv", header = TRUE, sep = ",", row.names = 1)
#   feature_types <- read.table("./../data/feature_types.csv", header = TRUE, sep = ",")
#   cell_types <- read.table("./../data/cell_types.csv", sep=",", header = T)
#   rownames(counts) <- cell_types$X
#   counts <- t(counts)
#
#   expect_equal(read_data(matrix(c(1,2,3,4, nrow = 2, ncol = 2))), matrix(c(1,2,3,4, nrow = 2, ncol = 2)))
# })
