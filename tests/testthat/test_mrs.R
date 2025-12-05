# Test mrs() function

test_that("mrs returns correct structure", {
  pdb_file <- system.file("extdata", "2acy.pdb", package = "mutenm")
  pdb <- bio3d::read.pdb(pdb_file)

  result <- mrs(pdb, nmut = 1, seed = 42)

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("dr2ij", "influence", "sensitivity", "params"))

  # Check dimensions
  nsites <- result$params$nsites
  expect_equal(dim(result$dr2ij), c(nsites, nsites))
  expect_length(result$influence, nsites)
  expect_length(result$sensitivity, nsites)

  # Check that influence and sensitivity are column/row means
  expect_equal(result$influence, colMeans(result$dr2ij))
  expect_equal(result$sensitivity, rowMeans(result$dr2ij))
})

test_that("mrs with different parameters works", {
  pdb_file <- system.file("extdata", "2acy.pdb", package = "mutenm")
  pdb <- bio3d::read.pdb(pdb_file)

  # Test with different node type
  result_cb <- mrs(pdb, node = "cb", nmut = 1, seed = 42)
  expect_type(result_cb, "list")
  expect_equal(result_cb$params$node, "cb")

  # Test with different model
  result_mw <- mrs(pdb, model = "ming_wall", nmut = 1, seed = 42)
  expect_type(result_mw, "list")
  expect_equal(result_mw$params$model, "ming_wall")
})

test_that("mrs is reproducible with seed", {
  pdb_file <- system.file("extdata", "2acy.pdb", package = "mutenm")
  pdb <- bio3d::read.pdb(pdb_file)

  result1 <- mrs(pdb, nmut = 1, seed = 123)
  result2 <- mrs(pdb, nmut = 1, seed = 123)

  expect_equal(result1$dr2ij, result2$dr2ij)
})
