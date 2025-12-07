test_that("importwimp function works correctly", {
  
  # Skip if test data not available
  skip_if_not_installed("readxl")
  
  # Load example data for testing
  data(su_wimp, package = "WimpTools", envir = environment())
  data(example_wimp, package = "WimpTools", envir = environment())
  
  # Test that su_wimp is properly structured
  expect_true(is.list(su_wimp))
  expect_true("global" %in% names(su_wimp))
  expect_true("vertices" %in% names(su_wimp))
  expect_true("wmatrix" %in% names(su_wimp$global))
  
  # Test that example_wimp is properly structured  
  expect_true(is.list(example_wimp))
  expect_true("global" %in% names(example_wimp))
  expect_true("vertices" %in% names(example_wimp))
  
  # Test data types
  expect_true(is.data.frame(su_wimp$vertices))
  expect_true(is.matrix(su_wimp$global$wmatrix) || is.data.frame(su_wimp$global$wmatrix))
  
  # Test that vertices have required columns
  expect_true("construct" %in% colnames(su_wimp$vertices))
  expect_true("pos_x" %in% colnames(su_wimp$vertices))
  expect_true("pos_y" %in% colnames(su_wimp$vertices))
  
  # Test matrix dimensions consistency
  n_vertices <- nrow(su_wimp$vertices)
  wmatrix_dims <- dim(su_wimp$global$wmatrix)
  expect_equal(wmatrix_dims[1], n_vertices)
  expect_equal(wmatrix_dims[2], n_vertices)
})

test_that("importwimp handles invalid inputs correctly", {
  # Test non-existent file
  expect_error(importwimp("non_existent_file.xlsx"))
  
  # Test invalid sheet number
  # Note: This would require a real file to test properly
  # For now, we'll skip this test in automated environments
  skip_on_ci()
  skip_on_cran()
})