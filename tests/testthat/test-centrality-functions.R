# Test suite for centrality functions

test_that("centrality functions work correctly", {
  # Load test data  
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test degree_index
  expect_no_error(degree_result <- degree_index(su_wimp))
  expect_true(is.matrix(degree_result))
  expect_true(ncol(degree_result) == 3)
  expect_true(all(c("Out", "In", "All") %in% colnames(degree_result)))
  expect_true(all(degree_result >= 0))
  
  # Test close_index
  expect_no_error(close_result <- close_index(su_wimp))
  expect_true(is.matrix(close_result) || is.numeric(close_result))
  
  # Test betw_index
  expect_no_error(betw_result <- betw_index(su_wimp))
  expect_true(is.matrix(betw_result) || is.numeric(betw_result))
  
  # Test eigen_index
  expect_no_error(eigen_result <- eigen_index(su_wimp))
  expect_true(is.matrix(eigen_result) || is.numeric(eigen_result))
  
  # Test density_index
  expect_no_error(density_result <- density_index(su_wimp))
  expect_true(is.numeric(density_result))
  expect_true(length(density_result) == 1)
  expect_true(density_result >= 0 && density_result <= 1)
  
  # Test auc_index
  expect_no_error(auc_result <- auc_index(su_wimp))
  expect_true(is.numeric(auc_result))
  expect_true(length(auc_result) == 1)
})

test_that("centrality functions handle different graph types", {
  # Load test data
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test with different modes for degree_index (if applicable)
  expect_no_error(degree_index(su_wimp))
  
  # Test close_index with different modes
  expect_no_error(close_index(su_wimp, mode = "in"))
  expect_no_error(close_index(su_wimp, mode = "out"))
  expect_no_error(close_index(su_wimp, mode = "all"))
})

test_that("centrality functions handle errors", {
  # Test error handling for invalid inputs
  expect_error(degree_index(NULL))
  expect_error(close_index(list()))
  expect_error(betw_index("not_a_wimp"))
  expect_error(eigen_index(123))
  expect_error(density_index(data.frame()))
  expect_error(auc_index(c(1, 2, 3)))
})