# Test suite for adjustment functions

test_that("construct_index works correctly", {
  # Load test data  
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test construct_index
  expect_no_error(construct_result <- construct_index(su_wimp))
  expect_true(is.data.frame(construct_result))
  expect_true("construct" %in% colnames(construct_result))
  expect_true(nrow(construct_result) > 0)
})

test_that("self_index works correctly", {
  # Load test data
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test self_index with different methods
  expect_no_error(ssi_result <- self_index(su_wimp, method = "ssi"))
  expect_true(is.numeric(ssi_result) || is.data.frame(ssi_result))
  
  expect_no_error(euc_result <- self_index(su_wimp, method = "euclidean"))
  expect_true(is.numeric(euc_result) || is.data.frame(euc_result))
  
  expect_no_error(int_result <- self_index(su_wimp, method = "intensity"))
  expect_true(is.numeric(int_result) || is.data.frame(int_result))
})

test_that("self_plot works correctly", {
  # Load test data
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test self_plot
  expect_no_error(plot_result <- self_plot(su_wimp))
  expect_true(inherits(plot_result, c("plotly", "htmlwidget")))
})

test_that("ssi_heatmap works correctly", {
  # Load test data
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test ssi_heatmap
  expect_no_error(heatmap_result <- ssi_heatmap(su_wimp))
  expect_true(inherits(heatmap_result, c("plotly", "htmlwidget")))
})

test_that("hypo_plot works correctly", {
  # Load test data
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test hypo_plot
  expect_no_error(hypo_result <- hypo_plot(su_wimp))
  expect_true(inherits(hypo_result, c("plotly", "htmlwidget")))
  
  # Test with different parameters
  expect_no_error(hypo_plot(su_wimp, text.size = 0.5))
  expect_no_error(hypo_plot(su_wimp, show.labels = FALSE))
})

test_that("adjustment functions handle errors", {
  # Test error handling for invalid inputs
  expect_error(construct_index(NULL))
  expect_error(self_index(list()))
  expect_error(self_plot("not_a_wimp"))
  expect_error(ssi_heatmap(123))
  expect_error(hypo_plot(data.frame()))
})