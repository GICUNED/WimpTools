# Test suite for change implications functions

test_that("if_index works correctly", {
  # Load test data  
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test if_index
  expect_no_error(if_result <- if_index(su_wimp))
  expect_true(is.data.frame(if_result) || is.matrix(if_result))
  expect_true(nrow(if_result) > 0)
  
  # Test with different standardization methods
  expect_no_error(if_index(su_wimp, std = "adjacent"))
  expect_no_error(if_index(su_wimp, std = "normalised"))
})

test_that("if_plot works correctly", {
  # Load test data
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test if_plot
  expect_no_error(if_plot_result <- if_plot(su_wimp))
  expect_true(inherits(if_plot_result, c("plotly", "htmlwidget")))
  
  # Test with different parameters
  expect_no_error(if_plot(su_wimp, show = "all"))
  expect_no_error(if_plot(su_wimp, center = "data"))
  expect_no_error(if_plot(su_wimp, text.size = 0.8))
})

test_that("if_barchart works correctly", {
  # Load test data
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  
  # Test if_barchart
  expect_no_error(if_bar_result <- if_barchart(su_wimp))
  expect_true(inherits(if_bar_result, c("plotly", "htmlwidget")))
  
  # Test with different parameters
  expect_no_error(if_barchart(su_wimp, show = "all"))
})

test_that("change implications functions handle errors", {
  # Test error handling for invalid inputs
  expect_error(if_index(NULL))
  expect_error(if_plot(list()))
  expect_error(if_barchart("not_a_wimp"))
})