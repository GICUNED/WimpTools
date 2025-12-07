test_that("graph functions work correctly", {
  
  # Load test data
  data(su_wimp, package = "WimpTools", envir = environment())
  
  # Test digraph
  digraph_result <- digraph(su_wimp)
  expect_true("plotly" %in% class(digraph_result) || "htmlwidget" %in% class(digraph_result))
  
  # Test digraph with different parameters
  digraph_custom <- digraph(su_wimp, title = "Test Digraph", threshold = 0.6)
  expect_true("plotly" %in% class(digraph_custom) || "htmlwidget" %in% class(digraph_custom))
  
  # Test idealdigraph  
  idealdigraph_result <- idealdigraph(su_wimp)
  expect_true("plotly" %in% class(idealdigraph_result) || "htmlwidget" %in% class(idealdigraph_result))
  
  # Test simdigraph
  simdigraph_result <- simdigraph(su_wimp)
  expect_true("plotly" %in% class(simdigraph_result) || "htmlwidget" %in% class(simdigraph_result))
  
  # Test inout_digraph
  inout_result <- inout_digraph(su_wimp)
  expect_true("plotly" %in% class(inout_result) || "htmlwidget" %in% class(inout_result))
})

test_that("graph functions handle edge cases", {
  
  # Test with minimal wimp object
  minimal_wimp <- list(
    global = list(
      wmatrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
      ideal_similarity = 0.62
    ),
    vertices = data.frame(
      construct = c("A", "B"),
      pos_x = c(0.3, 0.7),
      pos_y = c(0.5, 0.5)
    )
  )
  
  # Test that functions don't crash with minimal data
  expect_no_error(digraph(minimal_wimp))
  expect_no_error(idealdigraph(minimal_wimp))
  expect_no_error(simdigraph(minimal_wimp))
  expect_no_error(inout_digraph(minimal_wimp))
  
  # Test invalid inputs
  expect_error(digraph(NULL))
  expect_error(idealdigraph(list()))
  expect_error(simdigraph("not_a_wimp"))
  expect_error(inout_digraph(matrix(1:4, nrow = 2)))
})

test_that("graph threshold parameters work correctly", {
  
  # Load test data
  data(su_wimp, package = "WimpTools", envir = environment())
  
  # Test different threshold values
  expect_no_error(digraph(su_wimp, threshold = 0.1))
  expect_no_error(digraph(su_wimp, threshold = 0.5))
  expect_no_error(digraph(su_wimp, threshold = 0.9))
  
  # Test edge case thresholds
  expect_no_error(digraph(su_wimp, threshold = 0))
  expect_no_error(digraph(su_wimp, threshold = 1))
  
  # Test invalid thresholds should not crash but might give warnings
  expect_no_error(digraph(su_wimp, threshold = -0.1))
  expect_no_error(digraph(su_wimp, threshold = 1.1))
})