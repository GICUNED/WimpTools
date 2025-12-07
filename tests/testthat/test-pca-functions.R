test_that("PCA functions work correctly", {
  
  # Load test data
  data(su_wimp, package = "WimpTools", envir = environment())
  
  # Test wimp_biplot
  biplot_result <- wimp_biplot(su_wimp)
  expect_true("plotly" %in% class(biplot_result) || "htmlwidget" %in% class(biplot_result))
  
  # Test wimp_biplot with different parameters
  biplot_custom <- wimp_biplot(su_wimp, title = "Test Biplot")
  expect_true("plotly" %in% class(biplot_custom) || "htmlwidget" %in% class(biplot_custom))
})

test_that("PCA functions handle edge cases", {
  
  # Test with minimal wimp object
  minimal_wimp <- list(
    global = list(
      wmatrix = matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), nrow = 3)
    ),
    vertices = data.frame(
      construct = c("A", "B", "C"),
      pos_x = c(0.2, 0.5, 0.8),
      pos_y = c(0.3, 0.7, 0.5)
    )
  )
  
  # Test that function doesn't crash with minimal data
  expect_no_error(wimp_biplot(minimal_wimp))
  
  # Test invalid inputs
  expect_error(wimp_biplot(NULL))
  expect_error(wimp_biplot(list()))
  expect_error(wimp_biplot("not_a_wimp"))
})