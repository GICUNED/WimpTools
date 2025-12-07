test_that("monitoring functions work correctly", {
  
  # Load test data
  data(su_wimp, package = "WimpTools", envir = environment())
  
  # Test monitoring_ssi
  ssi_result <- monitoring_ssi(su_wimp)
  expect_true("plotly" %in% class(ssi_result) || "htmlwidget" %in% class(ssi_result))
  
  # Test monitoring_adj
  adj_result <- monitoring_adj(su_wimp)
  expect_true("plotly" %in% class(adj_result) || "htmlwidget" %in% class(adj_result))
  
  # Test monitoring_ph
  ph_result <- monitoring_ph(su_wimp)
  expect_true("plotly" %in% class(ph_result) || "htmlwidget" %in% class(ph_result))
})

test_that("monitoring functions handle edge cases", {
  
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
  expect_no_error(monitoring_ssi(minimal_wimp))
  expect_no_error(monitoring_adj(minimal_wimp))
  expect_no_error(monitoring_ph(minimal_wimp))
  
  # Test invalid inputs
  expect_error(monitoring_ssi(NULL))
  expect_error(monitoring_adj(list()))
  expect_error(monitoring_ph("not_a_wimp"))
})