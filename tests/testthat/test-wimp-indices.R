test_that("WIMP indices work correctly", {
  
  # Load test data
  data(su_wimp, package = "WimpTools", envir = environment())
  
  # Test construct_index
  construct_result <- construct_index(su_wimp)
  expect_true(is.data.frame(construct_result))
  expect_true("construct" %in% colnames(construct_result))
  expect_true("construct_index" %in% colnames(construct_result))
  expect_true(all(construct_result$construct_index >= 0))
  expect_true(all(construct_result$construct_index <= 1))
})

test_that("WIMP indices handle edge cases", {
  
  # Test with minimal wimp object
  minimal_wimp <- list(
    global = list(
      wmatrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    ),
    vertices = data.frame(
      construct = c("A", "B"),
      pos_x = c(0.3, 0.7),
      pos_y = c(0.5, 0.5)
    )
  )
  
  # Test that function doesn't crash with minimal data
  expect_no_error(construct_index(minimal_wimp))
  
  # Test invalid inputs
  expect_error(construct_index(NULL))
  expect_error(construct_index(list()))
  expect_error(construct_index("not_a_wimp"))
})