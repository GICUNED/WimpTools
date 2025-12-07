test_that("plot optimization functions work correctly", {
  
  # Test smart label positioning with simple data
  x_coords <- c(0.2, 0.5, 0.8)
  y_coords <- c(0.3, 0.7, 0.4)
  labels <- c("Label A", "Label B", "Label C")
  
  # Test that the function runs without error
  expect_no_error({
    result <- .smart_label_positions(x_coords, y_coords, labels)
  })
  
  # Test the result structure
  result <- .smart_label_positions(x_coords, y_coords, labels)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
  expect_true(all(c("x", "y", "label", "xshift", "yshift") %in% colnames(result)))
  
  # Test that shifts are numeric
  expect_true(is.numeric(result$xshift))
  expect_true(is.numeric(result$yshift))
  
  # Test with edge cases
  expect_no_error(.smart_label_positions(c(), c(), c()))
  expect_no_error(.smart_label_positions(c(0.5), c(0.5), c("Single")))
})

test_that("plot optimization handles invalid inputs", {
  
  # Test with mismatched vector lengths
  expect_error(.smart_label_positions(c(0.1, 0.2), c(0.3), c("A", "B")))
  
  # Test with NA coordinates
  x_with_na <- c(0.2, NA, 0.8)
  y_with_na <- c(0.3, 0.7, 0.4)
  labels <- c("A", "B", "C")
  
  result <- .smart_label_positions(x_with_na, y_with_na, labels)
  expect_equal(nrow(result), 2)  # Should filter out NA
})

test_that("label overlap detection works correctly", {
  
  # Test no overlap
  no_overlap <- .labels_overlap(0.1, 0.1, 0.05, 0.03, 0.5, 0.5, 0.05, 0.03)
  expect_false(no_overlap)
  
  # Test clear overlap
  overlap <- .labels_overlap(0.1, 0.1, 0.05, 0.03, 0.12, 0.11, 0.05, 0.03)
  expect_true(overlap)
  
  # Test point coverage
  covers <- .label_covers_point(0.5, 0.5, 0.1, 0.05, 0.5, 0.5, 0.01)
  expect_true(covers)
  
  no_cover <- .label_covers_point(0.5, 0.5, 0.1, 0.05, 0.8, 0.8, 0.01)
  expect_false(no_cover)
})