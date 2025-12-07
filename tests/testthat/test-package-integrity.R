test_that("package data loads correctly", {
  
  # Test that su_wimp data loads
  data(su_wimp, package = "WimpTools", envir = environment())
  expect_true(exists("su_wimp"))
  expect_true(is.list(su_wimp))
  
  # Test that example_wimp data loads  
  data(example_wimp, package = "WimpTools", envir = environment())
  expect_true(exists("example_wimp"))
  expect_true(is.list(example_wimp))
  
  # Test data structure consistency
  expect_true("global" %in% names(su_wimp))
  expect_true("vertices" %in% names(su_wimp))
  expect_true("global" %in% names(example_wimp))
  expect_true("vertices" %in% names(example_wimp))
})

test_that("S3 methods work correctly", {
  
  # Load test data
  data(su_wimp, package = "WimpTools", envir = environment())
  
  # Test print method
  expect_output(print(su_wimp))
  
  # Test summary method  
  expect_output(summary(su_wimp))
  expect_true(is.list(summary(su_wimp)))
})

test_that("package namespace is correct", {
  
  # Test that main functions are exported
  expect_true(exists("importwimp"))
  expect_true(exists("degree_index"))
  expect_true(exists("adj_plot"))
  expect_true(exists("hypo_plot"))
  expect_true(exists("digraph"))
  expect_true(exists("wimp_biplot"))
  expect_true(exists("construct_index"))
  
  # Test that internal functions are not exported
  expect_false(exists(".smart_label_positions", envir = globalenv()))
})

test_that("package dependencies work", {
  
  # Test required packages
  expect_true(requireNamespace("plotly", quietly = TRUE))
  expect_true(requireNamespace("ggplot2", quietly = TRUE))
  expect_true(requireNamespace("readxl", quietly = TRUE))
  expect_true(requireNamespace("dplyr", quietly = TRUE))
})