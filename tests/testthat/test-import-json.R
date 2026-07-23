library(testthat)
library(jsonlite)
library(WimpTools)

test_that("importwimp_json works correctly with valid JSON string", {
  # 1. Standard grid JSON with [min, max] scale array
  json_str <- '{
    "metadata": {
      "ID": "TST",
      "age": "25",
      "gender": "male"
    },
    "scale": [1, 7],
    "constructs": [
      {"left_pole": "Anxious", "right_pole": "Calm", "self": 3, "ideal": 6},
      {"left_pole": "Lazy", "right_pole": "Hard Worker", "self": 2, "ideal": 5}
    ],
    "hypo_matrix": [
      [3, 4],
      [2, 5]
    ]
  }'

  expect_no_error({
    wimp <- importwimp_json(json_str)
  })

  expect_true(inherits(wimp, "wimp"))
  expect_equal(wimp$global$ID, "TST")
  expect_equal(wimp$global$age, "25")
  expect_equal(wimp$global$gender, "male")
  expect_equal(wimp$global$scale, c(1, 7))
  expect_equal(wimp$global$n_constructs, 2)
  
  # Vertices df structure
  expect_s3_class(wimp$vertices, "data.frame")
  expect_equal(nrow(wimp$vertices), 2)
  expect_true(all(c("left_pole", "right_pole", "self", "ideal", "self_pole", "ideal_pole", "congruence") %in% colnames(wimp$vertices)))
  
  # Check vertices values (scale min=1, max=7, center=4, half_range=3)
  # self=3 -> (3-4)/3 = -0.333
  expect_equal(wimp$vertices$self[1], -1/3, tolerance = 1e-4)
  # ideal=6 -> (6-4)/3 = 0.666
  expect_equal(wimp$vertices$ideal[1], 2/3, tolerance = 1e-4)
  
  # Edges structure
  expect_s3_class(wimp$edges, "data.frame")
  expect_true(all(c("id", "from", "to", "weight") %in% colnames(wimp$edges)))
})

test_that("importwimp_json handles scale object, preference column and other extra attributes", {
  # 2. Scale as object and extra columns
  json_str <- '{
    "scale": {"min": 1, "max": 7},
    "constructs": [
      {
        "left_pole": "Anxious",
        "right_pole": "Calm",
        "self": 2,
        "ideal": 6,
        "preference": 5,
        "category": "Emotional"
      },
      {
        "left_pole": "Lazy",
        "right_pole": "Hard Worker",
        "self": 2,
        "ideal": 6,
        "preference": 4,
        "category": "Work"
      }
    ],
    "hypo_matrix": [
      [2, 3],
      [4, 2]
    ]
  }'

  wimp <- importwimp_json(json_str)

  expect_equal(wimp$global$scale, c(1, 7))
  expect_true("preference" %in% colnames(wimp$vertices))
  expect_true("category" %in% colnames(wimp$vertices))
  
  # preference=5 -> normalized to (5-4)/3 = 0.333
  expect_equal(wimp$vertices$preference[1], 1/3, tolerance = 1e-4)
  # category should be a character column
  expect_equal(wimp$vertices$category[1], "Emotional")
  expect_equal(wimp$vertices$category[2], "Work")
})

test_that("importwimp_json accepts pre-parsed lists", {
  grid_list <- list(
    scale = c(1, 7),
    constructs = data.frame(
      left_pole = c("Anxious", "Lazy"),
      right_pole = c("Calm", "Hard Worker"),
      self = c(3, 2),
      ideal = c(6, 5),
      stringsAsFactors = FALSE
    ),
    hypo_matrix = matrix(c(3, 2, 4, 5), nrow = 2, ncol = 2)
  )

  expect_no_error({
    wimp <- importwimp_json(grid_list)
  })

  expect_true(inherits(wimp, "wimp"))
  expect_equal(wimp$global$scale, c(1, 7))
  expect_equal(wimp$global$n_constructs, 2)
})

test_that("importwimp_json triggers warnings for out-of-scale ratings", {
  # Ratings 0 and 8 are out of bounds for [1, 7]
  json_str <- '{
    "scale": [1, 7],
    "constructs": [
      {"left_pole": "Anxious", "right_pole": "Calm", "self": 0, "ideal": 6},
      {"left_pole": "Lazy", "right_pole": "Hard Worker", "self": 3, "ideal": 8}
    ],
    "hypo_matrix": [
      [3, 4],
      [2, 5]
    ]
  }'

  expect_warning(importwimp_json(json_str), "ratings are out of scale bounds")
})

test_that("importwimp_json handles error conditions gracefully", {
  # 1. Invalid JSON syntax
  expect_error(importwimp_json("{invalid_json}"), "lexical error")

  # 2. Missing required keys (no hypo_matrix)
  json_no_hypo <- '{
    "scale": [1, 7],
    "constructs": [
      {"left_pole": "Anxious", "right_pole": "Calm", "self": 3, "ideal": 6}
    ]
  }'
  expect_error(importwimp_json(json_no_hypo), "Required field missing")

  # 3. Dimension mismatch (3 constructs, but 2x2 matrix)
  json_mismatch <- '{
    "scale": [1, 7],
    "constructs": [
      {"left_pole": "Anxious", "right_pole": "Calm", "self": 3, "ideal": 6},
      {"left_pole": "Lazy", "right_pole": "Hard Worker", "self": 2, "ideal": 5},
      {"left_pole": "Sad", "right_pole": "Happy", "self": 4, "ideal": 6}
    ],
    "hypo_matrix": [
      [3, 4],
      [2, 5]
    ]
  }'
  expect_error(importwimp_json(json_mismatch), "hypo_matrix must be a square matrix of size 3 x 3")
})

test_that("JSON-imported wimp is compatible with downstream package functions", {
  # Test compatibility with downstream package analysis functions
  json_str <- '{
    "scale": [1, 7],
    "constructs": [
      {"left_pole": "Anxious", "right_pole": "Calm", "self": 3, "ideal": 6, "preference": 6},
      {"left_pole": "Lazy", "right_pole": "Hard Worker", "self": 2, "ideal": 5, "preference": 5}
    ],
    "hypo_matrix": [
      [3, 4],
      [2, 5]
    ]
  }'

  wimp <- importwimp_json(json_str)

  # Check that digraph works
  expect_no_error({
    dg <- digraph(wimp)
  })
  expect_true(inherits(dg, "visNetwork"))

  # Check that self_index works
  expect_no_error({
    si <- self_index(wimp, method = "ssi")
  })
  expect_true(inherits(si, "self_index"))

  # Check that construct_index works
  expect_no_error({
    ci <- construct_index(wimp)
  })
  expect_true(is.matrix(ci))
})
