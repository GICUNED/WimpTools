test_that("importrepgrid_json reads a RepGrid record and keeps the metadata", {
  skip_if_not_installed("OpenRepGrid")
  json <- '{"id": "abc", "type": "repgrid", "title": "Demo", "status": "complete",
    "patientId": null, "notes": null, "params": null,
    "data": {"scaleMin": 1, "scaleMax": 5,
      "elements": ["Self", "Mother", "Ideal"],
      "constructs": [
        {"left": "calm", "right": "anxious", "ratings": [2, 4, 1]},
        {"left": "open", "right": "closed", "ratings": [3, 5, 2]},
        {"left": "active", "right": "passive", "ratings": [1, 2, 1]}]}}'
  rg <- importrepgrid_json(json)
  expect_s4_class(rg, "repgrid")
  expect_equal(unname(OpenRepGrid::getScale(rg)), c(1, 5))
  expect_equal(dim(OpenRepGrid::ratings(rg)), c(3L, 3L))
  expect_equal(unname(OpenRepGrid::ratings(rg)[2, ]), c(3, 5, 2))
  expect_equal(OpenRepGrid::elements(rg), c("Self", "Mother", "Ideal"))
  expect_equal(OpenRepGrid::constructs(rg)$leftpole, c("calm", "open", "active"))
  expect_equal(rg@meta$title, "Demo")
  expect_null(rg@meta$patientId)
})

test_that("importrepgrid_json accepts files and parsed lists", {
  skip_if_not_installed("OpenRepGrid")
  rec <- list(type = "repgrid", data = list(
    scaleMin = 1, scaleMax = 3, elements = list("a", "b"),
    constructs = list(list(left = "l1", right = "r1", ratings = list(1, 3)),
                      list(left = "l2", right = "r2", ratings = list(2, 2)))))
  tmp <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE), tmp)
  expect_equal(OpenRepGrid::ratings(importrepgrid_json(tmp)),
               OpenRepGrid::ratings(importrepgrid_json(rec)))
})

test_that("importrepgrid_json validates the record", {
  skip_if_not_installed("OpenRepGrid")
  expect_error(importrepgrid_json('{"type": "wimpgrid", "data": {}}'), "importwimp")
  expect_error(importrepgrid_json('{"type": "repgrid"}'), "no `data`")
  bad_len <- list(type = "repgrid", data = list(scaleMin = 1, scaleMax = 5,
    elements = list("a", "b"), constructs = list(
      list(left = "l", right = "r", ratings = list(1)),
      list(left = "l", right = "r", ratings = list(1, 2)))))
  expect_error(importrepgrid_json(bad_len), "ratings but there are")
  out <- list(type = "repgrid", data = list(scaleMin = 1, scaleMax = 5,
    elements = list("a", "b"), constructs = list(
      list(left = "l", right = "r", ratings = list(1, 9)),
      list(left = "l", right = "r", ratings = list(1, 2)))))
  expect_error(importrepgrid_json(out), "outside the scale")
  expect_error(importrepgrid_json("{not json"), "not a valid JSON")
})

test_that("an imported grid works with the widgets", {
  skip_if_not_installed("OpenRepGrid")
  x <- OpenRepGrid::feixas2004
  rec <- list(type = "repgrid", title = "t", data = list(
    scaleMin = 1, scaleMax = 7,
    elements = as.list(OpenRepGrid::elements(x)),
    constructs = lapply(seq_len(nrow(OpenRepGrid::ratings(x))), function(i) list(
      left = OpenRepGrid::constructs(x)$leftpole[i],
      right = OpenRepGrid::constructs(x)$rightpole[i],
      ratings = as.list(unname(OpenRepGrid::ratings(x)[i, ]))))))
  rg <- importrepgrid_json(rec)
  expect_equal(unname(OpenRepGrid::ratings(rg)), unname(OpenRepGrid::ratings(x)))
  expect_s3_class(repgrid_biplot(rg), "plotly")
  expect_s3_class(repgrid_dilemmas(rg), "plotly")
  expect_equal(nrow(repgrid_indices(rg)$global), 15)
})
