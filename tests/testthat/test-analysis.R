test_that("identifyHVGs marks the requested number of HVGs", {
  object <- make_test_object()
  object <- normalizeCPM(object)
  object <- logTransform(object)

  object <- identifyHVGs(
    object,
    assay_name = "logcounts",
    n_top = 2,
    method = "variance"
  )

  rd <- SummarizedExperiment::rowData(object)

  expect_true("is_hvg" %in% colnames(rd))
  expect_true("hvg_rank" %in% colnames(rd))
  expect_equal(sum(rd$is_hvg), 2)
  expect_equal(object@params$identifyHVGs$n_selected, 2)
})

test_that("identifyHVGs validates arguments", {
  object <- make_test_object()

  expect_error(identifyHVGs(object, assay_name = "missing"), "assay_name")

  object <- normalizeCPM(object)
  object <- logTransform(object)

  expect_error(identifyHVGs(object, n_top = 0), "n_top")
  expect_error(identifyHVGs(object, min_mean = -1), "min_mean")
  expect_error(identifyHVGs(object, min_mean = 2, max_mean = 1), "max_mean")
  expect_error(identifyHVGs(object, method = "bad_method"), "arg")
})

test_that(".select_clustering_rows selects HVGs by rank", {
  object <- make_test_object()
  SummarizedExperiment::rowData(object)$is_hvg <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  SummarizedExperiment::rowData(object)$hvg_rank <- c(3, NA, 1, 2, NA)

  selected <- .select_clustering_rows(object, use_hvg = TRUE, n_hvg = 2)

  expect_equal(selected, c(3, 4))
  expect_equal(.select_clustering_rows(object, use_hvg = FALSE), seq_len(nrow(object)))
})

test_that(".select_clustering_rows errors when HVGs are unavailable", {
  object <- make_test_object()

  expect_error(.select_clustering_rows(object), "HVG column not found")

  SummarizedExperiment::rowData(object)$is_hvg <- rep(FALSE, nrow(object))

  expect_error(.select_clustering_rows(object), "No HVGs found")
})
