test_that("ScDeconvExperiment creates a valid object", {
  object <- make_test_object()

  expect_s4_class(object, "ScDeconvExperiment")
  expect_true(methods::validObject(object))
  expect_equal(nrow(object), 5)
  expect_equal(ncol(object), 4)
  expect_true("counts" %in% SummarizedExperiment::assayNames(object))
  expect_type(object@params, "list")
  expect_type(object@history, "character")
})

test_that("ScDeconvExperiment rejects invalid input", {
  expect_error(
    ScDeconvExperiment(matrix(1:4, nrow = 2)),
    "'sce' must be a SingleCellExperiment object.",
    fixed = TRUE
  )

  sce_without_counts <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = matrix(1:4, nrow = 2))
  )

  expect_error(
    ScDeconvExperiment(sce_without_counts),
    "counts"
  )
})
