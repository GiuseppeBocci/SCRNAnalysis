test_that("filterLowExpressionGenes removes genes below min_cells", {
  object <- make_test_object()

  filtered <- filterLowExpressionGenes(
    object,
    min_count = 1,
    min_cells = 3,
    store_metrics = TRUE
  )

  expect_lt(nrow(filtered), nrow(object))
  expect_true(all(rowSums(SummarizedExperiment::assay(filtered, "counts") >= 1) >= 3))
  expect_true("gene_filter_keep" %in% colnames(SummarizedExperiment::rowData(filtered)))
  expect_equal(filtered@params$filterLowExpressionGenes$min_cells, 3)
})

test_that("filterLowExpressionGenes validates assay and percentile", {
  object <- make_test_object()

  expect_error(filterLowExpressionGenes(object, assay_name = "missing"), "assay_name")
  expect_error(filterLowExpressionGenes(object, lower_percentile = -0.1), "lower_percentile")
  expect_error(filterLowExpressionGenes(object, lower_percentile = 1.1), "lower_percentile")
})

test_that("filterCellsByDetectedGenes stores metrics and filters cells", {
  object <- make_test_object()

  filtered <- filterCellsByDetectedGenes(
    object,
    min_count = 1,
    lower_percentile = 0,
    upper_percentile = 1,
    store_metrics = TRUE
  )

  expect_equal(ncol(filtered), ncol(object))
  expect_true("detected_genes" %in% colnames(SummarizedExperiment::colData(filtered)))
  expect_true("total_counts" %in% colnames(SummarizedExperiment::colData(filtered)))
  expect_true("detected_genes_filter_keep" %in% colnames(SummarizedExperiment::colData(filtered)))
  expect_equal(filtered@params$filterCellsByDetectedGenes$cells_before, ncol(object))
})

test_that("filterCellsByDetectedGenes validates percentiles", {
  object <- make_test_object()

  expect_error(filterCellsByDetectedGenes(object, assay_name = "missing"), "assay_name")
  expect_error(filterCellsByDetectedGenes(object, lower_percentile = -0.1), "Percentiles")
  expect_error(filterCellsByDetectedGenes(object, upper_percentile = 1.1), "Percentiles")
  expect_error(
    filterCellsByDetectedGenes(object, lower_percentile = 0.8, upper_percentile = 0.2),
    "lower_percentile"
  )
})
