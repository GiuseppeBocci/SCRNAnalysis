test_that("compareClusterCellTypes stores count and proportion tables", {
  object <- make_comparison_object()

  object <- compareClusterCellTypes(object)
  comparison <- S4Vectors::metadata(object)$cluster_cell_type_comparison

  expect_type(comparison, "list")
  expect_true(all(c("counts", "row_proportions", "column_proportions", "long_table") %in% names(comparison)))
  expect_equal(sum(comparison$counts), ncol(object))
  expect_true("compareClusterCellTypes" %in% names(object@params))
})

test_that("compareClusterCellTypes validates required columns", {
  object <- make_test_object()

  expect_error(compareClusterCellTypes(object), "cluster_column")

  SummarizedExperiment::colData(object)$cluster <- factor(c("A", "A", "B", "B"))

  expect_error(compareClusterCellTypes(object), "cell_type_column")
})


