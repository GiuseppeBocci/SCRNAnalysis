test_that("normalizeCPM stores CPM-normalized counts", {
  object <- make_test_object()
  counts <- SummarizedExperiment::assay(object, "counts")
  expected_library_size <- colSums(counts)
  expected <- t(t(counts) / expected_library_size) * 1e6

  object <- normalizeCPM(object)

  expect_true("normcounts" %in% SummarizedExperiment::assayNames(object))
  expect_equal(SummarizedExperiment::assay(object, "normcounts"), expected)
  expect_equal(SummarizedExperiment::colData(object)$library_size, expected_library_size)
  expect_equal(object@params$normalization$method, "CPM")
  expect_true("normalizeCPM" %in% object@history)
})

test_that("normalizeToTotal uses the requested target total", {
  object <- make_test_object()
  counts <- SummarizedExperiment::assay(object, "counts")
  target_total <- 1000
  expected <- t(t(counts) / colSums(counts)) * target_total

  object <- normalizeToTotal(object, target_total = target_total)

  expect_equal(SummarizedExperiment::assay(object, "normcounts"), expected)
  expect_equal(object@params$normalization$method, "target_total")
  expect_equal(object@params$normalization$target_total, target_total)
})

test_that("normalizeToTotal validates target_total", {
  object <- make_test_object()

  expect_error(normalizeToTotal(object, target_total = 0), "target_total")
  expect_error(normalizeToTotal(object, target_total = -1), "target_total")
  expect_error(normalizeToTotal(object, target_total = c(1, 2)), "target_total")
})

test_that("logTransform creates logcounts assay", {
  object <- make_test_object()
  object <- normalizeToTotal(object, target_total = 1000)

  normcounts <- SummarizedExperiment::assay(object, "normcounts")
  expected <- log(normcounts + 1, base = 2)

  object <- logTransform(object)

  expect_true("logcounts" %in% SummarizedExperiment::assayNames(object))
  expect_equal(SummarizedExperiment::assay(object, "logcounts"), expected)
  expect_equal(object@params$log_transform$input_assay, "normcounts")
  expect_true("logTransform" %in% object@history)
})

test_that("logTransform validates input", {
  object <- make_test_object()

  expect_error(logTransform(object), "Assay not found")

  object <- normalizeCPM(object)

  expect_error(logTransform(object, pseudo_count = -1), "pseudo_count")
  expect_error(logTransform(object, base = 1), "base")
  expect_error(logTransform(object, base = 0), "base")
})

test_that("normalizeMedianLibrarySize scales to median library size", {
  object <- make_test_object()
  counts <- SummarizedExperiment::assay(object, "counts")
  library_size <- colSums(counts)
  median_library_size <- stats::median(library_size)
  expected <- t(t(counts) / library_size) * median_library_size

  object <- normalizeMedianLibrarySize(object)

  expect_equal(SummarizedExperiment::assay(object, "normcounts"), expected)
  expect_equal(object@params$normalization$method, "median_library_size")
  expect_equal(object@params$normalization$target_total, median_library_size)
})

test_that("normalization rejects zero library size cells", {
  counts <- matrix(
    c(
      1,0,
      2,0,
      3,0
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  colnames(counts) <- c("Cell1", "Cell2")

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  object <- ScDeconvExperiment(sce)

  expect_error(normalizeCPM(object), "library size equal to zero")
})

test_that("normalizeToTotal preserves sparse count assays", {
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 3, 2, 4, 0), nrow = 3),
    sparse = TRUE
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))

  object <- ScDeconvExperiment(
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  )
  normalized <- normalizeToTotal(object, target_total = 1e4)

  expect_s4_class(SummarizedExperiment::assay(normalized, "normcounts"), "dgCMatrix")
  expect_identical(
    dimnames(SummarizedExperiment::assay(normalized, "normcounts")),
    dimnames(counts)
  )
  expect_equal(
    Matrix::colSums(SummarizedExperiment::assay(normalized, "normcounts")),
    stats::setNames(rep(1e4, ncol(normalized)), colnames(normalized))
  )
})
