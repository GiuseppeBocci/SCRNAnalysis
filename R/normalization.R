###########
# Helpers #
###########

# Internal helper: get raw counts
.get_counts <- function(x) {
  if (!"counts" %in% SummarizedExperiment::assayNames(x)) {
    stop("The object must contain a 'counts' assay.", call. = FALSE)
  }

  counts <- SummarizedExperiment::assay(x, "counts")

  if (is.null(dim(counts)) || length(dim(counts)) != 2) {
    expected_length <- nrow(x) * ncol(x)

    if (length(counts) != expected_length) {
      stop(
        "The 'counts' assay must be a two-dimensional matrix-like object.",
        call. = FALSE
      )
    }

    counts <- matrix(
      counts,
      nrow = nrow(x),
      ncol = ncol(x),
      dimnames = list(rownames(x), colnames(x))
    )
  }

  counts
}

# Internal helper: compute library size
.compute_library_size <- function(counts) {
  if (is.null(dim(counts)) || length(dim(counts)) != 2) {
    stop(
      "The 'counts' assay must be a two-dimensional matrix-like object. ",
      "This can happen after subsetting without drop = FALSE.",
      call. = FALSE
    )
  }

  # Use Matrix::colSums rather than base::colSums so sparse Matrix assays
  # such as dgCMatrix objects extracted from Seurat are handled correctly.
  library_size <- Matrix::colSums(counts)
  library_size <- as.numeric(library_size)
  names(library_size) <- colnames(counts)

  if (any(library_size == 0)) {
    stop("Some cells have library size equal to zero.", call. = FALSE)
  }

  library_size
}

# Internal helper: normalize each column of a matrix-like count assay.
#
# Right multiplication by a diagonal matrix preserves sparse Matrix classes,
# unlike the transpose/division idiom, which can dispatch to methods that
# simplify a sparse assay to a vector.
.scale_columns <- function(counts, scale_factors) {
  if (length(scale_factors) != ncol(counts)) {
    stop("'scale_factors' must contain one value for each cell.", call. = FALSE)
  }

  scaled <- if (methods::is(counts, "Matrix")) {
    counts %*% Matrix::Diagonal(x = as.numeric(scale_factors))
  } else {
    sweep(counts, MARGIN = 2, STATS = scale_factors, FUN = "*")
  }

  dimnames(scaled) <- dimnames(counts)
  scaled
}



#' Normalize counts to CPM
#'
#' Normalizes raw counts by dividing each cell by its library size
#' and multiplying by 1,000,000.
#'
#' @param x A ScDeconvExperiment object.
#'
#' @return A ScDeconvExperiment object with CPM values stored in
#'   the `"normcounts"` assay.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(
#'   rpois(1000, lambda = 5),
#'   nrow = 100,
#'   ncol = 10
#' )
#' rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#' colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = counts)
#' )
#'
#' se <- ScDeconvExperiment(sce)
#'
#' normalized_se <- normalizeCPM(se)
#'
#' SummarizedExperiment::assayNames(normalized_se)
#'
#' @export
setGeneric("normalizeCPM", function(x) {
  standardGeneric("normalizeCPM")
})

#' @rdname normalizeCPM
#' @export
setMethod(
  "normalizeCPM",
  signature(x = "ScDeconvExperiment"),
  function(x) {
    counts <- .get_counts(x)
    library_size <- .compute_library_size(counts)

    normcounts <- .scale_columns(counts, 1e6 / library_size)

    SummarizedExperiment::assay(x, "normcounts") <- normcounts
    SummarizedExperiment::colData(x)$library_size <- library_size

    x@params$normalization <- list(
      method = "CPM",
      target_total = 1e6
    )
    x@history <- c(x@history, "normalizeCPM")
    methods::validObject(x)
    x
  }
)



#' Normalize counts to a user-defined total
#'
#' Normalizes raw counts by dividing each cell by its library size
#' and multiplying by a user-defined target total.
#'
#' @param x A ScDeconvExperiment object.
#' @param target_total Numeric scaling factor. Defaults to `1e6`,
#'   corresponding to CPM normalization.
#'
#' @return A ScDeconvExperiment object with normalized values stored in
#'   the `"normcounts"` assay.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(
#'   rpois(1000, lambda = 5),
#'   nrow = 100,
#'   ncol = 10
#' )
#' rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#' colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = counts)
#' )
#'
#' se <- ScDeconvExperiment(sce)
#'
#' normalized_se <- normalizeToTotal(
#'   se,
#'   target_total = 1e4
#' )
#'
#' SummarizedExperiment::assayNames(normalized_se)
#'
#' @export
setGeneric("normalizeToTotal", function(x, target_total = 1e6) {
  standardGeneric("normalizeToTotal")
})

#' @rdname normalizeToTotal
#' @export
setMethod(
  "normalizeToTotal",
  signature(x = "ScDeconvExperiment"),
  function(x, target_total = 1e6) {
    if (!is.numeric(target_total) || length(target_total) != 1 || target_total <= 0) {
      stop("'target_total' must be a single positive numeric value.")
    }

    counts <- .get_counts(x)
    library_size <- .compute_library_size(counts)

    normcounts <- .scale_columns(counts, target_total / library_size)

    SummarizedExperiment::assay(x, "normcounts") <- normcounts
    SummarizedExperiment::colData(x)$library_size <- library_size

    x@params$normalization <- list(
      method = "target_total",
      target_total = target_total
    )

    x@history <- c(x@history, "normalizeToTotal")
    methods::validObject(x)
    x
  }
)


#' Log-transform normalized counts
#'
#' Applies a log transformation to an assay and stores the result
#' in the `"logcounts"` assay.
#'
#' @param x A ScDeconvExperiment object.
#' @param assay_name Name of the assay to transform. Defaults to `"normcounts"`.
#' @param pseudo_count Numeric value added before log transformation.
#'   Defaults to `1`.
#' @param base Logarithm base. Defaults to `2`.
#'
#' @return A ScDeconvExperiment object with log-transformed values stored in
#'   the `"logcounts"` assay.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(
#'   rpois(1000, lambda = 5),
#'   nrow = 100,
#'   ncol = 10
#' )
#' rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#' colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'
#' normcounts <- t(t(counts) / colSums(counts)) * 1e4
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = counts,
#'     normcounts = normcounts
#'   )
#' )
#'
#' se <- ScDeconvExperiment(sce)
#'
#' log_transformed_se <- logTransform(
#'   se,
#'   assay_name = "normcounts",
#'   pseudo_count = 1,
#'   base = 2
#' )
#'
#' SummarizedExperiment::assayNames(log_transformed_se)
#'
#' @export
setGeneric(
  "logTransform",
  function(x, assay_name = "normcounts", pseudo_count = 1, base = 2) {
    standardGeneric("logTransform")
  }
)

#' @rdname logTransform
#' @export
setMethod(
  "logTransform",
  signature(x = "ScDeconvExperiment"),
  function(x, assay_name = "normcounts", pseudo_count = 1, base = 2) {
    if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
      stop("Assay not found: ", assay_name)
    }

    if (!is.numeric(pseudo_count) || length(pseudo_count) != 1 || pseudo_count < 0) {
      stop("'pseudo_count' must be a single non-negative numeric value.")
    }

    if (!is.numeric(base) || length(base) != 1 || base <= 0 || base == 1) {
      stop("'base' must be a positive numeric value different from 1.")
    }

    values <- SummarizedExperiment::assay(x, assay_name)

    logcounts <- log(values + pseudo_count, base = base)

    SummarizedExperiment::assay(x, "logcounts") <- logcounts

    x@params$log_transform <- list(
      input_assay = assay_name,
      output_assay = "logcounts",
      pseudo_count = pseudo_count,
      base = base
    )
    x@history <- c(x@history, "logTransform")
    methods::validObject(x)
    x
  }
)



#' Normalize counts to the median library size
#'
#' Normalizes raw counts by dividing each cell by its library size
#' and multiplying by the median library size across cells.
#'
#' @param x A ScDeconvExperiment object.
#'
#' @return A ScDeconvExperiment object with normalized values stored in
#'   the `"normcounts"` assay.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(
#'   rpois(1000, lambda = 5),
#'   nrow = 100,
#'   ncol = 10
#' )
#' rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#' colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = counts)
#' )
#'
#' se <- ScDeconvExperiment(sce)
#'
#' normalized_se <- normalizeMedianLibrarySize(se)
#'
#' SummarizedExperiment::assayNames(normalized_se)
#'
#' @export
setGeneric("normalizeMedianLibrarySize", function(x) {
  standardGeneric("normalizeMedianLibrarySize")
})


#' @rdname normalizeMedianLibrarySize
#' @export
setMethod(
  "normalizeMedianLibrarySize",
  signature(x = "ScDeconvExperiment"),
  function(x) {
    counts <- .get_counts(x)
    library_size <- .compute_library_size(counts)

    median_library_size <- stats::median(library_size)

    normcounts <- .scale_columns(counts, median_library_size / library_size)

    SummarizedExperiment::assay(x, "normcounts") <- normcounts
    SummarizedExperiment::colData(x)$library_size <- library_size

    x@params$normalization <- list(
      method = "median_library_size",
      target_total = median_library_size
    )
    x@history <- c(x@history, "normalizeMedianLibrarySize")
    methods::validObject(x)
    x
  }
)
