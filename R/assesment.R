#' Filter genes with very low expression
#'
#' Removes genes expressed in too few cells.
#'
#' @param object A ScDeconvExperiment object.
#' @param assay_name Name of the assay to use. Default is "counts".
#' @param min_count Minimum count required to consider a gene expressed
#'   in a cell.
#' @param lower_percentile Percentile used to define the minimum number
#'   of cells in which a gene should be detected.
#' @param min_cells Optional absolute minimum number of cells. If NULL,
#'   it is estimated from lower_percentile.
#' @param store_metrics Logical. Whether to store gene-level metrics in rowData.
#'
#' @return The filtered ScDeconvExperiment object.
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
#' filtered_se <- filterLowExpressionGenes(
#'   se,
#'   min_count = 1,
#'   lower_percentile = 0.01
#' )
#'
#' filtered_se@params$filterLowExpressionGenes
#'
#' @export
setGeneric(
  "filterLowExpressionGenes",
  function(object,
           assay_name = "counts",
           min_count = 1,
           lower_percentile = 0.01,
           min_cells = NULL,
           store_metrics = TRUE) {
    standardGeneric("filterLowExpressionGenes")
  }
)

#' @rdname filterLowExpressionGenes
#' @export
setMethod(
  "filterLowExpressionGenes",
  "ScDeconvExperiment",
  function(object,
           assay_name = "counts",
           min_count = 1,
           lower_percentile = 0.01,
           min_cells = NULL,
           store_metrics = TRUE) {

    if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
      stop("'assay_name' not found in object.")
    }
    if (lower_percentile < 0 || lower_percentile > 1) {
      stop("'lower_percentile' must be between 0 and 1.")
    }

    counts <- SummarizedExperiment::assay(object, assay_name)
    detected_cells <- Matrix::rowSums(counts >= min_count)
    total_counts <- Matrix::rowSums(counts)
    mean_counts <- Matrix::rowMeans(counts)

    if (is.null(min_cells)) {
      min_cells <- stats::quantile(
        detected_cells,
        probs = lower_percentile,
        na.rm = TRUE,
        names = FALSE
      )
      min_cells <- max(1, ceiling(min_cells))
    }

    genes_to_keep <- detected_cells >= min_cells

    if (store_metrics) {
      SummarizedExperiment::rowData(object)$detected_cells <- detected_cells
      SummarizedExperiment::rowData(object)$total_counts <- total_counts
      SummarizedExperiment::rowData(object)$mean_counts <- mean_counts
      SummarizedExperiment::rowData(object)$gene_filter_keep <- genes_to_keep
    }

    n_before <- nrow(object)
    object <- object[genes_to_keep, , drop = FALSE]
    n_after <- nrow(object)

    object@params$filterLowExpressionGenes <- list(
      assay_name = assay_name,
      min_count = min_count,
      lower_percentile = lower_percentile,
      min_cells = min_cells,
      genes_before = n_before,
      genes_after = n_after
    )

    object@history <- c(
      object@history,
      paste0(
        "Filtered low-expression genes: ",
        n_before - n_after,
        " removed, ",
        n_after,
        " retained."
      )
    )

    methods::validObject(object)
    object
  }
)


#' Filter low-quality cells by detected genes
#'
#' Removes cells with too few or too many detected genes using
#' percentile-based thresholds.
#'
#' @param object A ScDeconvExperiment object.
#' @param assay_name Name of the assay to use. Default is "counts".
#' @param min_count Minimum count required to consider a gene expressed.
#' @param lower_percentile Lower percentile for detected genes per cell.
#' @param upper_percentile Upper percentile for detected genes per cell.
#' @param store_metrics Logical. Whether to store cell-level metrics in colData.
#'
#' @return The filtered ScDeconvExperiment object.
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
#' filtered_se <- filterCellsByDetectedGenes(
#'   se,
#'   min_count = 1,
#'   lower_percentile = 0.01,
#'   upper_percentile = 0.99
#' )
#'
#' filtered_se@params$filterCellsByDetectedGenes
#'
#' @export
setGeneric(
  "filterCellsByDetectedGenes",
  function(object,
           assay_name = "counts",
           min_count = 1,
           lower_percentile = 0.01,
           upper_percentile = 0.99,
           store_metrics = TRUE) {
    standardGeneric("filterCellsByDetectedGenes")
  }
)

#' @rdname filterCellsByDetectedGenes
#' @export
setMethod(
  "filterCellsByDetectedGenes",
  "ScDeconvExperiment",
  function(object,
           assay_name = "counts",
           min_count = 1,
           lower_percentile = 0.01,
           upper_percentile = 0.99,
           store_metrics = TRUE) {

    if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
      stop("'assay_name' not found in object.")
    }
    if (lower_percentile < 0 || lower_percentile > 1 ||
        upper_percentile < 0 || upper_percentile > 1) {
      stop("Percentiles must be between 0 and 1.")
    }
    if (lower_percentile >= upper_percentile) {
      stop("'lower_percentile' must be smaller than 'upper_percentile'.")
    }

    counts <- SummarizedExperiment::assay(object, assay_name)

    if (is.null(dim(counts))) {
      counts <- matrix(
        counts,
        nrow = nrow(object),
        ncol = ncol(object),
        dimnames = list(rownames(object), colnames(object))
      )
    }
    detected_genes <- Matrix::colSums(counts >= min_count)
    total_counts <- Matrix::colSums(counts)

    lower_threshold <- stats::quantile(
      detected_genes,
      probs = lower_percentile,
      na.rm = TRUE,
      names = FALSE
    )
    upper_threshold <- stats::quantile(
      detected_genes,
      probs = upper_percentile,
      na.rm = TRUE,
      names = FALSE
    )

    cells_to_keep <- detected_genes >= lower_threshold & detected_genes <= upper_threshold

    if (store_metrics) {
      SummarizedExperiment::colData(object)$detected_genes <- detected_genes
      SummarizedExperiment::colData(object)$total_counts <- total_counts
      SummarizedExperiment::colData(object)$detected_genes_filter_keep <- cells_to_keep
    }

    n_before <- ncol(object)
    object <- object[, cells_to_keep, drop = FALSE]
    n_after <- ncol(object)

    object@params$filterCellsByDetectedGenes <- list(
      assay_name = assay_name,
      min_count = min_count,
      lower_percentile = lower_percentile,
      upper_percentile = upper_percentile,
      lower_threshold = lower_threshold,
      upper_threshold = upper_threshold,
      cells_before = n_before,
      cells_after = n_after
    )

    object@history <- c(
      object@history,
      paste0(
        "Filtered cells by detected genes: ",
        n_before - n_after,
        " removed, ",
        n_after,
        " retained."
      )
    )

    methods::validObject(object)
    object
  }
)
