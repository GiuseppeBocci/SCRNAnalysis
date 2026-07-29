#' Identify highly variable genes
#'
#' Identifies highly variable genes using gene-level variability statistics.
#'
#' @param object A ScDeconvExperiment object.
#' @param assay_name Name of the assay to use. Default is "logcounts".
#' @param n_top Number of highly variable genes to select.
#' @param method Variability score to use. One of "variance", "dispersion", or "cv2".
#' @param min_mean Minimum average expression required for a gene to be eligible.
#' @param max_mean Maximum average expression allowed for a gene to be eligible.
#' @param store_metrics Logical. Whether to store HVG metrics in rowData.
#'
#' @return The ScDeconvExperiment object with HVG information stored in rowData.
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
#' logcounts <- log1p(counts)
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = counts,
#'     logcounts = logcounts
#'   ),
#'   colData = S4Vectors::DataFrame(
#'     cell_type = rep(c("A", "B"), each = 5)
#'   ),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(counts)
#'   )
#' )
#'
#' se <- ScDeconvExperiment(sce)
#'
#' se <- identifyHVGs(
#'   se,
#'   assay_name = "logcounts",
#'   n_top = 20,
#'   method = "dispersion"
#' )
#'
#' head(SummarizedExperiment::rowData(se)$is_hvg)
#' head(SummarizedExperiment::rowData(se)$hvg_rank)
#' head(SummarizedExperiment::rowData(se)$hvg_score)
#'
#' @export
setGeneric(
  "identifyHVGs",
  function(object,
           assay_name = "logcounts",
           n_top = 2000,
           method = "dispersion",
           min_mean = 0,
           max_mean = Inf,
           store_metrics = TRUE) {
    standardGeneric("identifyHVGs")
  }
)

#' @rdname identifyHVGs
#' @export
setMethod(
  "identifyHVGs",
  "ScDeconvExperiment",
  function(object,
           assay_name = "logcounts",
           n_top = 2000,
           method = "dispersion",
           min_mean = 0,
           max_mean = Inf,
           store_metrics = TRUE) {

    if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
      stop("'assay_name' not found in object.")
    }
    method <- match.arg(method,choices = c("variance", "dispersion", "cv2"))
    if (!is.numeric(n_top) || length(n_top) != 1 || n_top <= 0) {
      stop("'n_top' must be a positive number.")
    }
    if (min_mean < 0) {
      stop("'min_mean' must be non-negative.")
    }
    if (max_mean <= min_mean) {
      stop("'max_mean' must be greater than 'min_mean'.")
    }

    log_counts <- SummarizedExperiment::assay(object, assay_name)
    log_counts <- as.matrix(log_counts)
    gene_mean <- Matrix::rowMeans(log_counts)


    gene_variance <- apply(log_counts, 1, stats::var)

    gene_dispersion <- rep(NA_real_, length(gene_mean))
    valid_mean <- gene_mean > 0
    gene_dispersion[valid_mean] <- gene_variance[valid_mean] / gene_mean[valid_mean]

    gene_cv2 <- rep(NA_real_, length(gene_mean))
    gene_cv2[valid_mean] <- gene_variance[valid_mean] / gene_mean[valid_mean]^2


    if (method == "variance") {
      hvg_score <- gene_variance
    } else if (method == "dispersion") {
      hvg_score <- gene_dispersion
    } else {
      hvg_score <- gene_cv2
    }

    eligible <- gene_mean >= min_mean & gene_mean <= max_mean & is.finite(hvg_score)

    eligible_indices <- which(eligible)

    ordered_indices <- eligible_indices[order(hvg_score[eligible_indices], decreasing = TRUE)]

    n_select <- min(as.integer(n_top), length(ordered_indices))
    selected_indices <- ordered_indices[seq_len(n_select)]

    is_hvg <- rep(FALSE, nrow(object))
    is_hvg[selected_indices] <- TRUE
    hvg_rank <- rep(NA_integer_, nrow(object))
    hvg_rank[ordered_indices] <- seq_along(ordered_indices)

    if (store_metrics) {
      SummarizedExperiment::rowData(object)$hvg_mean <- gene_mean
      SummarizedExperiment::rowData(object)$hvg_variance <- gene_variance
      SummarizedExperiment::rowData(object)$hvg_dispersion <- gene_dispersion
      SummarizedExperiment::rowData(object)$hvg_cv2 <- gene_cv2
      SummarizedExperiment::rowData(object)$hvg_score <- hvg_score
      SummarizedExperiment::rowData(object)$hvg_rank <- hvg_rank
      SummarizedExperiment::rowData(object)$is_hvg <- is_hvg
    }

    object@params$identifyHVGs <- list(
      assay_name = assay_name,
      n_top = n_top,
      n_selected = sum(is_hvg),
      method = method,
      min_mean = min_mean,
      max_mean = max_mean
    )

    object@history <- c(
      object@history,
      paste0(
        "Identified ",
        sum(is_hvg),
        " highly variable genes using method '",
        method,
        "'."
      )
    )

    methods::validObject(object)
    object
  }
)


# Internal helper: select rows for clustering
.select_clustering_rows <- function(object,
                                    use_hvg = TRUE,
                                    hvg_column = "is_hvg",
                                    hvg_rank_column = "hvg_rank",
                                    n_hvg = NULL) {
  if (!use_hvg) {
    return(seq_len(nrow(object)))
  }

  rd <- SummarizedExperiment::rowData(object)

  if (!hvg_column %in% colnames(rd)) {
    stop("HVG column not found in rowData: ", hvg_column)
  }

  selected <- which(rd[[hvg_column]])

  if (length(selected) == 0) {
    stop("No HVGs found in rowData column: ", hvg_column)
  }

  if (!is.null(n_hvg)) {
    if (hvg_rank_column %in% colnames(rd)) {
      ranks <- rd[[hvg_rank_column]][selected]
      selected <- selected[order(ranks, na.last = NA)]
    }

    selected <- utils::head(selected, n_hvg)
  }

  selected
}


#' Perform clustering of cells
#'
#' Performs cell clustering using PCA followed by hierarchical, k-means,
#' Louvain, or Leiden clustering.
#'
#' @param object A ScDeconvExperiment object.
#' @param assay_name Name of the assay to use. Default is `"logcounts"`.
#' @param algorithm Clustering algorithm. One of `"louvain"`, `"leiden"`,
#'   `"kmeans"`, or `"hierarchical"`.
#' @param use_hvg Logical. Whether to use highly variable genes.
#' @param hvg_column Name of the rowData column marking HVGs.
#' @param hvg_rank_column Name of the rowData column containing HVG ranks.
#' @param n_hvg Optional number of HVGs to use.
#' @param n_pcs Number of principal components.
#' @param n_clusters Number of clusters for k-means and hierarchical clustering.
#' @param k Number of nearest neighbors for graph-based clustering.
#' @param resolution Resolution parameter for Louvain or Leiden clustering.
#' @param cluster_column Name of the colData column where clusters are stored.
#' @param pca_name Name of the reducedDim slot where PCA coordinates are stored.
#'
#' @return The ScDeconvExperiment object with cluster labels stored in colData.
#'
#' @examples
#' if (requireNamespace("scater", quietly = TRUE)) {
#'   set.seed(1)
#'   counts <- matrix(
#'     rpois(1000, lambda = 5),
#'     nrow = 100,
#'     ncol = 10
#'   )
#'   rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#'   colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'   
#'   logcounts <- log1p(counts)
#'   
#'   sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(
#'       counts = counts,
#'       logcounts = logcounts
#'     ),
#'     colData = S4Vectors::DataFrame(
#'       cell_type = rep(c("A", "B"), each = 5)
#'     ),
#'     rowData = S4Vectors::DataFrame(
#'       gene_id = rownames(counts)
#'     )
#'   )
#'   
#'   se <- ScDeconvExperiment(sce)
#'
#'   se <- performClustering(
#'     se,
#'     assay_name = "logcounts",
#'     algorithm = "kmeans",
#'     use_hvg = FALSE,
#'     n_pcs = 5,
#'     n_clusters = 2
#'   )
#'
#'   head(SummarizedExperiment::colData(se)$cluster)
#' }
#'
#' @export
setGeneric(
  "performClustering",
  function(object,
           assay_name = "logcounts",
           algorithm = c("louvain", "leiden", "kmeans", "hierarchical"),
           use_hvg = TRUE,
           hvg_column = "is_hvg",
           hvg_rank_column = "hvg_rank",
           n_hvg = NULL,
           n_pcs = 20,
           n_clusters = 10,
           k = 20,
           resolution = 1,
           cluster_column = "cluster",
           pca_name = "PCA") {
    standardGeneric("performClustering")
  }
)

#' @rdname performClustering
#' @export
setMethod(
  "performClustering",
  "ScDeconvExperiment",
  function(object,
           assay_name = "logcounts",
           algorithm = c("louvain", "leiden", "kmeans", "hierarchical"),
           use_hvg = TRUE,
           hvg_column = "is_hvg",
           hvg_rank_column = "hvg_rank",
           n_hvg = NULL,
           n_pcs = 20,
           n_clusters = 10,
           k = 20,
           resolution = 1,
           cluster_column = "cluster",
           pca_name = "PCA") {

    algorithm <- match.arg(algorithm)

    selected_rows <- .select_clustering_rows(
      object = object,
      use_hvg = use_hvg,
      hvg_column = hvg_column,
      hvg_rank_column = hvg_rank_column,
      n_hvg = n_hvg
    )

    object <- scater::runPCA(
      object,
      exprs_values = assay_name,
      subset_row = selected_rows,
      ncomponents = n_pcs,
      name = pca_name
    )

    embedding <- SingleCellExperiment::reducedDim(object, pca_name)

    if (algorithm == "hierarchical") {
      n_clusters <- min(n_clusters, ncol(object))
      hc <- stats::hclust(stats::dist(embedding),method = "ward.D2")
      clusters <- stats::cutree(hc, k = n_clusters)

    } else if (algorithm == "kmeans") {
      n_clusters <- min(n_clusters, ncol(object))
      km <- stats::kmeans(
        embedding,
        centers = n_clusters
      )
      clusters <- km$cluster

    } else if (algorithm == "louvain") {
      graph <- scran::buildSNNGraph(
        object,
        use.dimred = pca_name,
        k = k
      )

      community <- igraph::cluster_louvain(
        graph,
        resolution = resolution
      )
      clusters <- igraph::membership(community)
    } else if (algorithm == "leiden") {
      graph <- scran::buildSNNGraph(
        object,
        use.dimred = pca_name,
        k = k
      )

      community <- igraph::cluster_leiden(
        graph,
        resolution_parameter = resolution
      )
      clusters <- igraph::membership(community)
    }

    SummarizedExperiment::colData(object)[[cluster_column]] <- factor(clusters)

    object@params$performClustering <- list(
      assay_name = assay_name,
      algorithm = algorithm,
      use_hvg = use_hvg,
      hvg_column = hvg_column,
      hvg_rank_column = hvg_rank_column,
      n_hvg = n_hvg,
      n_pcs = n_pcs,
      n_clusters = length(unique(clusters)),
      k = k,
      resolution = resolution,
      cluster_column = cluster_column,
      pca_name = pca_name
    )

    object@history <- c(
      object@history,
      paste0(
        "Performed ",
        algorithm,
        " clustering using ",
        length(selected_rows),
        " genes."
      )
    )

    methods::validObject(object)
    object
  }
)
