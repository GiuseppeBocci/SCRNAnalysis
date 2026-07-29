###########
# Helpers #
###########

utils::globalVariables(c(
  "cluster",
  "cell_type",
  "value",
  "n",
  "x",
  "y",
  "xend",
  "yend",
  "label",
  "stratum"
))

.print_plot_error <- function(message) {
  print(paste0("Error: ", message))
  invisible(NULL)
}

.get_comparison_metadata <- function(object, metadata_name = "cluster_cell_type_comparison") {
  md <- S4Vectors::metadata(object)

  if (!metadata_name %in% names(md)) {
    return(.print_plot_error(
      paste0("metadata entry '", metadata_name, "' not found. Run compareClusterCellTypes() first.")
    ))
  }

  md[[metadata_name]]
}

.get_comparison_matrix <- function(comparison, value = "row_proportions") {
  value <- match.arg(value, c("counts", "row_proportions", "column_proportions"))
  as.matrix(comparison[[value]])
}


#' Plot cluster-cell type heatmap
#'
#' @param object A ScDeconvExperiment object.
#' @param metadata_name Metadata entry containing the comparison result.
#' @param value Matrix to plot. One of "counts", "row_proportions",
#'   or "column_proportions".
#'
#' @return A ggplot object, or NULL if metadata is missing.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   counts <- matrix(
#'     rpois(1000, lambda = 5),
#'     nrow = 100,
#'     ncol = 10
#'   )
#'   rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#'   colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'
#'   sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts)
#'   )
#'
#'   se <- ScDeconvExperiment(sce)
#'
#'   SummarizedExperiment::colData(se)$cluster <- rep(
#'     c("Cluster1", "Cluster2"),
#'     each = 5
#'   )
#'   SummarizedExperiment::colData(se)$cell_type <- rep(
#'     c("TypeA", "TypeB"),
#'     times = 5
#'   )
#'
#'   comparison_se <- compareClusterCellTypes(
#'     se,
#'     cluster_column = "cluster",
#'     cell_type_column = "cell_type"
#'   )
#'
#'   plotClusterCellTypeHeatmap(
#'     comparison_se,
#'     metadata_name = "cluster_cell_type_comparison",
#'     value = "row_proportions"
#'   )
#' }
#'
#' @export
setGeneric(
  "plotClusterCellTypeHeatmap",
  function(object, metadata_name = "cluster_cell_type_comparison", value = "row_proportions") {
    standardGeneric("plotClusterCellTypeHeatmap")
  }
)

#' @rdname plotClusterCellTypeHeatmap
#' @export
setMethod(
  "plotClusterCellTypeHeatmap",
  "ScDeconvExperiment",
  function(object, metadata_name = "cluster_cell_type_comparison", value = "row_proportions") {
    comparison <- .get_comparison_metadata(object, metadata_name)
    if (is.null(comparison)) return(invisible(NULL))

    mat <- .get_comparison_matrix(comparison, value)
    plot_data <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(plot_data) <- c("cluster", "cell_type", "value")

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = cell_type, y = cluster, fill = value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = round(value, 3)), size = 3) +
      ggplot2::labs(
        title = paste0("Cluster-cell type heatmap (", value, ")"),
        x = comparison$cell_type_column,
        y = comparison$cluster_column,
        fill = value
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    # print(p)
    p
  }
)


#' Plot dendrogram from cluster-cell type comparison
#'
#' @param object A ScDeconvExperiment object.
#' @param metadata_name Metadata entry containing the comparison result.
#' @param value Matrix to use. One of "counts", "row_proportions",
#'   or "column_proportions".
#' @param margin What to cluster. One of "clusters" or "cell_types".
#' @param hclust_method Linkage method passed to stats::hclust.
#'
#' @return A ggplot object, or NULL if metadata is missing.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("ggdendro", quietly = TRUE)) {
#'   set.seed(1)
#'   counts <- matrix(
#'     rpois(1000, lambda = 5),
#'     nrow = 100,
#'     ncol = 10
#'   )
#'   rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#'   colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'
#'   sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts)
#'   )
#'
#'   se <- ScDeconvExperiment(sce)
#'
#'   SummarizedExperiment::colData(se)$cluster <- rep(
#'     c("Cluster1", "Cluster2"),
#'     each = 5
#'   )
#'   SummarizedExperiment::colData(se)$cell_type <- rep(
#'     c("TypeA", "TypeB"),
#'     times = 5
#'   )
#'
#'   comparison_se <- compareClusterCellTypes(
#'     se,
#'     cluster_column = "cluster",
#'     cell_type_column = "cell_type"
#'   )
#'
#'   plotClusterCellTypeDendrogram(
#'     comparison_se,
#'     metadata_name = "cluster_cell_type_comparison",
#'     value = "row_proportions",
#'     margin = "clusters"
#'   )
#' }
#'
#' @export
setGeneric(
  "plotClusterCellTypeDendrogram",
  function(object,
           metadata_name = "cluster_cell_type_comparison",
           value = "row_proportions",
           margin = "clusters",
           hclust_method = "ward.D2") {
    standardGeneric("plotClusterCellTypeDendrogram")
  }
)

#' @rdname plotClusterCellTypeDendrogram
#' @export
setMethod(
  "plotClusterCellTypeDendrogram",
  "ScDeconvExperiment",
  function(object,
           metadata_name = "cluster_cell_type_comparison",
           value = "row_proportions",
           margin = "clusters",
           hclust_method = "ward.D2") {
    margin <- match.arg(margin, c("clusters", "cell_types"))

    comparison <- .get_comparison_metadata(object, metadata_name)
    if (is.null(comparison)) return(invisible(NULL))

    mat <- .get_comparison_matrix(comparison, value)
    clustering_matrix <- if (margin == "clusters") mat else t(mat)

    hc <- stats::hclust(stats::dist(clustering_matrix), method = hclust_method)
    dendro_data <- ggdendro::dendro_data(hc, type = "rectangle")

    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = ggdendro::segment(dendro_data),
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
      ) +
      ggplot2::geom_text(
        data = ggdendro::label(dendro_data),
        ggplot2::aes(x = x, y = y, label = label),
        angle = 90,
        hjust = 1,
        size = 3
      ) +
      ggplot2::labs(
        title = paste0("Dendrogram of ", margin, " (", value, ")"),
        x = NULL,
        y = "Distance"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )

    # print(p)
    p
  }
)


#' Plot alluvial comparison between clusters and cell types
#'
#' @param object A ScDeconvExperiment object.
#' @param metadata_name Metadata entry containing the comparison result.
#'
#' @return A ggplot object, or NULL if metadata is missing.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("ggalluvial", quietly = TRUE)) {
#'   set.seed(1)
#'   counts <- matrix(
#'     rpois(1000, lambda = 5),
#'     nrow = 100,
#'     ncol = 10
#'   )
#'   rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
#'   colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
#'
#'   sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts)
#'   )
#'
#'   se <- ScDeconvExperiment(sce)
#'
#'   SummarizedExperiment::colData(se)$cluster <- rep(
#'     c("Cluster1", "Cluster2"),
#'     each = 5
#'   )
#'   SummarizedExperiment::colData(se)$cell_type <- rep(
#'     c("TypeA", "TypeB"),
#'     times = 5
#'   )
#'
#'   comparison_se <- compareClusterCellTypes(
#'     se,
#'     cluster_column = "cluster",
#'     cell_type_column = "cell_type"
#'   )
#'
#'   plotClusterCellTypeAlluvial(
#'     comparison_se,
#'     metadata_name = "cluster_cell_type_comparison"
#'   )
#' }
#'
#' @export
setGeneric(
  "plotClusterCellTypeAlluvial",
  function(object, metadata_name = "cluster_cell_type_comparison") {
    standardGeneric("plotClusterCellTypeAlluvial")
  }
)

#' @rdname plotClusterCellTypeAlluvial
#' @export
setMethod(
  "plotClusterCellTypeAlluvial",
  "ScDeconvExperiment",
  function(object, metadata_name = "cluster_cell_type_comparison") {
    comparison <- .get_comparison_metadata(object, metadata_name)
    if (is.null(comparison)) return(invisible(NULL))

    plot_data <- comparison$long_table
    plot_data <- plot_data[plot_data$n > 0, , drop = FALSE]

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(axis1 = cluster, axis2 = cell_type, y = n)
    ) +
      ggalluvial::geom_alluvium(
        ggplot2::aes(fill = cluster),
        width = 1 / 12
      ) +
      ggalluvial::geom_stratum(width = 1 / 12) +
      ggalluvial::stat_stratum(
        ggplot2::aes(label = ggplot2::after_stat(stratum)),
        geom = "text"
      ) +
      ggplot2::scale_x_discrete(
        limits = c("Cluster", "Cell type"),
        expand = c(0.05, 0.05)
      ) +
      ggplot2::labs(
        x = NULL,
        y = "Number of cells",
        title = "Cluster-cell type alluvial plot"
      ) +
      ggplot2::theme_minimal()

    # print(p)
    p
  }
)
