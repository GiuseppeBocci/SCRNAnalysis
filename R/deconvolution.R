#' Perform cell type annotation
#'
#' Assigns a cell type to each cell using either SingleR or Seurat.
#'
#' The reference object must already contain known cell type labels in its
#' metadata. For SingleR, the reference should be a SummarizedExperiment or
#' SingleCellExperiment object. For Seurat, the reference should be a Seurat
#' object.
#'
#' @param object A ScDeconvExperiment object.
#' @param reference Annotated reference object.
#' @param flavour Annotation method. Either "singler" or "seurat".
#' @param reference_label_column Column in the reference metadata containing
#'   known cell type labels.
#' @param assay_name Assay from object used by SingleR. Default is "logcounts".
#' @param ref_assay_name Assay from reference used by SingleR. Default is
#'   "logcounts".
#' @param label_column Column name in colData where predicted labels are stored.
#' @param dims Dimensions used by Seurat label transfer.
#'
#' @return The ScDeconvExperiment object with predicted cell types in colData.
#'
#' @examples
#' if (requireNamespace("SingleR", quietly = TRUE)) {
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
#'     )
#'   )
#'
#'   se <- ScDeconvExperiment(sce)
#'
#'   ref_counts <- matrix(
#'     rpois(2000, lambda = 5),
#'     nrow = 100,
#'     ncol = 20
#'   )
#'   rownames(ref_counts) <- rownames(counts)
#'   colnames(ref_counts) <- paste0("RefCell", seq_len(ncol(ref_counts)))
#'
#'   ref_logcounts <- log1p(ref_counts)
#'
#'   ref_colData <- S4Vectors::DataFrame(
#'     cell_type = rep(c("TypeA", "TypeB"), each = 10)
#'   )
#'
#'   reference <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(logcounts = ref_logcounts),
#'     colData = ref_colData
#'   )
#'
#'   annotated_se <- performDeconvolution(
#'     object = se,
#'     reference = reference,
#'     flavour = "singler",
#'     reference_label_column = "cell_type"
#'   )
#'
#'   SummarizedExperiment::colData(annotated_se)$cell_type
#' }
#'
#' @export
setGeneric(
  "performDeconvolution",
  function(object,
           reference,
           flavour = "singler",
           reference_label_column = "cell_type",
           assay_name = "logcounts",
           ref_assay_name = "logcounts",
           label_column = "cell_type",
           dims = 1:30) {
    standardGeneric("performDeconvolution")
  }
)

#' @rdname performDeconvolution
#' @export
setMethod(
  "performDeconvolution",
  "ScDeconvExperiment",
  function(object,
           reference,
           flavour = "singler",
           reference_label_column = "cell_type",
           assay_name = "logcounts",
           ref_assay_name = "logcounts",
           label_column = "cell_type",
           dims = 1:30) {

    if (!flavour %in% c("singler", "seurat")) {
      stop("'flavour' must be either 'singler' or 'seurat'.")
    }

    score_column <- paste0(label_column, "_score")

    if (flavour == "singler") {
      reference_labels <- SummarizedExperiment::colData(reference)[[
        reference_label_column
      ]]

      pred <- SingleR::SingleR(
        test = object,
        ref = reference,
        labels = reference_labels,
        assay.type.test = assay_name,
        assay.type.ref = ref_assay_name
      )

      assigned_labels <- pred$labels
      assigned_labels[is.na(assigned_labels)] <- "Unknown"

      prediction_score <- apply(
        pred$scores,
        1,
        max,
        na.rm = TRUE
      )

      raw_result <- pred
    }

    if (flavour == "seurat") {
      query <- Seurat::CreateSeuratObject(
        counts = SummarizedExperiment::assay(object, "counts")
      )
      query <- Seurat::NormalizeData(
        query,
        verbose = FALSE
      )
      query <- Seurat::FindVariableFeatures(
        query,
        verbose = FALSE
      )
      query <- Seurat::ScaleData(
        query,
        verbose = FALSE
      )
      query <- Seurat::RunPCA(
        query,
        verbose = FALSE
      )

      anchors <- Seurat::FindTransferAnchors(
        reference = reference,
        query = query,
        dims = dims,
        verbose = FALSE
      )

      refdata <- reference[[reference_label_column]][, 1]

      pred <- Seurat::TransferData(
        anchorset = anchors,
        refdata = refdata,
        dims = dims,
        verbose = FALSE
      )
      pred <- as.data.frame(pred)

      if (all(colnames(object) %in% rownames(pred))) {
        pred <- pred[colnames(object), , drop = FALSE]
      }

      assigned_labels <- pred$predicted.id
      assigned_labels[is.na(assigned_labels)] <- "Unknown"
      prediction_score <- pred$prediction.score.max

      raw_result <- pred
    }

    SummarizedExperiment::colData(object)[[label_column]] <- factor(
      assigned_labels
    )

    SummarizedExperiment::colData(object)[[score_column]] <- prediction_score

    md <- S4Vectors::metadata(object)

    md$deconvolution <- list(
      flavour = flavour,
      reference_label_column = reference_label_column,
      label_column = label_column,
      score_column = score_column,
      result = raw_result
    )

    S4Vectors::metadata(object) <- md

    object@params$performDeconvolution <- list(
      flavour = flavour,
      reference_label_column = reference_label_column,
      label_column = label_column,
      score_column = score_column,
      cells_annotated = ncol(object)
    )
    object@history <- c(
      object@history,
      paste0(
        "Assigned cell types using ",
        flavour,
        ". Labels stored in colData column '",
        label_column,
        "'."
      )
    )

    methods::validObject(object)
    object
  }
)


#' Compare clusters with assigned cell types
#'
#' Creates a contingency table comparing clustering labels with assigned
#' cell types. The result is stored in metadata and printed to screen.
#'
#' @param object A ScDeconvExperiment object.
#' @param cluster_column Column in colData containing cluster labels.
#' @param cell_type_column Column in colData containing assigned cell types.
#' @param metadata_name Name used to store the result in metadata.
#' @param include_na Logical. Whether to include NA values in the table.
#'
#' @return The ScDeconvExperiment object with the comparison tables stored in metadata.
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
#' SummarizedExperiment::colData(se)$cluster <- rep(
#'   c("Cluster1", "Cluster2"),
#'   each = 5
#' )
#' SummarizedExperiment::colData(se)$cell_type <- rep(
#'   c("TypeA", "TypeB"),
#'   times = 5
#' )
#'
#' comparison_se <- compareClusterCellTypes(
#'   se,
#'   cluster_column = "cluster",
#'   cell_type_column = "cell_type"
#' )
#'
#' S4Vectors::metadata(
#'   comparison_se
#' )$cluster_cell_type_comparison$long_table
#'
#' @export
setGeneric(
  "compareClusterCellTypes",
  function(object,
           cluster_column = "cluster",
           cell_type_column = "cell_type",
           metadata_name = "cluster_cell_type_comparison",
           include_na = FALSE) {
    standardGeneric("compareClusterCellTypes")
  }
)

#' @rdname compareClusterCellTypes
#' @export
setMethod(
  "compareClusterCellTypes",
  "ScDeconvExperiment",
  function(object,
           cluster_column = "cluster",
           cell_type_column = "cell_type",
           metadata_name = "cluster_cell_type_comparison",
           include_na = FALSE) {

    cd <- SummarizedExperiment::colData(object)

    if (!cluster_column %in% colnames(cd)) {
      stop("'cluster_column' not found in colData.")
    }
    if (!cell_type_column %in% colnames(cd)) {
      stop("'cell_type_column' not found in colData.")
    }

    clusters <- cd[[cluster_column]]
    cell_types <- cd[[cell_type_column]]

    use_na <- if (include_na) "ifany" else "no"

    count_table <- table(
      cluster = clusters,
      cell_type = cell_types,
      useNA = use_na
    )

    row_proportions <- prop.table(
      count_table,
      margin = 1
    )
    column_proportions <- prop.table(
      count_table,
      margin = 2
    )

    long_table <- as.data.frame(
      count_table,
      stringsAsFactors = FALSE
    )

    colnames(long_table) <- c(
      "cluster",
      "cell_type",
      "n"
    )

    long_table$row_fraction <- as.vector(row_proportions)
    long_table$column_fraction <- as.vector(column_proportions)

    long_table <- long_table[
      order(
        long_table$cluster,
        -long_table$n,
        long_table$cell_type
      ),
      ,
      drop = FALSE
    ]

    rownames(long_table) <- NULL

    md <- S4Vectors::metadata(object)

    md[[metadata_name]] <- list(
      counts = count_table,
      row_proportions = row_proportions,
      column_proportions = column_proportions,
      long_table = long_table,
      cluster_column = cluster_column,
      cell_type_column = cell_type_column
    )

    S4Vectors::metadata(object) <- md

    print(long_table)

    object@params$compareClusterCellTypes <- list(
      cluster_column = cluster_column,
      cell_type_column = cell_type_column,
      metadata_name = metadata_name
    )
    object@history <- c(
      object@history,
      paste0(
        "Compared clusters from column '",
        cluster_column,
        "' with cell types from column '",
        cell_type_column,
        "'."
      )
    )

    methods::validObject(object)
    object
  }
)
