## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

required_packages <- c(
  "SCRNAnalysis",
  "Seurat",
  "SeuratObject",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "scater",
  "ggplot2",
  "ggdendro",
  "ggalluvial"
)

packages_available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
workflow_eval <- all(packages_available)

## ----missing-packages, echo=FALSE, results='asis'-----------------------------
if (!workflow_eval) {
  missing_packages <- names(packages_available)[!packages_available]
  cat(
    "> Some chunks are not evaluated because the following packages are missing: ",
    paste(missing_packages, collapse = ", "),
    ".\n\n",
    sep = ""
  )
}

## ----packages, eval=workflow_eval, message=FALSE, warning=FALSE---------------
library(SCRNAnalysis)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(S4Vectors)

## ----load-seurat-data, eval=workflow_eval, message=FALSE, warning=FALSE-------
pbmc_small <- SeuratObject::pbmc_small
reference <- Seurat::UpdateSeuratObject(pbmc_small)
reference

## ----inspect-reference-labels, eval=workflow_eval-----------------------------
head(reference[[]])
table(reference$letter.idents)

## ----preprocess-reference, eval=workflow_eval, message=FALSE, warning=FALSE----
reference <- Seurat::NormalizeData(reference, verbose = FALSE)
reference <- Seurat::FindVariableFeatures(reference, verbose = FALSE)
reference <- Seurat::ScaleData(reference, verbose = FALSE)
reference <- Seurat::RunPCA(reference, verbose = FALSE)

## ----create-object, eval=workflow_eval, message=FALSE, warning=FALSE----------
get_seurat_assay_data <- function(x, assay = SeuratObject::DefaultAssay(x), layer = "counts") {
  tryCatch(
    SeuratObject::GetAssayData(x, assay = assay, layer = layer),
    error = function(e) SeuratObject::GetAssayData(x, assay = assay, slot = layer)
  )
}

counts <- get_seurat_assay_data(reference, assay = "RNA", layer = "counts")
query_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
query <- ScDeconvExperiment(query_sce)
query

dim(SummarizedExperiment::assay(query, "counts"))

## ----quality-control, eval=workflow_eval--------------------------------------
query <- filterLowExpressionGenes(
  query,
  assay_name = "counts",
  min_count = 1,
  lower_percentile = 0.01
)
query <- filterCellsByDetectedGenes(
  query,
  assay_name = "counts",
  min_count = 1,
  lower_percentile = 0.01,
  upper_percentile = 0.99
)
query@params$filterLowExpressionGenes
query@params$filterCellsByDetectedGenes
dim(SummarizedExperiment::assay(query, "counts"))

## ----qc-metadata, eval=workflow_eval------------------------------------------
head(SummarizedExperiment::rowData(query))
head(SummarizedExperiment::colData(query))

## ----normalization, eval=workflow_eval----------------------------------------
query <- normalizeToTotal(query, target_total = 1e4)
query <- logTransform(query, assay_name = "normcounts", pseudo_count = 1, base = 2)

SummarizedExperiment::assayNames(query)
head(SummarizedExperiment::colData(query)$library_size)

## ----hvg, eval=workflow_eval--------------------------------------------------
n_hvg <- min(100, nrow(query))
query <- identifyHVGs(query, assay_name = "logcounts", n_top = n_hvg, method = "dispersion")

table(SummarizedExperiment::rowData(query)$is_hvg)
head(SummarizedExperiment::rowData(query)[order(SummarizedExperiment::rowData(query)$hvg_rank), ])

## ----clustering, eval=workflow_eval, message=FALSE, warning=FALSE-------------
n_pcs <- min(5, ncol(query) - 1, sum(SummarizedExperiment::rowData(query)$is_hvg) - 1)
n_pcs <- max(2, n_pcs)

query <- performClustering(
  object = query,
  assay_name = "logcounts",
  algorithm = "kmeans",
  use_hvg = TRUE,
  n_hvg = n_hvg,
  n_pcs = n_pcs,
  n_clusters = 7,
  cluster_column = "cluster",
  pca_name = "PCA"
)

table(SummarizedExperiment::colData(query)$cluster)

## ----seurat-annotation, eval=workflow_eval, message=FALSE, warning=FALSE------
dims <- seq_len(min(10, ncol(SeuratObject::Embeddings(reference, "pca"))))

query <- performDeconvolution(
  object = query,
  reference = reference,
  flavour = "seurat",
  reference_label_column = "letter.idents",
  label_column = "seurat_label",
  dims = dims
)

table(SummarizedExperiment::colData(query)$seurat_label)
summary(SummarizedExperiment::colData(query)$seurat_label_score)

## ----cluster-cell-type-comparison, eval=workflow_eval-------------------------
query <- compareClusterCellTypes(
  object = query,
  cluster_column = "cluster",
  cell_type_column = "seurat_label",
  metadata_name = "cluster_cell_type_comparison"
)

comparison <- S4Vectors::metadata(query)$cluster_cell_type_comparison
head(comparison$long_table)

## ----heatmap, eval=workflow_eval, message=FALSE, warning=FALSE----------------
plotClusterCellTypeHeatmap(query, metadata_name = "cluster_cell_type_comparison", value = "row_proportions")

## ----dendrogram, eval=workflow_eval, message=FALSE, warning=FALSE-------------
plotClusterCellTypeDendrogram(query, metadata_name = "cluster_cell_type_comparison", value = "row_proportions", margin = "clusters")

## ----alluvial, eval=workflow_eval, message=FALSE, warning=FALSE---------------
plotClusterCellTypeAlluvial(query, metadata_name = "cluster_cell_type_comparison")

## ----workflow-history, eval=workflow_eval-------------------------------------
query@history
str(query@params, max.level = 2)

## ----singler-alternative, eval=FALSE------------------------------------------
# reference_sce <- Seurat::as.SingleCellExperiment(reference)
# 
# if (!"logcounts" %in% SummarizedExperiment::assayNames(reference_sce)) {
#   ref_logcounts <- get_seurat_assay_data(reference, assay = "RNA", layer = "data")
#   SummarizedExperiment::assay(reference_sce, "logcounts") <- ref_logcounts
# }
# 
# SummarizedExperiment::colData(reference_sce)$cell_type <- reference$letter.idents
# 
# query <- performDeconvolution(
#   object = query,
#   reference = reference_sce,
#   flavour = "singler",
#   reference_label_column = "cell_type",
#   assay_name = "logcounts",
#   ref_assay_name = "logcounts",
#   label_column = "singler_label"
# )
# 
# table(SummarizedExperiment::colData(query)$singler_label)

