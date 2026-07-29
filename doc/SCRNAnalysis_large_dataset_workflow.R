## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

required_packages <- c(
  "SCRNAnalysis", "scRNAseq", "SingleCellExperiment",
  "SummarizedExperiment", "S4Vectors", "scater", "SingleR",
  "ggplot2", "ggdendro", "ggalluvial"
)
packages_available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

# ZeiselBrainData is downloaded through ExperimentHub on first use. Keep the
# package build self-contained; set this environment variable to render all
# computations and figures locally.
run_large_example <- identical(
  Sys.getenv("SCRNANALYSIS_RUN_LARGE_VIGNETTE"), "true"
)
workflow_eval <- all(packages_available) && run_large_example

## ----missing-packages, echo=FALSE, results='asis'-----------------------------
if (!all(packages_available)) {
  cat(
    "> Install the optional packages required for this workflow: ",
    paste(names(packages_available)[!packages_available], collapse = ", "),
    ".\n\n",
    sep = ""
  )
}

## ----enable-large-workflow, eval=FALSE----------------------------------------
# Sys.setenv(SCRNANALYSIS_RUN_LARGE_VIGNETTE = "true")
# devtools::build_vignettes(".")

## ----load-data, eval=workflow_eval, message=FALSE, warning=FALSE--------------
# library(SCRNAnalysis)
# library(SingleCellExperiment)
# library(SummarizedExperiment)
# library(S4Vectors)
# 
# brain <- scRNAseq::ZeiselBrainData()
# labels <- SummarizedExperiment::colData(brain)$level1class
# 
# keep <- !is.na(labels)
# brain <- brain[, keep]
# labels <- droplevels(factor(labels[keep]))
# 
# table(labels)
# dim(brain)

## ----split-data, eval=workflow_eval-------------------------------------------
# set.seed(2026)
# query_cells <- sample(seq_len(ncol(brain)), size = floor(0.4 * ncol(brain)))
# 
# query_sce <- brain[, query_cells]
# reference_sce <- brain[, -query_cells]
# SummarizedExperiment::colData(reference_sce)$cell_type <- labels[-query_cells]
# 
# # SingleR expects a log-expression assay in the reference.
# reference_sce <- scater::logNormCounts(reference_sce)
# query <- ScDeconvExperiment(query_sce)

## ----preprocess-query, eval=workflow_eval-------------------------------------
# query <- filterLowExpressionGenes(
#   query,
#   assay_name = "counts",
#   min_count = 1,
#   lower_percentile = 0.01
# )
# query <- filterCellsByDetectedGenes(
#   query,
#   assay_name = "counts",
#   min_count = 1,
#   lower_percentile = 0.01,
#   upper_percentile = 0.99
# )
# query <- normalizeToTotal(query, target_total = 1e4)
# query <- logTransform(query, assay_name = "normcounts")
# 
# n_hvg <- min(2000, nrow(query))
# query <- identifyHVGs(query, assay_name = "logcounts", n_top = n_hvg)

## ----cluster-and-annotate, eval=workflow_eval, message=FALSE, warning=FALSE----
# n_pcs <- min(30, ncol(query) - 1, sum(SummarizedExperiment::rowData(query)$is_hvg) - 1)
# n_pcs <- max(2, n_pcs)
# 
# query <- performClustering(
#   query,
#   assay_name = "logcounts",
#   algorithm = "kmeans",
#   use_hvg = TRUE,
#   n_hvg = n_hvg,
#   n_pcs = n_pcs,
#   n_clusters = 10,
#   cluster_column = "cluster"
# )
# 
# query <- performDeconvolution(
#   query,
#   reference = reference_sce,
#   flavour = "singler",
#   reference_label_column = "cell_type",
#   assay_name = "logcounts",
#   ref_assay_name = "logcounts",
#   label_column = "predicted_cell_type"
# )

## ----compare, eval=workflow_eval----------------------------------------------
# query <- compareClusterCellTypes(
#   query,
#   cluster_column = "cluster",
#   cell_type_column = "predicted_cell_type",
#   metadata_name = "large_dataset_comparison"
# )

## ----heatmap, eval=workflow_eval, message=FALSE, warning=FALSE----------------
# plotClusterCellTypeHeatmap(
#   query,
#   metadata_name = "large_dataset_comparison",
#   value = "row_proportions"
# )

## ----dendrogram, eval=workflow_eval, message=FALSE, warning=FALSE-------------
# plotClusterCellTypeDendrogram(
#   query,
#   metadata_name = "large_dataset_comparison",
#   value = "row_proportions",
#   margin = "clusters"
# )

## ----alluvial, eval=workflow_eval, message=FALSE, warning=FALSE---------------
# plotClusterCellTypeAlluvial(
#   query,
#   metadata_name = "large_dataset_comparison"
# )

## ----history, eval=workflow_eval----------------------------------------------
# query@history
# str(query@params, max.level = 2)

