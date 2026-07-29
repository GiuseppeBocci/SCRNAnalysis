make_test_object <- function() {
  counts <- matrix(
    c(
      10, 0, 5, 0,
      3,  1, 0, 0,
      0,  4, 8, 1,
      7,  0, 0, 2,
      2,  2, 2, 2
    ),
    nrow = 5,
    byrow = TRUE
  )
  rownames(counts) <- c("MT-Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
  colnames(counts) <- paste0("Cell", 1:4)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )
  ScDeconvExperiment(sce)
}

make_comparison_object <- function() {
  object <- make_test_object()
  SummarizedExperiment::colData(object)$cluster <- factor(c("A", "A", "B", "B"))
  SummarizedExperiment::colData(object)$cell_type <- factor(c("T", "T", "B", "T"))
  object
}
