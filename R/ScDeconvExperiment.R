#' ScDeconvExperiment class
#'
#' S4 class extending SingleCellExperiment with analysis parameters
#' and workflow history for single cell RNA-seq clustering and
#' deconvolution-based annotation.
#'
#' @slot params A list of analysis parameters.
#' @slot history A character vector describing performed analysis steps.
#'
#' @importFrom methods setClass setGeneric setMethod signature new validObject is
#' @importFrom SummarizedExperiment assayNames
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @exportClass ScDeconvExperiment
setClass(
  "ScDeconvExperiment",
  contains = "SingleCellExperiment",
  slots = list(
    params = "list",
    history = "character"
  ),
  prototype = list(
    params = list(),
    history = character()
  ),
  validity = function(object) {
    errors <- character()

    if (!"counts" %in% SummarizedExperiment::assayNames(object)) {
      errors <- c(errors, "The object must contain a 'counts' assay.")
    }
    if (nrow(object) == 0 || ncol(object) == 0) {
      errors <- c(errors, "The object must not be empty.")
    }
    if (length(errors) > 0) {
      return(errors)
    }

    return(TRUE)
  }
)



#' Create a ScDeconvExperiment object
#'
#' Creates a ScDeconvExperiment object from a SingleCellExperiment object.
#' The input object must contain a `"counts"` assay.
#'
#' @param sce A SingleCellExperiment object containing a `"counts"` assay.
#' @param params A list of analysis parameters. Defaults to an empty list.
#' @param history A character vector describing the analysis steps already
#'   performed. Defaults to an empty character vector.
#'
#' @return A ScDeconvExperiment object.
#'
#' @examples
#' counts <- matrix(
#'   c(1, 0, 3, 4, 2, 0),
#'   nrow = 3,
#'   ncol = 2
#' )
#' rownames(counts) <- c("GeneA", "GeneB", "GeneC")
#' colnames(counts) <- c("Cell1", "Cell2")
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = counts)
#' )
#'
#' obj <- ScDeconvExperiment(sce)
#' obj
#'
#' @export
ScDeconvExperiment <- function(sce, params = list(), history = character()) {
  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("'sce' must be a SingleCellExperiment object.")
  }

  object <- methods::new(
    "ScDeconvExperiment",
    sce,
    params = params,
    history = history
  )

  methods::validObject(object)

  object
}
