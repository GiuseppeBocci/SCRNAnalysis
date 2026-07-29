# SCRNAnalysis

`SCRNAnalysis` is an R package for streamlined and reproducible single-cell RNA-seq analysis. It defines the `ScDeconvExperiment` class, extending `SingleCellExperiment`, and provides tools for:

* gene and cell quality control;
* count normalization and log transformation;
* highly variable gene selection;
* dimensionality reduction and clustering;
* cell-type annotation with Seurat label transfer or SingleR;
* comparison and visualization of clusters and predicted cell types;
* tracking workflow parameters and analysis history.

## Installation

Install the development version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github(
    "GiuseppeBocci/SCRNAnalysis",
    build_vignettes = TRUE
)
```

Then load the package:

```r
library(SCRNAnalysis)
```

## Documentation

The main vignette provides a lightweight end-to-end example using Seurat's `pbmc_small` dataset:

[Read the main workflow vignette](https://github.com/GiuseppeBocci/SCRNAnalysis/blob/main/vignettes/SCRNAnalysis_workflow_en.Rmd)

A more realistic workflow using a larger annotated mouse brain dataset is also available:

[Read the larger-dataset vignette](https://github.com/GiuseppeBocci/SCRNAnalysis/blob/main/vignettes/SCRNAnalysis_large_dataset_workflow.Rmd)

Installed vignettes can also be opened from R:

```r
browseVignettes("SCRNAnalysis")
```

## BiocCheck

The package follows a Bioconductor-oriented structure and can be assessed with [`BiocCheck`](https://bioconductor.org/packages/BiocCheck). Run `R CMD check` first, then execute `BiocCheck` from the directory containing the package:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("BiocCheck")
BiocCheck::BiocCheck("SCRNAnalysis")
```

`BiocCheck` reports Bioconductor-specific errors, warnings, and notes that may need to be addressed before submission.

## License

This project is released under the MIT License.
