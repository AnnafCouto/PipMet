# PipMet (Legacy Version)

> **Note:** This repository contains the original version of the PipMet package developed during 2021-2022. 
> For the latest stable version, updated documentation, and the official release associated with the 2026 publication, please visit the official repository: [github.com/PipMet/PipMet](https://github.com/PipMet/PipMet).

## Publication
This work is part of the research published in:
*Brenelli, T. L., Couto, A. C. F., Aricetti, J.. (2026). PipMet: Pipeline for processing GC-MS Metabolomics data and statistical graphics visualization [Link] (https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.15000601/v2)

<!-- badges: start -->

The **PipMet** package was developed to perform end-to-end processing of
metabolomic-based GC-MS data, with automated generation of high-quality
figures throughout the workflow. All user inputs are obtained through
pop-up windows.  


<!-- badges: end -->

## Installation

### Dependencies installation

Before installing `PipMet`, make sure that all required dependencies are
installed.  
The following code will automatically check for and install all
necessary **Bioconductor** and **CRAN** packages:

```r
# Install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor dependencies
BiocManager::install(c(
  "xcms", "MSnbase", "CluMSID", "metaMS", "BiocParallel",
  "Biobase", "ProtGenerics", "CAMERA", "NormalyzerDE"
), ask = FALSE, update = TRUE)

# Install CRAN dependencies
cran_pkgs <- c("svDialogs", "pheatmap", "ddpcr", "webchem", "fritools", "pracma")
installed <- cran_pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(cran_pkgs[!installed])

```
### PipMet installation

You can install the released version of PipMet from
[GitHub](https://github.com) with:

``` r {eval=FALSE}
devtools::install_github("AnnafCouto/PipMet")
```

The user may also install the package through Bioconductor repository:

``` r {eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PipMet")
```

## Example

The package is constituted of two main functions with pre-set parameters
and algorithms for GC-MS data processing. The `workData()` function
reads, treats and process GC-MS sample data, with metabolite
identification and quantification. The second one, `workLib()` provides
a workflow for an internal library creation to be uploaded into NIST MS
Search software for spectra annotation.

The package was thought to be as friendly-user as possible. Therefore, when information is needed,
pop-ups will appear to collect input.

``` r {eval=FALSE}
library(PipMet)
result <- workData(
   sample_dir = system.file("extdata", package = "PipMet"),
   metadata = system.file("extdata", "metadata.csv", package = "PipMet"),
   extension = ".mzXML",
   myDir = '~/',
   example = TRUE,
   pictures = TRUE
)
```

``` r
library(PipMet)
workLib(
   extension = ".mzML",
   myDir = '~/',
   example = TRUE,
)
```

Set ‘pictures = TRUE’ to generate pictures throughout the code.

For more information, see the package
[vignette](vignettes/workData.Rmd).
