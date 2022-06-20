
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PipMet

<!-- badges: start -->

The PipMet package was developed to process begin-to-end
metabolomic-based GC-MS data, with the automatized generation of
pictures throughout the processing in high quality. Every input from the
user is taken through pop-up windows.This package was firstly developed
for the Brazilian National Biorenewables Laboratory (LNBR) from the
Brazilian National Center for Research in Energy and Materials (CNPEM).

<!-- badges: end -->

## Installation

You can install the released version of PipMet from
[GitHub](https://github.com) with:

``` r
devtools::install_github("AnnafCouto/PipMet")
```

## Example - Data processing

The package is constituted of one main function with pre-set parameters
and algorithms. To start the processing:

``` r
library(PipMet)
#> Warning in fun(libname, pkgname): mzR has been built against a different Rcpp version (1.0.7)
#> than is installed on your system (1.0.8). This might lead to errors
#> when loading mzR. If you encounter such issues, please send a report,
#> including the output of sessionInfo() to the Bioc support forum at 
#> https://support.bioconductor.org/. For details see also
#> https://github.com/sneumann/mzR/wiki/mzR-Rcpp-compiler-linker-issue.
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
#> 
#> Attaching package: 'PipMet'
#> The following object is masked from 'package:stats':
#> 
#>     heatmap
```
The package was thought to be the more friendly-user as possible, even though it runs in code lines. Therefore, when information is needed, pop-ups will appear so that the user can input the information in them.

``` r
result <- GC_dataProcess(
   sample_dir = system.file("extdata", package = "PipMet"),
   metadata = system.file("extdata", "metadata.csv", package = "PipMet"),
   extensao = ".mzML",
   myDir = '~/',
   example = TRUE,
   pictures = TRUE
)
```
To run your own samples, set 'example' to FALSE. You can leave sample_dir, metadata, myDir and extension empty (NULL) so the package will ask you in pop-ups.

## Example - Internal library building

The package provides an option to build your own internal library of standards from GC-MS files. The user must provide or fill a .csv table with the compound's informations such as name, retention time, filename, formula, exact mass and identifiers such as CAS, ChemSpider ou InChIKey numbers. 

``` r
GC_databaseProcess(
  sample_dir = system.file("extdata", package = "PipMet"),
  lib_metadata = system.file("extdata", "lib_metadata.csv", package = "PipMet"),
  extensao = ".mzXML",
  myDir = "~/",
  example = TRUE
}
```
