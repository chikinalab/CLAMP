# PLIER2 <img src="man/figures/plier2.png" width="121px" height="140px" align="right" style="padding-left:10px;background-color:white;" />

<!-- badges: start -->
[![GitHub issues](https://img.shields.io/github/issues/mchikina/PLIER2)](https://github.com/mchikina/mchikina/PLIER2)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check-bioc](https://github.com/mchikina/PLIER2/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/mchikina/PLIER2/actions)
[![R-CMD-check](https://github.com/chikinalab/PLIER2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chikinalab/PLIER2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Bioconductor release status

|      Branch      |    R CMD check   | Last updated |
|:----------------:|:----------------:|:------------:|
| [_devel_](http://bioconductor.org/packages/devel/bioc/html/PLIER2.html) | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/PLIER2.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/PLIER2) | ![](http://bioconductor.org/shields/lastcommit/devel/bioc/PLIER2.svg) |
| [_release_](http://bioconductor.org/packages/release/bioc/html/PLIER2.html) | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/PLIER2.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/PLIER2) | ![](http://bioconductor.org/shields/lastcommit/release/bioc/PLIER2.svg) |

The goal of PLIER2 is to provide an easy-to-use package to extract interpretable latent variables from large transcriptomic datasets using biological priors.

## Local development via Conda

We keep a fully specified environment file at `envs/plier2.yaml`. From your package root, create and activate it like so:

```bash
conda env create -f envs/plier2.yaml

conda activate plier2

library(remotes)

REPO_PATH <- "~/path/to/PLIER2"  # adjust

remotes::install_local(
  REPO_PATH,
  force        = TRUE,
  dependencies = FALSE
)

library(PLIER2)
packageVersion("PLIER2")
```

## Installation

You can install the latest release of `PLIER2` from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("PLIER2")
```

If you want to test the development version, you can install it from the github repository:

``` r
BiocManager::install("mchikina/PLIER2")
```

Now you can load the package using:

``` r
library(PLIER2)
```

## Basic usage

For detailed instructions on how to use PLIER2, please see the [vignette](https://mchikina.github.io/PLIER2/articles/PLIER2.html).

``` r
library(PLIER2)
#some example
```

## Code of Conduct

Please note that the PLIER2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
