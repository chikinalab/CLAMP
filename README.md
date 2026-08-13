# CLAMP <img src="man/figures/clamp.png" alt="CLAMP logo" width="121px" height="140px" align="right" style="padding-left:10px;background-color:white;" />

<!-- badges: start -->
[![GitHub issues](https://img.shields.io/github/issues/chikinalab/CLAMP)](https://github.com/chikinalab/CLAMP/issues)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![BiocCheck](https://github.com/chikinalab/CLAMP/actions/workflows/BiocCheck.yaml/badge.svg)](https://github.com/chikinalab/CLAMP/actions/workflows/BiocCheck.yaml)
[![R-CMD-check](https://github.com/chikinalab/CLAMP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chikinalab/CLAMP/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Bioconductor release status

|      Branch      |    R CMD check   | Last updated |
|:----------------:|:----------------:|:------------:|
| [_devel_](https://bioconductor.org/packages/devel/bioc/html/CLAMP.html) | [![Bioconductor-devel Build Status](https://bioconductor.org/shields/build/devel/bioc/CLAMP.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/CLAMP) | ![Bioconductor-devel last commit](https://bioconductor.org/shields/lastcommit/devel/bioc/CLAMP.svg) |
| [_release_](https://bioconductor.org/packages/release/bioc/html/CLAMP.html) | [![Bioconductor-release Build Status](https://bioconductor.org/shields/build/release/bioc/CLAMP.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/CLAMP) | ![Bioconductor-release last commit](https://bioconductor.org/shields/lastcommit/release/bioc/CLAMP.svg) |

The goal of CLAMP (**C**urated **L**atent-variable **A**nalysis with **M**olecular **P**riors) is to provide an easy-to-use package to extract interpretable latent variables from large transcriptomic datasets using biological priors.

## Installation

You can install the latest release of `CLAMP` from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("CLAMP")
```

If you want to test the development version, you can install it from the github repository:

``` r
BiocManager::install("chikinalab/CLAMP")
```

Now you can load the package using:

``` r
library(CLAMP)
```

## Basic usage

Full documentation is available at [chikinalab.org/CLAMP](https://chikinalab.org/CLAMP/). For detailed instructions on how to use CLAMP, please see the [vignette](https://chikinalab.github.io/CLAMP/articles/get_started.html).

``` r
library(CLAMP)

set.seed(1)

# Load example dataset (genes × samples)
data("dataWholeBlood")

# CPM-normalize, filter low-expressed genes, log2 and z-score normalizations
dataWholeBlood_cpm <- cpmCLAMP(dataWholeBlood)

prep <- preprocessCLAMP(
    dataWholeBlood_cpm,
    mean_cutoff = 0.5,
    var_cutoff = 0.1
)

Y_z <- zscoreCLAMP(
    prep$Y_filtered,
    prep$rowStats
)

# Compute truncated SVD and select number of latent variables
svd_k <- select_svd_k(Y_z)

svd_res <- compute_svd(
    Y_z,
    k = svd_k
)

clamp_k <- select_clamp_k(
    svd_res,
    n_samples = ncol(Y_z),
    svd_k = svd_k
)

# Run CLAMPbase to initialize latent variables
baseRes <- CLAMPbase(
    Y = Y_z,
    svdres = svd_res,
    clamp_k = clamp_k
)

# Load pre-fetched gene set libraries bundled with the package
gmtList <- list(
    CellMarkers = readRDS(
        system.file("extdata", "CellMarker_2024.rds", package = "CLAMP")
    ),
    KEGG = readRDS(
        system.file("extdata", "KEGG_2021_Human.rds", package = "CLAMP")
    )
)

# Combine gene set libraries into a single sparse matrix
pathMatCell <- gmtListToSparseMat(gmtList)

# Load additional xCell reference matrix
data("xCell")

# Match pathways to the filtered gene space used by CLAMP
matchedPathsWB <- getMatchedPathwayMatList(
    pathMatCell,
    xCell,
    new.genes = rownames(Y_z),
    min.genes = 2
)

# Run CLAMPfull using CLAMPbase initialization and matched priors
fullRes <- CLAMPfull(
    Y = Y_z,
    priorMat = matchedPathsWB,
    clamp.base.result = baseRes,
    svdres = svd_res,
    clamp_k = clamp_k,
    use_cpp = TRUE
)

# Inspect outputs
dim(fullRes$Z)
dim(fullRes$B)
dim(fullRes$U)

# Inspect significant pathway annotations
summary_df <- as.data.frame(fullRes$summary)

sig_summary <- summary_df[
    summary_df$FDR < 0.05 & summary_df$AUC > 0.7,
    c("LV", "pathway", "FDR", "AUC")
]

sig_summary <- sig_summary[order(sig_summary$FDR), ]

head(sig_summary, 20)
```

## Code of Conduct

Please note that the CLAMP project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
