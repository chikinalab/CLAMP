# Project new data into CLAMP latent space

Computes the latent loadings `B` for new gene expression data using the
latent variables `Z` from a fitted CLAMP model. This allows transfer of
the learned latent structure to new datasets with matched genes.

## Usage

``` r
projectCLAMP(
  CLAMPres,
  newdata,
  scale = 1,
  ncores = 1,
  align = TRUE,
  verbose = TRUE
)
```

## Arguments

- CLAMPres:

  A result object from
  [`CLAMPfull()`](https://chikinalab.org/CLAMP/reference/CLAMPfull.md)
  or
  [`CLAMPbase()`](https://chikinalab.org/CLAMP/reference/CLAMPbase.md),
  containing at least `Z` and `L2`.

- newdata:

  A gene expression matrix (genes x samples) to be projected. Can be a
  standard matrix, sparse matrix, or FBM/big.matrix.

- scale:

  Optional numeric multiplier for the L2 regularization terms. Default
  is 1.

- ncores:

  Number of cores to use for parallel computation (only used if newdata
  is an FBM). Default is 1.

- align:

  Logical; if `TRUE` (default), row names shared by `CLAMPres$Z` and
  `newdata` are used to align both matrices to the same genes in the
  same order before projection. If row names are not available,
  dimensions must already match.

- verbose:

  Logical; if `TRUE` (default), report how many common rows are used for
  projection.

## Value

A matrix `B` of projected latent loadings (LVs x samples) for the new
dataset.

## Details

This function uses ridge-regularized least squares to compute
`B = solve(Z'Z + L2 * I) * Z'Y`, where `Z` is the latent matrix from the
trained CLAMP model and `Y` is the new dataset. If `newdata` is a
Filebacked Big Matrix (FBM) and does not need row-name subsetting, the
computation is optimized using
[`bigstatsr::big_cprodMat()`](https://privefl.github.io/bigstatsr/reference/big_cprodMat.html).

## Examples

``` r
# fit a tiny CLAMP model for projection
Y0 <- matrix(rnorm(5 * 3),
    nrow = 5,
    dimnames = list(paste0("Gene", 1:5), paste0("S", 1:3))
)
base <- CLAMPbase(Y0, clamp_k = 2, max.iter = 1, trace = FALSE)
#> ****
#> Computing SVD
#> CLAMP k is set to 2
#> L1 is set to 0.660478480751283
#> L2 is set to 1.98143544225385
# new data can be provided in a different row order
newY <- matrix(rnorm(5 * 2),
    nrow = 5,
    dimnames = list(rev(rownames(Y0)), paste0("N", 1:2))
)
projB <- projectCLAMP(base, newdata = newY)
#> 5 common rows found
# check dimensions: 2 latent vars x 2 samples
dim(projB)
#> [1] 2 2
```
