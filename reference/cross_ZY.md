# Cross-product Z^T Y with FBM or dense matrices

Computes \\Z^T Y\\, handling FBM objects from bigstatsr as well as base
R matrices.

## Usage

``` r
cross_ZY(Y, Z)
```

## Arguments

- Y:

  Gene expression matrix (genes x samples), dense or FBM.

- Z:

  Latent variable matrix (genes x k).

## Value

A numeric matrix giving Z^T Y.

## Examples

``` r
set.seed(123)

genes <- 40
samples <- 10
k <- 4

Y <- matrix(rnorm(genes * samples), nrow = genes)
Z <- matrix(rnorm(genes * k), nrow = genes)

# Compute Z^T Y
res1 <- cross_ZY(Y, Z)
dim(res1) # k × samples
#> [1]  4 10
```
