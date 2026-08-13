# Ridge regression update for B

Solves \\B = (Z^T Z + L2)^{-1} Z^T Y\\ with ridge penalty matrix L2.

## Usage

``` r
ridge_B(Y, Z, L2k)
```

## Arguments

- Y:

  Gene expression matrix (genes x samples), dense or FBM.

- Z:

  Latent variable matrix (genes x k).

- L2k:

  Ridge penalty matrix (k x k).

## Value

A numeric matrix of size k x samples.

## Examples

``` r
set.seed(123)

genes <- paste0("Gene", 1:50)
samples <- paste0("S", 1:20)
k <- 5

Y <- matrix(rnorm(50 * 20), nrow = 50, dimnames = list(genes, samples))
Z <- matrix(rnorm(50 * k),
    nrow = 50,
    dimnames = list(genes, paste0("LV", 1:k))
)

lambda <- 0.1
L2k <- diag(lambda, k)

# Solve for B = (Z'Z + L2)^(-1) Z'Y
B <- ridge_B(Y, Z, L2k)
```
