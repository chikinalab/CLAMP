# Ridge-regularized pseudoinverse via SVD

Computes a stable pseudoinverse of a symmetric positive semi-definite
matrix using singular value decomposition (SVD) and ridge
regularization. This is useful when the matrix is ill-conditioned or
rank-deficient.

## Usage

``` r
pinv.ridge(m, alpha = 0)
```

## Arguments

- m:

  A symmetric numeric matrix (e.g., from
  [`crossprod()`](https://rdrr.io/pkg/Matrix/man/matmult-methods.html)).

- alpha:

  Non-negative scalar specifying the ridge penalty. A small positive
  value stabilizes the inversion by shrinking large singular values.

## Value

A numeric matrix representing the ridge-regularized pseudoinverse of
`m`.
