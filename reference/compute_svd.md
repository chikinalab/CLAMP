# Compute a truncated SVD for a CLAMP input matrix

Dispatches to the appropriate SVD backend based on the class of `Y`:
[`bigstatsr::big_SVD`](https://privefl.github.io/bigstatsr/reference/big_SVD.html)
for `FBM` objects,
[`irlba::irlba`](https://rdrr.io/pkg/irlba/man/irlba.html) for sparse
`dgCMatrix` objects, and
[`rsvd::rsvd`](https://rdrr.io/pkg/rsvd/man/rsvd.html) otherwise. Used
by the CLAMP solvers so that the SVD step is handled in one place.

## Usage

``` r
compute_svd(Y, k = NULL)
```

## Arguments

- Y:

  A matrix-like object (dense matrix, `dgCMatrix`, or `FBM`).

- k:

  Integer number of components to compute. If `NULL` (the default),
  `select_svd_k(Y)` is used.

## Value

A list with `d`, `u`, `v` components (structure depends on the backend
but these three fields are always present).

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(100), nrow = 20, ncol = 5)
res <- compute_svd(Y, k = 3)
length(res$d)
#> [1] 3
```
