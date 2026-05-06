# Select default number of components for a CLAMP solver SVD

Returns the default number of components to compute in a truncated SVD
for the given input matrix. Used by the CLAMP solvers when `svd_k` is
not provided explicitly. Other SVD contexts use their own heuristics.

## Usage

``` r
select_svd_k(Y)
```

## Arguments

- Y:

  A matrix-like object (dense matrix, `dgCMatrix`, or `FBM`).

## Value

An integer: `max(2, floor((min(nrow(Y), ncol(Y)) - 1) / 4))`.

## Examples

``` r
select_svd_k(matrix(0, nrow = 100, ncol = 20))
#> [1] 4
```
