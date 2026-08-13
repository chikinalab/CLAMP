# Compute counts-per-million (CPM) for CLAMP pipelines

Compute counts-per-million (CPM) for CLAMP pipelines

## Usage

``` r
cpmCLAMP(counts)
```

## Arguments

- counts:

  A numeric matrix or data.frame of raw counts (genes x samples).

## Value

A numeric matrix of CPM values (same dimensions), ready for CLAMP input.

## Examples

``` r
mat <- matrix(seq_len(12), nrow = 3)
cpmCLAMP(mat)
#>          [,1]     [,2]     [,3]     [,4]
#> [1,] 166666.7 266666.7 291666.7 303030.3
#> [2,] 333333.3 333333.3 333333.3 333333.3
#> [3,] 500000.0 400000.0 375000.0 363636.4
```
