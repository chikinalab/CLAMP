# Clean a Filebacked Big Matrix (FBM) by log-transforming and handling NAs

This function inspects an FBM to determine if log-transformation is
needed (based on value range) and whether NA values are present. If the
maximum value is \>= 100, it applies a log2(x + 1) transformation
in-place. If any NA values are detected, they are replaced with 0.

## Usage

``` r
cleanFBM(fbm, ncores = 1)
```

## Arguments

- fbm:

  A `bigmemory::FBM` or
  [`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
  object.

- ncores:

  Integer; number of cores to use for parallel operations (default 1).

## Value

A list with:

- max_value:

  The maximum value encountered in the FBM (after log transformation if
  applied).

- had_na:

  Logical indicating whether any NA values were found and filled.

## Details

Modifies the FBM in place. Uses
[`bigstatsr::big_apply()`](https://privefl.github.io/bigstatsr/reference/big_apply.html)
to process in parallel-safe chunks.

## Examples

``` r
fbm <- bigstatsr::FBM(3, 4, init = matrix(
    c(0, 1, 2, NA, 100, 200, 300, 400, 5, 6, 7, 8),
    nrow = 3
))
cleanFBM(fbm, ncores = 1)
#> Applying log2 transformation
#> Filling NAs with 0
#> $max_value
#> [1] 400
#> 
#> $had_na
#> [1] TRUE
#> 
```
