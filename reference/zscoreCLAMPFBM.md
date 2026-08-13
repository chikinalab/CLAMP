# Z-score a filtered FBM in-place

Standardizes each row of an FBM using provided row means and variances.

## Usage

``` r
zscoreCLAMPFBM(fbm_filtered, rowStats, chunk_size = 1000, ncores = 1)
```

## Arguments

- fbm_filtered:

  A bigstatsr::FBM produced by preprocessCLAMPFBM().

- rowStats:

  A list with row_means and row_variances from that FBM.

- chunk_size:

  Columns per block (default 1000).

- ncores:

  Integer; number of cores to use for parallel operations (default 1).

## Value

A normalized FBM with z-scored rows.

## Examples

``` r
library(bigstatsr)
fbm <- FBM(
    nrow = 2, ncol = 4,
    init = matrix(seq_len(8), nrow = 2)
)
stats <- list(
    row_means     = rowMeans(fbm[]),
    row_variances = apply(fbm[], 1, var)
)
zscoreCLAMPFBM(fbm, stats, chunk_size = 2)
#> Applying Z-score transformation
```
