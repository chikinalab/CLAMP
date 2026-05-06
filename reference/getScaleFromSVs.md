# Estimate noise scale from singular values with linear tail extrapolation

This function estimates a characteristic scale from a vector of singular
values by fitting a linear model to the tail and extrapolating to length
`n`. The median of the extrapolated values is returned as the estimate.
If no sufficiently linear tail is detected (based on `min_r2`) or the
extrapolation produces negative values, a fallback estimate is returned
using the 75\\

## Usage

``` r
getScaleFromSVs(sv, n, min_r2 = 0.95)
```

## Arguments

- sv:

  Numeric vector of singular values sorted in decreasing order.

- n:

  Integer, total length to which the linear tail is extrapolated.

- min_r2:

  Numeric, minimum R-squared value required for accepting the linear
  tail fit. Defaults to `0.95`.

## Value

A numeric scalar giving the estimated scale. If linear extrapolation is
unreliable, returns `sv[ceiling(0.75 * length(sv))]`.

## Examples

``` r
sv <- exp(-seq(0, 5, length.out = 50)) + rnorm(50, 0, 0.01)
getScaleFromSVs(sv, n = 100)
#> $scale
#> [1] 0.01510799
#> 
```
