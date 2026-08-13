# Find the location of the maximum of a smoothing spline

Find the location of the maximum of a smoothing spline

## Usage

``` r
findSplineMax(x, y, n = 1000, spar = NULL)
```

## Arguments

- x:

  Numeric vector of x values.

- y:

  Numeric vector of y values (same length as `x`).

- n:

  Integer, number of grid points to evaluate (default 1000).

- spar:

  Smoothing parameter passed to
  [`stats::smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html).

## Value

A named list with components:

- x:

  The x coordinate at which the spline reaches its maximum.

- y:

  The corresponding maximum fitted y value.

## Examples

``` r
x <- seq_len(10)
y <- sin(x) + rnorm(10, 0, 0.1)
findSplineMax(x, y)
#> $x
#> [1] 7.963964
#> 
#> $y
#> [1] 1.15936
#> 
```
