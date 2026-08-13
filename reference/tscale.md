# Row-wise scaling (mean 0, sd 1)

Standardizes each row of a numeric matrix to have mean 0 and standard
deviation 1. Missing values are ignored in the computation of the mean
and standard deviation.

## Usage

``` r
tscale(x)
```

## Arguments

- x:

  A numeric matrix. Each row will be scaled independently.

## Value

A numeric matrix of the same dimensions as `x`, where each row has mean
0 and standard deviation 1 (ignoring `NA`s). If a row has zero variance,
it is returned unchanged.

## See also

[`base::scale()`](https://rdrr.io/r/base/scale.html)

## Examples

``` r
mat <- matrix(seq_len(9), nrow = 3)
tscale(mat)
#>           [,1] [,2]     [,3]
#> [1,] -1.224745    0 1.224745
#> [2,] -1.224745    0 1.224745
#> [3,] -1.224745    0 1.224745
```
