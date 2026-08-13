# Squash extreme z-scores

Squash extreme z-scores

## Usage

``` r
squashZscore(zdata, maxScore = 2)
```

## Arguments

- zdata:

  Numeric vector or matrix of z-scores.

- maxScore:

  Numeric scalar, maximum absolute score (default 2).

## Value

A numeric object of same dimensions as `zdata`, with values shrunk by a
hyperbolic tangent transformation.

## Examples

``` r
z <- rnorm(10, 0, 5)
squashZscore(z)
#>  [1]  0.8842618 -1.9666024  0.1065201 -1.9890084 -1.9821643  1.9995401
#>  [7]  1.9182935 -1.9985126 -1.0657078  0.6332484
```
