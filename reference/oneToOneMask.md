# One-to-one masking of maximum associations

Selects the highest-scoring one-to-one pairs between rows and columns of
a matrix, similar to a greedy bipartite matching. All other entries are
set to a sentinel value (default `-100`).

## Usage

``` r
oneToOneMask(cc)
```

## Arguments

- cc:

  A numeric matrix of association scores.

## Value

A numeric matrix of the same dimensions as `cc`, where only the selected
one-to-one maxima are retained and all other entries are set to `-100`.

## Examples

``` r
set.seed(1)
m <- matrix(runif(16), 4, 4)
oneToOneMask(m)
#>              [,1]         [,2]          [,3]         [,4]
#> [1,] -100.0000000 -100.0000000 -100.00000000    0.6870228
#> [2,] -100.0000000 -100.0000000    0.06178627 -100.0000000
#> [3,] -100.0000000    0.9446753 -100.00000000 -100.0000000
#> [4,]    0.9082078 -100.0000000 -100.00000000 -100.0000000
```
