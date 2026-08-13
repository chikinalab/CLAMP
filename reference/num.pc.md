# Estimate number of principal components via elbow or permutation method

Estimate number of principal components via elbow or permutation method

## Usage

``` r
num.pc(data, method = c("elbow", "permutation"), B = 20, seed = NULL)
```

## Arguments

- data:

  Either a matrix (e.g. z-scored data) or an SVD result (list with \$d).

- method:

  One of "elbow" (fast) or "permutation" (slower).

- B:

  Number of permutations (for method = "permutation").

- seed:

  Seed for reproducibility.

## Value

       Estimated number of PCs.

## Examples

``` r
# generate a small random matrix: 5 features x 10 samples
mat <- matrix(rnorm(5 * 10), nrow = 5)
# fast elbow estimate
num.pc(mat, method = "elbow")
#> Computing SVD
#> [1] 2
# slower permutation estimate (use fewer perms for example speed)
num.pc(mat, method = "permutation", B = 5)
#> Computing SVD
#> [1] 0
```
