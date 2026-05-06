# Winsorize matrix columns by capping the top-k values

For each column, replaces values greater than the k-th largest entry
with that threshold.

## Usage

``` r
winsor_topk(M, k)
```

## Arguments

- M:

  A numeric matrix.

- k:

  Integer; number of top elements to cap. Must be \>= 1 and \<= nrow(M).

## Value

A numeric matrix of the same dimensions as `M`, winsorized per column.

## Examples

``` r
set.seed(123)
M <- matrix(rnorm(20 * 5, mean = 0, sd = 2),
    nrow = 20,
    dimnames = list(paste0("Gene", 1:20), paste0("S", 1:5))
)

# Display column maxima before winsorization
apply(M, 2, max)
#>       S1       S2       S3       S4       S5 
#> 3.573826 2.507630 4.337912 4.100169 4.374666 

# Winsorize each column by capping top 3 values
M_winsor <- winsor_topk(M, k = 3)
```
