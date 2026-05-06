# Compute all-vs-all AUC matrix

Calculates the area under the ROC curve (AUC) for all pairs of columns
between a prediction matrix `B` and binary targets in `target`. Each
column of `B` is ranked, and AUC is computed based on how well the ranks
separate positive vs. negative samples in each target column.

## Usage

``` r
allAgainstAllAUCs(B, target)
```

## Arguments

- B:

  A numeric matrix of predictions (samples × features).

- target:

  A binary matrix of the same number of rows as `B` (samples × targets),
  where 1 indicates positive and 0 indicates negative.

## Value

A numeric matrix of AUC values (features × targets).

## Examples

``` r
set.seed(1)
B <- matrix(rnorm(100), nrow = 20)
target <- matrix(sample(0:1, 40, replace = TRUE), nrow = 20)
allAgainstAllAUCs(B, target)
#>           [,1]      [,2]
#> [1,] 0.4848485 0.4895833
#> [2,] 0.5959596 0.6354167
#> [3,] 0.1919192 0.5520833
#> [4,] 0.4545455 0.5000000
#> [5,] 0.2929293 0.2500000
```
