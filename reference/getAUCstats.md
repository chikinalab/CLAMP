# Count number of latent variables exceeding AUC thresholds

Given a summary data frame from cross-validation, reports the number of
latent variables with maximum AUC exceeding 0.7, 0.8, and 0.9.

## Usage

``` r
getAUCstats(summary)
```

## Arguments

- summary:

  A data frame with `LV index` and `AUC` columns.

## Value

A named numeric vector with counts for thresholds 0.7, 0.8, and 0.9.
