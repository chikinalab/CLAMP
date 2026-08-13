# Get maximum AUC per latent variable

Summarizes a cross-validation results data frame to extract the highest
AUC value associated with each latent variable (LV).

## Usage

``` r
getMaxAUC(summary, verbose = FALSE)
```

## Arguments

- summary:

  A data frame (e.g., from
  [`crossVal()`](https://chikinalab.org/CLAMP/reference/crossVal.md)) .

- verbose:

  Logical; if `TRUE`, prints counts of LVs exceeding AUC thresholds.
  Default is `FALSE`.

## Value

A data frame with columns LV index and max_AUC.
