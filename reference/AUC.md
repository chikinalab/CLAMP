# Compute AUC using Wilcoxon rank-sum test

Computes the area under the ROC curve (AUC) by applying a Wilcoxon
rank-sum test between predicted values for positive and negative labels.
This is equivalent to computing the Mann-Whitney U statistic.

## Usage

``` r
AUC(labels, values)
```

## Arguments

- labels:

  A numeric or logical vector indicating class labels. Values \> 0 are
  treated as positive.

- values:

  A numeric vector of prediction scores corresponding to `labels`.

## Value

A list with:

- `auc`:

  Estimated AUC, or 0.5 if one class is missing

- `pval`:

  Wilcoxon test p-value, or `NA` if one class is missing
