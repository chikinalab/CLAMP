# Binarize matrix by top-k values per column

Keeps only the top `top` values in each column of a matrix, setting
others to 0.

## Usage

``` r
binarizeTop(Z, top, keepVals = TRUE)
```

## Arguments

- Z:

  A numeric matrix.

- top:

  Number of top entries to keep in each column.

- keepVals:

  If `TRUE`, retains original values above the cutoff; otherwise, sets
  them to 1.

## Value

A modified matrix with only top entries retained per column.
