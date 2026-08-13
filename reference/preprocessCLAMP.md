# Preprocess an expression matrix for CLAMP

Cleans an expression matrix, filters genes by mean expression and
variance, and returns the filtered matrix and per-gene statistics. To
match
[`preprocessCLAMPFBM()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMPFBM.md),
values are transformed with `log2(Y + 1)` when the maximum value is at
least 100, and missing values are replaced with zero.

## Usage

``` r
preprocessCLAMP(Y, mean_cutoff = 0, var_cutoff = 0)
```

## Arguments

- Y:

  Numeric matrix of gene expression (rows = genes, cols = samples)

- mean_cutoff:

  Numeric. Minimum row-mean required to keep a gene (default 0).

- var_cutoff:

  Numeric. Minimum row-variance required to keep a gene (default 0).

## Value

A list with components:

- Y_filtered: filtered matrix (genes x samples)

- rowStats: data.frame with columns mean and variance for each kept gene

- kept_rows: integer vector of the original row indices that were kept

## Examples

``` r
# construct a small example matrix
mat <- matrix(
    c(
        1, 5, 10,
        2, 6, 11,
        3, 7, 12,
        4, 8, 13
    ),
    nrow = 4, byrow = FALSE,
    dimnames = list(paste0("gene", seq_len(4)), paste0("sample", seq_len(3)))
)

# keep genes with mean >= 6 and variance >= 2
res <- preprocessCLAMP(mat, mean_cutoff = 6, var_cutoff = 2)
#> Already on log scale or all NA
#> No NA values found
```
