# Z-score a filtered expression matrix for CLAMP

Centers each gene to mean 0 and scales to unit variance.

## Usage

``` r
zscoreCLAMP(Y_filtered, rowStats)
```

## Arguments

- Y_filtered:

  Numeric matrix (genes x samples) returned by preprocessCLAMP

- rowStats:

  Data frame with numeric columns `mean` and `variance`, row-named to
  match `rownames(Y_filtered)`

## Value

Numeric matrix of the same dimensions as Y_filtered, with each row
centered and scaled.

## Examples

``` r
# simple 2 genes x 3 samples matrix
Y <- matrix(
    c(
        2, 4, 6, # gene1 counts
        8, 10, 12 # gene2 counts
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("gene1", "gene2"), paste0("sample", seq_len(3)))
)

# compute per‐gene mean and variance
rowStats <- data.frame(
    mean = rowMeans(Y),
    variance = apply(Y, 1, var),
    row.names = rownames(Y)
)

# z‐score each row
Y_z <- zscoreCLAMP(Y, rowStats)
```
