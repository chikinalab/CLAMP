# Preprocess a bigstatsr FBM for CLAMP

Makes a writable copy of the input FBM, cleans it (log2 transform if
needed, fill NAs), filters rows by mean/variance, and returns the
filtered FBM plus stats and indices.

## Usage

``` r
preprocessCLAMPFBM(
  fbm,
  mean_cutoff = NULL,
  var_cutoff = NULL,
  backingfile = NULL,
  block_size = 1000,
  ncores = 1,
  log2_transform = TRUE
)
```

## Arguments

- fbm:

  A bigstatsr::FBM (genes x samples), possibly read-only.

- mean_cutoff:

  Numeric or NULL. Minimum row mean to keep (NULL = no mean filter).

- var_cutoff:

  Numeric or NULL. Minimum row variance to keep (NULL = no var filter).

- backingfile:

  Character or NULL. Base name for the *copy* FBM and filtered FBM on
  disk. If NULL, defaults to paste0(fbm\$backingfile, '\_preproc') and
  '\_filtered'.

- block_size:

  Number of rows to process at a time when copying data. Default is
  1000.

- ncores:

  Integer; number of cores to use for parallel operations (default 1).

- log2_transform:

  Logical; enable automatic `log2(x + 1)` transformation when the
  maximum value is at least 100 (default TRUE). Set to FALSE to skip
  transformation. Missing values are still replaced with zero.

## Value

A list with:

- fbm_filtered:

  The filtered FBM (writable).

- rowStats:

  List with row_means & row_variances for fbm_filtered.

- kept_rows:

  Integer vector of original row indices that were retained.

## Examples

``` r
library(bigstatsr)
# create a toy matrix and back it with an FBM
mat <- matrix(
    c(
        1, 2, 3, # geneA
        10, 20, 30, # geneB
        100, 200, 300 # geneC
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("geneA", "geneB", "geneC"), paste0("s", seq_len(3)))
)
fbm <- FBM(nrow(mat), ncol(mat), init = mat)

# preprocess without filtering (all genes kept)
res_all <- preprocessCLAMPFBM(fbm)
#> Applying log2 transformation
#> No NA values found
```
