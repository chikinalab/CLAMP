# Compute CPM on a file-backed matrix for CLAMP (in-place)

Compute CPM on a file-backed matrix for CLAMP (in-place)

## Usage

``` r
cpmCLAMPFBM(fbm_counts, block_size = 1000, ncores = 1)
```

## Arguments

- fbm_counts:

  A bigstatsr::FBM of raw counts (genes x samples).

- block_size:

  Integer; columns per block (default 1000).

- ncores:

  Integer; number of cores to use for parallel operations (default 1).

## Value

Invisibly returns the modified FBM (now holding CPM values).

## Examples

``` r
library(bigstatsr)
mat <- matrix(c(10, 20, 30, 40, 50, 60),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), paste0("sample", seq_len(3)))
)
fbm <- FBM(nrow(mat), ncol(mat), init = mat)
cpmCLAMPFBM(fbm, block_size = 1)
```
