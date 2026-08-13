# Convert a list of GMT gene sets to a sparse matrix

Converts a list of named gene sets (e.g., from
[`getGMT()`](https://chikinalab.org/CLAMP/reference/getGMT.md)) into a
sparse binary matrix where rows are genes, columns are gene sets, and
entries are 1 if the gene is in the set.

## Usage

``` r
gmtListToSparseMat(gmtList)
```

## Arguments

- gmtList:

  A nested list of gene sets. Outer names are gene set names; each entry
  is a character vector of gene names.

## Value

A sparse binary matrix with genes as rows and gene sets as columns.

## Examples

``` r
# define a simple nested GMT list
gmt1 <- list(
    PathwayA = c("Gene1", "Gene2", "Gene3"),
    PathwayB = c("Gene2", "Gene4")
)
gmt2 <- list(
    PathwayC = c("Gene1", "Gene4"),
    PathwayD = c("Gene3", "Gene5")
)
# combine into a nested list
nestedList <- list(gmt1 = gmt1, gmt2 = gmt2)
# convert to sparse matrix
sparseMat <- gmtListToSparseMat(nestedList)
```
