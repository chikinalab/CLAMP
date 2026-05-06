# Subset and filter pathway matrix to match target genes

Filters a gene-by-pathway annotation matrix to retain only pathways with
sufficient overlap with a given gene set. The result is a sparse matrix
aligned to `new.genes`, with columns (pathways) retained only if they
have at least `min.genes` matched genes.

## Usage

``` r
getMatchedPathwayMat(pathMat, new.genes, min.genes = 10)
```

## Arguments

- pathMat:

  A sparse binary matrix of genes (rows) x pathways (columns).

- new.genes:

  Character vector of gene names to match.

- min.genes:

  Minimum number of overlapping genes required to keep a pathway.

## Value

A sparse matrix of dimensions `length(new.genes)` x filtered pathways.

## Examples

``` r
library(Matrix)
# create a toy gene-by-pathway sparse matrix
genes <- paste0("g", seq_len(6))
pathways <- c("Path1", "Path2", "Path3")
# Path1: g1, g2; Path2: g2, g3, g4; Path3: g5
pathMat <- sparseMatrix(
    i = c(1, 2, 2, 3, 4, 5),
    j = c(1, 1, 2, 2, 2, 3),
    dims = c(length(genes), length(pathways)),
    dimnames = list(genes, pathways)
)
new.genes <- genes
filtered <- getMatchedPathwayMat(pathMat, new.genes, min.genes = 2)
#> There are 6 genes in the intersection between data and prior
#> Removing 1 pathways
```
