# Subset and filter pathway matrix to match target genes

Filters a gene-by-pathway annotation matrix to retain only pathways with
sufficient overlap with a given gene set. The result is a sparse matrix
aligned to `new.genes`, with columns (pathways) retained only if they
have at least `min.genes` matched genes.

## Usage

``` r
getMatchedPathwayMatOld(pathMat, new.genes, min.genes = 10)
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
