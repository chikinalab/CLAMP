# Subset and filter multiple pathway matrices to match target genes

Filters gene-by-pathway annotation matrices to retain only pathways with
sufficient overlap with a given gene set. The result is a sparse matrix
aligned to `new.genes`, combining all inputs column-wise.

## Usage

``` r
getMatchedPathwayMat2(..., new.genes, min.genes = 10)
```

## Arguments

- ...:

  One or more sparse binary matrices (genes x pathways).

- new.genes:

  Character vector of gene names to match.

- min.genes:

  Minimum number of overlapping genes required to keep a pathway.

## Value

A sparse matrix with rows = `new.genes` and columns = filtered pathways
from all inputs.
