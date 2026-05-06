# Subset and filter multiple pathway matrices to match target genes

Filters one or more gene-by-pathway annotation matrices to retain only
pathways with sufficient overlap with a given target gene set. Each
input matrix is restricted to `new.genes`, and pathways with fewer than
`min.genes` overlapping genes are removed. The resulting matrices are
column-bound into a single sparse matrix aligned to `new.genes`.

## Usage

``` r
getMatchedPathwayMatList(..., new.genes, min.genes = 10)
```

## Arguments

- ...:

  One or more binary matrices (genes x pathways), either base `matrix`
  or sparse `Matrix` objects. Row names must be gene identifiers.

- new.genes:

  Character vector of gene names to align all pathway matrices to.

- min.genes:

  Integer; minimum number of overlapping genes required for a pathway to
  be retained. Default is 10.

## Value

A sparse binary matrix with rows equal to `new.genes` and columns equal
to the union of filtered pathways from all input matrices.

## Examples

``` r
set.seed(123)
library(Matrix)

# Simulate two small pathway matrices (genes × pathways)
genes <- paste0("Gene", 1:100)
pathways1 <- paste0("Path", 1:5)
pathways2 <- paste0("Path", 6:10)

mat1 <- matrix(sample(c(0, 1), 100 * 5, replace = TRUE, prob = c(0.9, 0.1)),
    nrow = 100, ncol = 5,
    dimnames = list(genes, pathways1)
)
mat2 <- matrix(sample(c(0, 1), 100 * 5, replace = TRUE, prob = c(0.9, 0.1)),
    nrow = 100, ncol = 5,
    dimnames = list(genes, pathways2)
)

# Define target genes (subset of total)
new.genes <- sample(genes, 50)

# Match and filter pathways with at least 5 genes
matched <- getMatchedPathwayMatList(mat1, mat2,
    new.genes = new.genes, min.genes = 5
)
#> There are 50 genes in the intersection between data and prior
#> Removing 4 pathways
#> There are 50 genes in the intersection between data and prior
#> Removing 4 pathways
```
