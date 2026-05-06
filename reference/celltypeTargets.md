# Cell-type deconvolution matrix

A numeric matrix of estimated cell-type proportions for whole-blood
samples. Rows correspond to sample IDs and columns to major immune cell
types. This dataset can be used for validation or illustrative purposes
in CLAMP analyses.

## Usage

``` r
data(celltypeTargets)
```

## Format

A numeric matrix with samples as rows and cell types as columns. Row
names are sample identifiers.

## Value

A numeric matrix of cell-type proportions.

## Examples

``` r
data(celltypeTargets)
```
