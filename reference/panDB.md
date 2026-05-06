# panDB gene-set database

A list of curated pathway and biological process gene sets used as prior
knowledge for CLAMP and other latent-variable models.

## Usage

``` r
data(panDB)
```

## Format

A named list of length M, where each element is a gene-set collection.

## Value

A list of gene-set collections used as priors for pathway-informed
modeling.

## Details

Each element corresponds to a functional collection (e.g., KEGG,
Reactome, GO), where each entry contains a character vector of gene
symbols.

## Examples

``` r
data(panDB)
```
