# Whole-blood reference expression matrix

A numeric matrix of whole-blood gene expression where rows correspond to
genes and columns to samples.

## Usage

``` r
data(dataWholeBlood)
```

## Format

A numeric matrix with G genes (rows) and N samples (columns). Row names
are gene symbols, and column names are sample IDs.

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130824>

## Value

A numeric matrix of expression values.

## Details

This object is a whole-blood RNA-seq expression matrix. The source data
are publicly available from NCBI GEO under accession
[GSE130824](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130824)
(Homo sapiens whole-blood RNA-seq, 36 samples). Raw sequencing data are
deposited in SRA, and the processed normalized expression file is
provided as `GSE130824_dataNormedFiltered.txt.gz`. The matrix bundled
with CLAMP is derived from that processed file and serves as a compact
example dataset for package demonstrations and unit tests.

## Examples

``` r
data(dataWholeBlood)
```
