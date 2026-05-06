# Read a GMT file into a list

Parses a local GMT file and returns a list of gene sets. Each gene set
is represented as a character vector of unique gene names.

## Usage

``` r
read_gmt(filename)
```

## Arguments

- filename:

  A `character(1)` string giving the path to a .gmt file.

## Value

A `list` where each element is a character vector of gene names, named
by the gene set ID.

## Examples

``` r
# Bioconductor requires runnable examples.
# We create a dummy GMT file for this example:
gmt_file <- tempfile(fileext = ".gmt")
writeLines(
    c(
        "PATHWAY_A\thttp://link.com\tGENE1\tGENE2\tGENE3",
        "PATHWAY_B\thttp://link.com\tGENE2\tGENE4"
    ),
    con = gmt_file
)

# Run the function
gs <- read_gmt(gmt_file)

# Inspect results
length(gs)
#> [1] 2
names(gs)
#> [1] "PATHWAY_A" "PATHWAY_B"
gs[["PATHWAY_A"]]
#> [1] "GENE1" "GENE2" "GENE3"
```
