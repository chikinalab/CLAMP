# Compare two sets of factor loadings or embeddings

Compares correspondence between two result matrices (e.g., factor
loadings) across multiple targets using correlation, AUC, or
t-statistics. Produces a paired comparison plot and summary statistics.

## Usage

``` r
compareBs(
  res1,
  res2,
  target,
  method = "p",
  xlab = "1",
  ylab = "2",
  stat.method = "t",
  oneToOne = TRUE
)
```

## Arguments

- res1, :

  res2 Result objects or matrices containing factor loadings. If a list,
  must contain element `B`; if of class `"rsvd"`, `t(res$v)` is used.

- res2:

  Second result object or matrix to compare against.

- target:

  Numeric matrix of target variables (samples × targets).

- method:

  Character, one of `"p"`, `"s"`, `"a"`, or `"t"`, indicating Pearson,
  Spearman, AUC, or t-statistic comparison.

- xlab, :

  ylab Labels for x and y axes in the plot.

- ylab:

  Label for y-axis (used in plot).

- stat.method:

  Statistical test to compare correlations (`"t"` or `"wilcox"`).

- oneToOne:

  Logical, whether to apply one-to-one masking of associations.

## Value

A list with:

- plot:

  A ggplot object comparing maximal correlations across targets.

- df:

  A data.frame with per-target metrics and (if available) top genes.

## Examples

``` r
set.seed(123)
# Simulated example with 50 genes × 20 samples
Y <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)

# Run two CLAMP-like decompositions (here using simple SVD)
svd1 <- rsvd::rsvd(Y, k = 5)
svd2 <- rsvd::rsvd(Y + matrix(rnorm(50 * 20, 0, 0.1), 50, 20), k = 5)

# Define a target variable (e.g., binary or continuous trait)
target <- matrix(rnorm(20 * 3), ncol = 3)
colnames(target) <- c("Trait1", "Trait2", "Trait3")

# Compare the two sets of embeddings
res <- compareBs(svd1, svd2, target,
    method = "p", xlab = "SVD1", ylab = "SVD2"
)
#> [1]  5 20
#> [1]  5 20
```
