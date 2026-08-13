# Differential latent-variable activity between sample groups

For each row of a CLAMP `B` matrix, compares mean activity in a
reference sample group against all other samples using a Wilcoxon
rank-sum test, with Benjamini–Hochberg FDR adjustment.

## Usage

``` r
differentialLVActivity(
  x,
  metadata,
  sample_col = "id",
  group_col = "type",
  reference
)
```

## Arguments

- x:

  A CLAMP result list (with element `B`) or a numeric matrix with latent
  variables in rows and samples in columns.

- metadata:

  Data frame with sample identifiers and group labels.

- sample_col:

  Name of the column in `metadata` holding sample identifiers matching
  `colnames(B)`.

- group_col:

  Name of the column in `metadata` holding group labels.

- reference:

  Label of the reference group. All other samples are treated as the
  comparison group.

## Value

A data frame with one row per latent variable, ordered by FDR. Columns:
`LV`, `Mean_Reference`, `Mean_Comparison`, `Mean_Diff`, `P_Value`,
`FDR`.

## Examples

``` r
B <- matrix(rnorm(30), nrow = 3)
rownames(B) <- paste0("LV", seq_len(3))
colnames(B) <- paste0("S", seq_len(10))
meta <- data.frame(
    id = colnames(B),
    type = rep(c("Control", "Case"), each = 5)
)
differentialLVActivity(B, meta, reference = "Control")
#>      LV Mean_Reference Mean_Comparison  Mean_Diff   P_Value       FDR
#> LV1 LV1     0.09618778     -0.41251441 -0.5087022 0.3095238 0.4642857
#> LV3 LV3    -0.63789922     -0.13613513  0.5017641 0.3095238 0.4642857
#> LV2 LV2     0.25218352     -0.02508714 -0.2772707 0.6904762 0.6904762
```
