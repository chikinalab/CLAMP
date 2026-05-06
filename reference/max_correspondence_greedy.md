# Greedy maximum correspondence from correlation matrix

Finds a one-to-one assignment (permutation) between rows and columns of
a square correlation matrix that maximizes the total correlation score,
using a greedy algorithm.

## Usage

``` r
max_correspondence_greedy(cor_mat)
```

## Arguments

- cor_mat:

  A square numeric matrix of pairwise correlations (rows = items, cols =
  items).

## Value

A vector of assignments (integer indices).
