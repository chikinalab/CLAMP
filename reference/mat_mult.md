# Matrix multiplication with support for FBM objects

Multiplies two matrices, using optimized multiplication if the first is
a Filebacked Big Matrix (FBM).

## Usage

``` r
mat_mult(mat1, mat2, ncores = 1)
```

## Arguments

- mat1:

  A matrix or an object of class `FBM`.

- mat2:

  A numeric matrix.

- ncores:

  Number of cores to use for parallel computation (only used if mat1 is
  an FBM). Default is 1.

## Value

Matrix product of `mat1` and `mat2`.

## Examples

``` r
set.seed(123)
mat1 <- matrix(rnorm(20), nrow = 4)
mat2 <- matrix(rnorm(15), nrow = 5)
res1 <- mat_mult(mat1, mat2)
res1
#>            [,1]       [,2]       [,3]
#> [1,]  0.6717269  1.1164009 -0.1310362
#> [2,]  1.4777363 -0.8350005 -2.5216939
#> [3,] -3.0540345 -0.5431782  1.7125552
#> [4,] -1.1756616 -3.7501369  1.9061432
```
