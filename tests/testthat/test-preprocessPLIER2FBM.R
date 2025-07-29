context("preprocessPLIER2FBM")

test_that("preprocessPLIER2FBM filters FBM and returns correct structure", {
    mat <- matrix(rnorm(100), 10, 10)
    fbm <- FBM(nrow(mat), ncol(mat))
    fbm[, ] <- mat
    out <- preprocessPLIER2FBM(fbm = fbm, mean_cutoff = 0, var_cutoff = 0)
    expect_type(out, "list")
    expect_s4_class(out$fbm_filtered, "FBM")
    expect_setequal(names(out$rowStats), c("row_means", "row_variances"))
    expect_true(is.numeric(out$kept_rows))
})

test_that("preprocessPLIER2FBM errors on non-FBM input", {
    expect_error(preprocessPLIER2FBM(fbm = matrix(1:4, 2, 2)))
})
