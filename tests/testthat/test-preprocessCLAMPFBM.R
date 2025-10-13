test_that("preprocessCLAMPFBM filters FBM and returns correct structure", {
    mat <- matrix(rnorm(100), 10, 10)
    fbm <- FBM(nrow(mat), ncol(mat))
    fbm[, ] <- mat
    out <- preprocessCLAMPFBM(fbm = fbm, mean_cutoff = 0, var_cutoff = 0)
    expect_type(out, "list")
    expect_s4_class(out$fbm_filtered, "FBM")
    expect_setequal(names(out$rowStats), c("row_means", "row_variances"))
    expect_true(is.numeric(out$kept_rows))
})

test_that("preprocessCLAMPFBM errors on non-FBM input", {
    expect_error(preprocessCLAMPFBM(fbm = matrix(1:4, 2, 2)))
})
