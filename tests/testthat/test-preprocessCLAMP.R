test_that("preprocessCLAMP filters and returns correct structure", {
    Y <- matrix(rnorm(100), 10, 10)
    out <- preprocessCLAMP(Y = Y, mean_cutoff = 0, var_cutoff = 0)
    expect_true(is.list(out))
    expect_true(is.matrix(out$Y_filtered))
    expect_true(is.data.frame(out$rowStats))
})

test_that("preprocessCLAMP errors on non-numeric input", {
    expect_error(preprocessCLAMP(
        Y = data.frame(a = letters[1:5]),
        mean_cutoff = 0, var_cutoff = 0
    ))
})
