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

test_that("preprocessCLAMP log-transforms and fills missing values", {
    Y <- matrix(
        c(0, 3, NA, 100, 127, 255),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("gene1", "gene2"), paste0("sample", 1:3))
    )
    expected <- log2(Y + 1)
    expected[is.na(expected)] <- 0

    expect_message(
        out <- preprocessCLAMP(Y, mean_cutoff = 0, var_cutoff = 0),
        "Applying log2 transformation"
    )
    expect_equal(out$Y_filtered, expected)
    expect_equal(out$rowStats$mean, unname(rowMeans(expected)))
    expect_equal(
        out$rowStats$variance,
        unname(apply(expected, 1, stats::var))
    )
})

test_that("preprocessCLAMP leaves log-scale input unchanged", {
    Y <- matrix(seq(0, 5), nrow = 2)

    expect_message(
        out <- preprocessCLAMP(Y, mean_cutoff = 0, var_cutoff = 0),
        "Already on log scale"
    )
    expect_equal(out$Y_filtered, Y)
})

test_that("preprocessCLAMP log-transforms at the threshold", {
    Y <- matrix(c(0, 99, 1, 100), nrow = 2)

    expect_message(
        out <- preprocessCLAMP(Y, mean_cutoff = 0, var_cutoff = 0),
        "Applying log2 transformation"
    )
    expect_equal(out$Y_filtered, log2(Y + 1))
})
