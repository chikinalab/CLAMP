test_that("zscoreCLAMPFBM applies Z‐score in‐place on an FBM", {
    mat <- matrix(c(2, 4, 6, 8), nrow = 2, byrow = TRUE)
    fbm <- bigstatsr::FBM(nrow(mat), ncol(mat), init = mat)
    stats <- list(
        row_means     = rowMeans(mat),
        row_variances = apply(mat, 1, var)
    )

    expect_message(
        zscoreCLAMPFBM(fbm, stats, chunk_size = ncol(mat)),
    )

    out <- fbm[, ]
    expected <- sweep(
        sweep(mat, 1, stats$row_means, FUN = "-"),
        1,
        sqrt(stats$row_variances),
        FUN = "/"
    )

    expect_equal(out, expected)
})
