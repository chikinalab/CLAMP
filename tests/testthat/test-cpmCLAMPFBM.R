test_that("cpmCLAMPFBM computes CPM in‐place on an FBM", {
    mat <- matrix(c(1, 3, 2, 4), nrow = 2)
    fbm <- bigstatsr::FBM(nrow(mat), ncol(mat), init = mat)
    lib <- colSums(mat)
    expected <- sweep(mat, 2, lib, "/") * 1e6

    cpmCLAMPFBM(fbm, block_size = 1)
    out <- fbm[, ]
    expect_equal(out, expected)
})
