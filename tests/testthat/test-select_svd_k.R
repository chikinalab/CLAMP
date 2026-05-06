test_that("select_svd_k returns floor((min(nrow, ncol) - 1) / 4)", {
    Y <- matrix(0, nrow = 100, ncol = 20)
    expect_equal(select_svd_k(Y), 4)
})
