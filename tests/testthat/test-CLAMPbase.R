test_that("CLAMPbase returns list with B and Z of correct dimensions", {
    mat <- matrix(rnorm(100), 10, 10)
    svdres <- rsvd(mat, k = 5)
    base <- CLAMPbase(Y = mat, k = 5, svdres = svdres, trace = FALSE)
    expect_type(base, "list")
    expect_true(all(c("B", "Z") %in% names(base)))
    expect_equal(dim(base$B), c(5, 10))
    expect_equal(dim(base$Z), c(10, 5))
})
