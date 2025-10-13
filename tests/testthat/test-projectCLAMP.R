test_that("projectCLAMP projects new data correctly", {
    mat <- matrix(rnorm(100), 10, 10)
    svdres <- rsvd(mat, k = 5)
    base <- CLAMPbase(Y = mat, k = 5, svdres = svdres, trace = FALSE)
    newdata <- matrix(rnorm(50), 10, 5)
    Bnew <- projectCLAMP(base, newdata)
    expect_true(is.matrix(Bnew))
    expect_equal(dim(Bnew), c(5, 5))
})
