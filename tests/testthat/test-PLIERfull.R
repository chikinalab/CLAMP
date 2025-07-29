context("PLIERfull")

test_that("PLIERfull returns list with B, Z, U when doCrossval=FALSE", {
    mat <- matrix(rnorm(100), 10, 10)
    svdres <- rsvd(mat, k = 5)
    base <- PLIERbase(Y = mat, k = 5, svdres = svdres, trace = FALSE)
    priorMat <- matrix(1, nrow(mat), 5)
    full <- PLIERfull(
        Y = mat, priorMat = priorMat, svdres = svdres,
        plier.base.result = base, k = 5,
        doCrossval = FALSE, trace = FALSE, max.U.updates = 0
    )
    expect_type(full, "list")
    expect_true(all(c("B", "Z", "U") %in% names(full)))
})
