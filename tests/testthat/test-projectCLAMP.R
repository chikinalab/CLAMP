test_that("projectCLAMP projects new data correctly", {
    mat <- matrix(rnorm(100), 10, 10)
    svdres <- rsvd(mat, k = 5)
    base <- CLAMPbase(Y = mat, clamp_k = 5, svdres = svdres, trace = FALSE)
    newdata <- matrix(rnorm(50), 10, 5)
    Bnew <- projectCLAMP(base, newdata, verbose = FALSE)
    expect_true(is.matrix(Bnew))
    expect_equal(dim(Bnew), c(5, 5))
})

test_that("projectCLAMP aligns common row names in model order", {
    Z <- matrix(
        c(1, 0, 2, 0, 1, 1),
        nrow = 3,
        dimnames = list(c("GeneA", "GeneB", "GeneC"), c("LV1", "LV2"))
    )
    newdata <- matrix(
        c(10, 20, 30, 40, 50, 60, 70, 80),
        nrow = 4,
        dimnames = list(c("GeneC", "GeneX", "GeneA", "GeneB"), c("S1", "S2"))
    )
    res <- list(Z = Z, L2 = 0.1)

    Bnew <- projectCLAMP(res, newdata, verbose = FALSE)
    expected <- solve(t(Z) %*% Z + 0.1 * diag(ncol(Z))) %*%
        t(Z) %*% newdata[rownames(Z), , drop = FALSE]

    expect_equal(Bnew, expected)
    expect_equal(rownames(Bnew), colnames(Z))
    expect_equal(colnames(Bnew), colnames(newdata))
})

test_that("projectCLAMP can require pre-aligned row order", {
    Z <- matrix(
        rnorm(12),
        nrow = 4,
        dimnames = list(paste0("Gene", 1:4), paste0("LV", 1:3))
    )
    newdata <- matrix(
        rnorm(8),
        nrow = 4,
        dimnames = list(rev(rownames(Z)), paste0("S", 1:2))
    )
    res <- list(Z = Z, L2 = 0.1)

    expect_error(
        projectCLAMP(res, newdata, align = FALSE, verbose = FALSE),
        "row names are not in the same order"
    )
})
