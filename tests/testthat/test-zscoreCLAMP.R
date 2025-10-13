test_that("zscoreCLAMP centers and scales each row", {
    Y <- matrix(c(2, 4, 6, 8),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("g1", "g2"), c("s1", "s2"))
    )

    mu <- rowMeans(Y)
    var_pop <- rowMeans((Y - mu)^2)

    rowStats <- data.frame(
        mean = mu,
        variance = var_pop,
        row.names = rownames(Y),
        check.names = FALSE
    )

    Z <- zscoreCLAMP(Y, rowStats = rowStats)

    expected <- sweep(Y, 1, mu, "-")
    expected <- sweep(expected, 1, sqrt(var_pop), "/")

    expect_equal(unname(as.matrix(Z)), unname(as.matrix(expected)))
})
