test_that("getChat inverts simple prior matrix correctly", {
    prior <- diag(1, 3)
    expect_message(
        Chat <- getChat(prior, scale = FALSE),
        "Inverting..."
    )
    expect_equal(dim(Chat), c(3, 3))
    expect_equal(as.numeric(diag(Chat)), rep(1 / 26, 3))
})
