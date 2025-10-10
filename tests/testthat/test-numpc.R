test_that("num.pc returns integer >= 1 for random data", {
    mat <- matrix(rnorm(50), 5, 10)
    k <- num.pc(list(d = svd(mat)$d))
    expect_true(is.numeric(k) && k >= 1)
})

test_that("num.pc errors on missing d element", {
    expect_error(num.pc(data = list()))
})
