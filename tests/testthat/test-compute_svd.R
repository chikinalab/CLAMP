test_that("compute_svd returns d, u, v for a dense matrix", {
  set.seed(1)
  Y <- matrix(rnorm(100), nrow = 20, ncol = 5)
  res <- compute_svd(Y, k = 3)
  expect_true(all(c("d", "u", "v") %in% names(res)))
  expect_length(res$d, 3)
  expect_equal(dim(res$u), c(20, 3))
  expect_equal(dim(res$v), c(5, 3))
})

test_that("compute_svd works on sparse dgCMatrix", {
  set.seed(1)
  Y <- Matrix::rsparsematrix(20, 10, density = 0.3)
  res <- compute_svd(Y, k = 3)
  expect_length(res$d, 3)
  expect_equal(dim(res$u), c(20, 3))
  expect_equal(dim(res$v), c(10, 3))
})

test_that("compute_svd works on FBM", {
  skip_if_not_installed("bigstatsr")
  set.seed(1)
  Y <- bigstatsr::FBM(nrow = 30, ncol = 10, init = rnorm(300))
  res <- compute_svd(Y, k = 3)
  expect_length(res$d, 3)
})

test_that("compute_svd defaults k to select_svd_k(Y) when k is NULL", {
  set.seed(1)
  Y <- matrix(rnorm(200), nrow = 20, ncol = 10)
  res <- compute_svd(Y)
  expect_length(res$d, select_svd_k(Y))
})

test_that("compute_svd recovers singular values on a known rank-k matrix", {
  set.seed(1)
  u <- qr.Q(qr(matrix(rnorm(50), 10, 5)))
  v <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
  d_true <- c(10, 5, 3, 1, 0.5)
  Y <- u %*% diag(d_true) %*% t(v)
  res <- compute_svd(Y, k = 5)
  expect_equal(sort(res$d, decreasing = TRUE), d_true, tolerance = 1e-6)
})
