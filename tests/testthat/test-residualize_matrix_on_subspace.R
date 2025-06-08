library(testthat)

describe("residualize_matrix_on_subspace helper", {

  test_that("residualization orthogonality for row orientation", {
    set.seed(42)
    C <- 5; k <- 3; m <- 2
    Y <- matrix(rnorm(C * k), nrow = C, ncol = k)
    X <- matrix(rnorm(m * k), nrow = m, ncol = k)
    Y_res <- residualize_matrix_on_subspace(Y, X)
    expect_equal(dim(Y_res), c(C, k))
    expect_true(max(abs(Y_res %*% t(X))) < 1e-6)
  })

  test_that("orientation switch when k > C gives same result", {
    set.seed(99)
    C <- 3; k <- 5; m <- 2
    Y <- matrix(rnorm(C * k), nrow = C, ncol = k)
    X <- matrix(rnorm(m * k), nrow = m, ncol = k)
    Y_res <- residualize_matrix_on_subspace(Y, X)

    qr_X_eff <- qr(t(X))
    ref <- t(qr.resid(qr_X_eff, t(Y)))

    expect_equal(Y_res, ref, tolerance = 1e-8)
    expect_true(max(abs(Y_res %*% t(X))) < 1e-6)
  })

})
