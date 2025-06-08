library(testthat)
library(Matrix)

describe("compute_graph_correlation", {

  test_that("handles empty graphs", {
    V <- 5
    W1 <- Matrix::Matrix(0, nrow = V, ncol = V, sparse = TRUE)
    W2 <- Matrix::Matrix(0, nrow = V, ncol = V, sparse = TRUE)
    expect_true(is.na(compute_graph_correlation(W1, W2)))
  })

  test_that("sampling uses set.seed reproducibly", {
    set.seed(123)
    V <- 20
    mat1 <- matrix(rnorm(V * V), V, V)
    mat1 <- (mat1 + t(mat1)) / 2
    diag(mat1) <- 0
    W1 <- Matrix::Matrix(mat1, sparse = TRUE)
    W1@x <- scale(W1@x)[,1]

    mat2 <- matrix(rnorm(V * V), V, V)
    mat2 <- (mat2 + t(mat2)) / 2
    diag(mat2) <- 0
    W2 <- Matrix::Matrix(mat2, sparse = TRUE)
    W2@x <- scale(W2@x)[,1]

    set.seed(1)
    rho1 <- compute_graph_correlation(W1, W2, max_edges = 10)
    set.seed(1)
    rho2 <- compute_graph_correlation(W1, W2, max_edges = 10)
    expect_equal(rho1, rho2)
    expect_true(is.numeric(rho1))
  })

  test_that("returns NA when edge variance is zero", {
    V <- 4
    mat1 <- matrix(0, V, V)
    mat1[upper.tri(mat1)] <- 1
    mat1 <- mat1 + t(mat1)
    diag(mat1) <- 0
    W1 <- Matrix::Matrix(mat1, sparse = TRUE)
    W1@x <- rep(1, length(W1@x))
    W2 <- Matrix::Matrix(0, nrow = V, ncol = V, sparse = TRUE)

    expect_warning(res <- compute_graph_correlation(W1, W2), "Zero variance")
    expect_true(is.na(res))
  })

})
