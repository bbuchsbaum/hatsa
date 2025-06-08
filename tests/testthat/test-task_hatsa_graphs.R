# test-task_hatsa_graphs.R

library(testthat)
library(Matrix)
skip_on_cran()

source("../../R/task_hatsa.R")

describe("compute_W_task_from_activations and related graph construction", {

test_that("compute_W_task_from_activations returns correct type and dimensions", {
  set.seed(1)
  C <- 5; Vp <- 4
  mat <- matrix(rnorm(C*Vp), nrow=C, ncol=Vp)
  pnames <- paste0("P", 1:Vp)
  W <- compute_W_task_from_activations(mat, pnames, 2, 2)
  expect_s4_class(W, "dgCMatrix")
  expect_equal(dim(W), c(Vp, Vp))
  expect_true(Matrix::isSymmetric(W))
})

test_that("compute_W_task_from_activations handles zero-variance columns", {
  mat <- matrix(1, nrow=5, ncol=4)
  pnames <- paste0("P", 1:4)
  expect_warning(W <- compute_W_task_from_activations(mat, pnames, 2, 2), NA)
  expect_s4_class(W, "dgCMatrix")
})

test_that("compute_W_task_from_activations handles C=1 and C=0", {
  mat1 <- matrix(rnorm(4), nrow=1)
  mat0 <- matrix(numeric(0), nrow=0, ncol=4)
  pnames <- paste0("P", 1:4)
  W1 <- compute_W_task_from_activations(mat1, pnames, 2, 2)
  W0 <- compute_W_task_from_activations(mat0, pnames, 2, 2)
  expect_equal(dim(W1), c(4,4))
  expect_equal(dim(W0), c(4,4))
  expect_true(all(W1 == 0))
  expect_true(all(W0 == 0))
})

test_that("compute_W_task_from_activations output is symmetric and sparse", {
  mat <- matrix(rnorm(20), nrow=5)
  pnames <- paste0("P", 1:4)
  W <- compute_W_task_from_activations(mat[,1:4], pnames, 2, 2)
  expect_true(Matrix::isSymmetric(W))
  expect_true(inherits(W, "sparseMatrix"))
})

test_that("compute_W_task_from_activations z-scores nonzero edges", {
  mat <- matrix(rnorm(20), nrow=5)
  pnames <- paste0("P", 1:4)
  W <- compute_W_task_from_activations(mat[,1:4], pnames, 2, 2)
  if (length(W@x) > 0) {
    expect_equal(mean(W@x), 0, tolerance=1e-8)
    expect_equal(sd(W@x), 1, tolerance=0.1)
  }
})

test_that("compute_W_task_from_activations supports custom similarity_method", {
  mat <- matrix(rnorm(20), nrow=5)
  pnames <- paste0("P", 1:4)
  simfun <- function(x) matrix(1, ncol(x), ncol(x))
  W <- compute_W_task_from_activations(mat[,1:4], pnames, 2, 2, similarity_method=simfun)
  expect_s4_class(W, "dgCMatrix")
})

test_that("compute_W_task_from_activations handles all-zero and all-NA input", {
  mat0 <- matrix(0, nrow=5, ncol=4)
  matNA <- matrix(NA_real_, nrow=5, ncol=4)
  pnames <- paste0("P", 1:4)
  W0 <- compute_W_task_from_activations(mat0, pnames, 2, 2)
  WNA <- compute_W_task_from_activations(matNA, pnames, 2, 2)
  expect_true(all(W0 == 0))
  expect_true(all(WNA == 0))
})

test_that("compute_W_task_from_activations errors on mismatched parcel_names", {
  mat <- matrix(rnorm(20), nrow=5)
  expect_error(compute_W_task_from_activations(mat[,1:4], paste0("P", 1:3), 2, 2))
})

# If compute_W_task_from_encoding is present, add a basic test
if (exists("compute_W_task_from_encoding", mode="function")) {
  test_that("compute_W_task_from_encoding returns correct type and dimensions", {
    mat <- matrix(rnorm(20), nrow=4, ncol=5)
    pnames <- paste0("P", 1:4)
    W <- compute_W_task_from_encoding(mat, pnames, 2, 2)
    expect_s4_class(W, "dgCMatrix")
    expect_equal(dim(W), c(4,4))
    expect_true(Matrix::isSymmetric(W))
  })
} })
