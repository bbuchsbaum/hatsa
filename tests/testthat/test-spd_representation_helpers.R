library(testthat)
library(Matrix)

# Helper to create a minimal hatsa_projector-like object for tests
create_basic_projector <- function(U_list, stacked, V, k) {
  N <- length(U_list)
  obj <- list(
    U_aligned_list = U_list,
    s = stacked,
    block_indices = split(seq_len(N * V), rep(seq_len(N), each = V)),
    parameters = list(N_subjects = N, V_p = V, k = k)
  )
  class(obj) <- "hatsa_projector"
  obj
}

test_that(".regularize_spd enforces positive definiteness", {
  S <- matrix(c(0, 0.5, 0.5, 0), 2, 2)
  res <- .regularize_spd(S, epsilon_spd = 0.1)
  expect_true(isSymmetric(res))
  eig <- eigen(res, symmetric = TRUE)$values
  expect_true(all(eig >= 0.1 - 1e-8))

  expect_error(.regularize_spd(matrix(c(1,2,3,4),2,2)), "numeric, symmetric")
  expect_equal(.regularize_spd(matrix(0,0,0)), matrix(0,0,0))
})


test_that(".densify_symmetrize_regularize handles sparse and large matrices", {
  M <- Matrix::sparseMatrix(i=c(1,2), j=c(2,1), x=c(1,3), dims=c(2,2))
  expect_warning(res <- .densify_symmetrize_regularize(M, 0.05), "not symmetric")
  expect_true(is.matrix(res) && isSymmetric(res))
  eig <- eigen(res, symmetric = TRUE)$values
  expect_true(all(eig > 0))

  M_large <- Matrix::sparseMatrix(i=1:3, j=1:3, x=1:3, dims=c(1001,1001))
  expect_warning(res_large <- .densify_symmetrize_regularize(M_large, 1e-4), "Coercing large sparse matrix")
  expect_equal(dim(res_large), c(1001,1001))
  eig_large <- eigen(res_large, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_large > 0))
})


test_that("get_spd_representations retrieves covariance-based SPD matrices", {
  set.seed(1)
  V <- 4; k <- 2
  U1 <- matrix(rnorm(V*k), V, k)
  U2 <- matrix(rnorm(V*k), V, k)
  stacked <- rbind(U1, U2)
  proj <- create_basic_projector(list(U1, NULL), stacked, V, k)

  res1 <- .get_subject_aligned_sketch(proj, 1)
  expect_equal(res1, U1)
  res2 <- .get_subject_aligned_sketch(proj, 2)
  expect_equal(res2, U2)
  expect_error(.get_subject_aligned_sketch(proj, 3), "Invalid subject_idx")

  all_spd <- get_spd_representations(proj, type="cov_coeffs")
  expect_length(all_spd, 2)
  expect_equal(names(all_spd), c("1","2"))
  expect_equal(all_spd[[1]], .regularize_spd(cov(U1)))
  expect_equal(all_spd[[2]], .regularize_spd(cov(U2)))

  spd_first <- get_spd_representations(proj, type="cov_coeffs", subject_idx = 1)
  expect_true(is.matrix(spd_first))
  expect_equal(spd_first, .regularize_spd(cov(U1)))
  expect_error(get_spd_representations(proj, subject_idx = 3), "Invalid subject_idx")
})

