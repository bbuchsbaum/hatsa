library(testthat)

# Tests for SPD matrix operations in riemannian_geometry.R

# Helper function to generate random SPD matrix
.make_spd <- function(p) {
  A <- matrix(rnorm(p * p), p, p)
  crossprod(A) + diag(p) * 0.5
}


test_that("matrix_logm_spd, matrix_expm_spd and matrix_sqrt_spd are consistent", {
  set.seed(1)
  S <- .make_spd(3)
  logS <- matrix_logm_spd(S)
  S_recon <- matrix_expm_spd(logS)
  expect_equal(S_recon, S, tolerance = 1e-6)

  S_sqrt <- matrix_sqrt_spd(S)
  expect_equal(S_sqrt %*% S_sqrt, S, tolerance = 1e-6)

  # zero-dimension edge case
  S0 <- matrix(0,0,0)
  expect_equal(matrix_logm_spd(S0), S0)
  expect_equal(matrix_expm_spd(S0), S0)
  expect_equal(matrix_sqrt_spd(S0), S0)
})


test_that("airm_distance matches log-Euclidean distance for commuting matrices", {
  S1 <- diag(c(2,3,4))
  S2 <- diag(c(1,5,2))
  dist_airm <- airm_distance(S1, S2)
  dist_loge <- riemannian_distance_spd(S1, S2)
  expect_equal(dist_airm, dist_loge, tolerance = 1e-8)
  expect_equal(airm_distance(S1, S1), 0)
})


test_that("logmap and expmap functions are inverses at base point", {
  set.seed(2)
  S1 <- .make_spd(3)
  S2 <- .make_spd(3)

  v_log <- logmap_spd_logeuclidean(S1, S2)
  S2_rec <- expmap_spd_logeuclidean(S1, v_log)
  expect_equal(S2_rec, S2, tolerance = 1e-6)

  v_airm <- logmap_spd_airm(S1, S2)
  S2_rec_airm <- expmap_spd_airm(S1, v_airm)
  expect_equal(S2_rec_airm, S2, tolerance = 1e-6)
})

