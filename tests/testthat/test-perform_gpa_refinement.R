library(testthat)
skip_on_cran()
skip_if_not_installed("expm")

random_SOk <- function(k) {
  M <- matrix(rnorm(k * k), k, k)
  qr_M <- qr(M)
  Q <- qr.Q(qr_M)
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  Q
}

set.seed(123)

## Basic convergence check (TCK-GPA-004)
test_that("perform_gpa_refinement aligns subjects to common template", {
  m <- 4; k <- 3
  A_base <- matrix(rnorm(m * k), nrow = m, ncol = k)
  R2 <- random_SOk(k)
  A_list <- list(A_base, A_base %*% R2)

  res <- hatsa:::perform_gpa_refinement(
    A_originals_list = A_list,
    n_refine = 5,
    k = k,
    m_parcel_rows = m
  )

  aligned1 <- A_list[[1]] %*% res$R_final_list[[1]]
  aligned2 <- A_list[[2]] %*% res$R_final_list[[2]]
  expect_equal(aligned1, res$T_anchor_final, tolerance = 1e-6)
  expect_equal(aligned2, res$T_anchor_final, tolerance = 1e-6)
})

## n_refine = 0 should compute rotations against mean template (TCK-GPA-005)
test_that("perform_gpa_refinement handles n_refine = 0 correctly", {
  m <- 4; k <- 3
  A_base <- matrix(rnorm(m * k), nrow = m, ncol = k)
  R2 <- random_SOk(k)
  A_list <- list(A_base, A_base %*% R2)
  mean_A <- Reduce(`+`, A_list) / length(A_list)

  res0 <- hatsa:::perform_gpa_refinement(
    A_originals_list = A_list,
    n_refine = 0,
    k = k,
    m_parcel_rows = m
  )

  expect_equal(res0$T_anchor_final, mean_A, tolerance = 1e-8)
  expect_equal(res0$R_final_list[[1]], solve_procrustes_rotation(A_list[[1]], mean_A), tolerance = 1e-6)
  expect_equal(res0$R_final_list[[2]], solve_procrustes_rotation(A_list[[2]], mean_A), tolerance = 1e-6)
})
