library(testthat)
skip_on_cran()
skip_if_not_installed("expm")

random_SOk <- function(k) {
  M <- matrix(rnorm(k * k), k, k)
  qr_M <- qr(M)
  Q <- qr.Q(qr_M)
  if (det(Q) < 0) Q[,1] <- -Q[,1]
  Q
}

set.seed(42)


test_that("solve_procrustes_rotation recovers known rotation", {
  m <- 5; k <- 3
  A_anchor <- matrix(rnorm(m * k), nrow = m, ncol = k)
  R_true <- random_SOk(k)
  T_anchor <- A_anchor %*% R_true
  R_est <- solve_procrustes_rotation(A_anchor, T_anchor)
  expect_equal(R_est, R_true, tolerance = 1e-6)
})


test_that("perform_geometric_gpa_refinement handles zero subjects and invalid matrices", {
  expect_warning(res_empty <- hatsa:::perform_geometric_gpa_refinement(
    A_originals_list = list(),
    n_refine = 3, k = 3, m_rows = 4, verbose = FALSE
  ), "empty")
  expect_length(res_empty$R_final_list, 0)
  expect_true(all(is.na(res_empty$T_anchor_final)))

  bad_list <- list(matrix(1:6, 2, 3))
  expect_warning(
    result <- hatsa:::perform_geometric_gpa_refinement(bad_list, n_refine = 2, k = 3, m_rows = 4, verbose = FALSE),
    "invalid"
  )
  expect_equal(length(result$R_final_list), 1)
  expect_true(all(result$R_final_list[[1]] == diag(3)))
})


test_that("perform_geometric_gpa_refinement converges for SVD and Riemannian modes", {
  m <- 5; k <- 3
  A_base <- matrix(rnorm(m * k), m, k)
  R2 <- random_SOk(k)
  A_list <- list(A_base, A_base %*% R2)

  res_svd <- hatsa:::perform_geometric_gpa_refinement(
    A_originals_list = A_list,
    n_refine = 6,
    k = k,
    m_rows = m,
    rotation_mode = "svd",
    verbose = FALSE
  )
  aligned1 <- A_list[[1]] %*% res_svd$R_final_list[[1]]
  aligned2 <- A_list[[2]] %*% res_svd$R_final_list[[2]]
  expect_equal(aligned1, res_svd$T_anchor_final, tolerance = 1e-5)
  expect_equal(aligned2, res_svd$T_anchor_final, tolerance = 1e-5)

  res_riem <- hatsa:::perform_geometric_gpa_refinement(
    A_originals_list = A_list,
    n_refine = 6,
    k = k,
    m_rows = m,
    rotation_mode = "riemannian",
    verbose = FALSE
  )
  aligned1_r <- A_list[[1]] %*% res_riem$R_final_list[[1]]
  aligned2_r <- A_list[[2]] %*% res_riem$R_final_list[[2]]
  expect_equal(aligned1_r, res_riem$T_anchor_final, tolerance = 1e-5)
  expect_equal(aligned2_r, res_riem$T_anchor_final, tolerance = 1e-5)
})
