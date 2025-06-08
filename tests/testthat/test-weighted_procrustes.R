library(testthat)
skip_on_cran()

random_SOk <- function(k) {
  M <- matrix(rnorm(k * k), k, k)
  qr_M <- qr(M)
  Q <- qr.Q(qr_M)
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  Q
}

set.seed(123)

test_that("weighted solver reduces residuals for targeted rows", {
  m_parcel <- 3
  m_task <- 2
  k <- 3
  N <- m_parcel + m_task

  T_target <- matrix(rnorm(N * k), nrow = N, ncol = k)
  R_parcel <- random_SOk(k)
  R_task <- random_SOk(k)

  A_source <- rbind(
    T_target[1:m_parcel, ] %*% t(R_parcel),
    T_target[(m_parcel + 1):N, ] %*% t(R_task)
  )

  R_unw <- solve_procrustes_rotation(A_source, T_target)
  weights <- list(parcel = 2, condition = 0.5)
  R_w <- solve_procrustes_rotation_weighted(
    A_source, T_target,
    m_parcel_rows = m_parcel,
    m_task_rows = m_task,
    fixed_omega_weights = weights
  )

  resid_unw <- A_source %*% R_unw - T_target
  resid_w <- A_source %*% R_w - T_target

  sse_parcel_unw <- sum(resid_unw[1:m_parcel, ]^2)
  sse_parcel_w <- sum(resid_w[1:m_parcel, ]^2)

  expect_lt(sse_parcel_w, sse_parcel_unw)
})

test_that("equal weights match unweighted solution", {
  m_parcel <- 4
  m_task <- 3
  k <- 4
  N <- m_parcel + m_task

  T_target <- matrix(rnorm(N * k), nrow = N, ncol = k)
  R1 <- random_SOk(k)
  R2 <- random_SOk(k)

  A_source <- rbind(
    T_target[1:m_parcel, ] %*% t(R1),
    T_target[(m_parcel + 1):N, ] %*% t(R2)
  )

  R_unw <- solve_procrustes_rotation(A_source, T_target)
  R_eq <- solve_procrustes_rotation_weighted(
    A_source, T_target,
    m_parcel_rows = m_parcel,
    m_task_rows = m_task,
    fixed_omega_weights = list(parcel = 1, condition = 1)
  )

  expect_equal(R_eq, R_unw, tolerance = 1e-6)
})

