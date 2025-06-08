library(testthat)


test_that("assume_orthonormal check warns or errors appropriately", {
  # Setup simple feature matrix and basis
  V_p <- 5
  C <- 3
  features <- matrix(rnorm(V_p * C), nrow = V_p)

  # Case 1: k_dims_basis > V_p_basis should error
  k_big <- V_p + 1
  U_big <- matrix(rnorm(V_p * k_big), nrow = V_p)
  expect_error(
    project_features_to_spectral_space(features, U_big, assume_orthonormal = TRUE)
  )

  # Case 2: basis not orthonormal triggers warning and falls back
  k <- 3
  U_non_ortho <- matrix(rnorm(V_p * k), nrow = V_p)
  expect_warning(
    res_warn <- project_features_to_spectral_space(features, U_non_ortho, assume_orthonormal = TRUE),
    "not orthonormal"
  )
  res_regular <- project_features_to_spectral_space(features, U_non_ortho, assume_orthonormal = FALSE)
  expect_equal(res_warn, res_regular)
})
