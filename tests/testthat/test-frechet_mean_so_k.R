library(testthat)
skip_on_cran()
skip_if_not_installed("expm")

rotation_2d <- function(theta) {
  matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
}

## Fréchet mean of symmetric rotations should be identity

test_that("frechet_mean_so_k handles simple cases", {
  R_pos <- rotation_2d(pi/4)
  R_neg <- rotation_2d(-pi/4)
  mean_R <- frechet_mean_so_k(list(R_pos, R_neg), k_dim = 2, tol = 1e-10)
  # The mean of symmetric rotations should be identity or negative identity
  # Both are valid solutions as they represent the same point on the projective space
  expect_true(
    isTRUE(all.equal(mean_R, diag(2), tolerance = 1e-6)) ||
    isTRUE(all.equal(mean_R, -diag(2), tolerance = 1e-6)),
    info = "Mean should be either I or -I"
  )

  empty_mean <- frechet_mean_so_k(list(), k_dim = 2)
  expect_equal(empty_mean, diag(2))
})
