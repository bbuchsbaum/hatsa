describe("compute_gev_spectrum_diagnostics", {
  test_that("compute_gev_spectrum_diagnostics computes expected statistics", {
  lambda <- c(-1.2, -0.5, 0.3, 0.8, 1.5)
  thresh <- 0.8
  stats <- compute_gev_spectrum_diagnostics(lambda, thresh)
  expect_equal(stats$n_eigenvalues, length(lambda))
  expect_equal(stats$min_eigenvalue, min(lambda))
  expect_equal(stats$max_eigenvalue, max(lambda))
  expect_equal(stats$n_below_thresh, sum(abs(lambda) < thresh))
  expect_equal(stats$n_above_thresh, sum(abs(lambda) >= thresh))
  expect_equal(stats$prop_below_thresh, stats$n_below_thresh / length(lambda))
  expect_equal(stats$prop_above_thresh, stats$n_above_thresh / length(lambda))
})

test_that("compute_gev_spectrum_diagnostics handles empty and invalid input", {
  empty_stats <- compute_gev_spectrum_diagnostics(numeric(0), 0.5)
  expect_equal(empty_stats$n_eigenvalues, 0)
  expect_true(all(is.na(c(empty_stats$min_eigenvalue, empty_stats$max_eigenvalue,
                          empty_stats$mean_eigenvalue, empty_stats$median_eigenvalue,
                          empty_stats$sd_eigenvalue, empty_stats$prop_below_thresh,
                          empty_stats$prop_above_thresh))))

  expect_error(compute_gev_spectrum_diagnostics("bad", 0.5), "numeric vector")
  expect_error(compute_gev_spectrum_diagnostics(c(1,2), c(0.5,0.2)), "single numeric")
})
})
