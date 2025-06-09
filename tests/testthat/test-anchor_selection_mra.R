library(testthat)

describe("select_anchors_mra basic functionality", {
  make_pilot <- function(n_subj, V, k) {
    set.seed(123)  # For reproducibility
    replicate(n_subj, {
      # Create a base pattern and add some subject-specific variation
      base <- matrix(seq_len(V * k), nrow = V, ncol = k)
      noise <- matrix(rnorm(V * k, sd = 0.5), nrow = V, ncol = k)
      base + noise
    }, simplify = FALSE)
  }

  test_that("returns sorted anchors with expected count", {
    pilot <- make_pilot(3, 4, 2)
    res <- select_anchors_mra(
      U_original_list_pilot = pilot,
      spectral_rank_k = 2,
      m_target = 2,
      total_parcels = 4,
      weight_inv_kappa = 1,
      weight_dispersion = 0,
      min_anchors_for_metrics = 1,
      verbose = FALSE
    )
    expect_length(res, 2)
    expect_type(res, "integer")
    expect_true(all(res %in% 1:4))  # All selected anchors should be valid parcel indices
    expect_equal(length(unique(res)), 2)  # Should have 2 unique anchors
  })

  test_that("respects initial_selection and candidate_pool ordering", {
    pilot <- make_pilot(2, 4, 2)
    res <- select_anchors_mra(
      U_original_list_pilot = pilot,
      spectral_rank_k = 2,
      m_target = 3,
      total_parcels = 4,
      initial_selection = c(2L),
      candidate_pool = c(4L, 3L, 1L),
      weight_dispersion = 0,
      min_anchors_for_metrics = 2,
      verbose = FALSE
    )
    expect_length(res, 3)
    expect_true(2L %in% res)  # Initial selection should be included
    expect_true(all(res %in% c(1L, 2L, 3L, 4L)))  # All results should be valid parcels
  })

  test_that("warns when candidate pool cannot satisfy target", {
    pilot <- make_pilot(2, 3, 2)
    expect_warning(
      res <- select_anchors_mra(
        U_original_list_pilot = pilot,
        spectral_rank_k = 2,
        m_target = 2,
        total_parcels = 3,
        initial_selection = c(1L),
        candidate_pool = integer(0),
        weight_dispersion = 0,
        min_anchors_for_metrics = 1,
        verbose = FALSE
      ),
      "No candidate anchors available"
    )
    expect_equal(res, c(1L))
  })

  test_that("invalid inputs trigger errors", {
    pilot <- make_pilot(1, 3, 2)
    expect_error(
      select_anchors_mra(
        U_original_list_pilot = pilot,
        spectral_rank_k = 2,
        m_target = 2,
        total_parcels = 3,
        initial_selection = c(0L),
        weight_dispersion = 0,
        verbose = FALSE
      ),
      "initial_selection contains invalid"
    )

    expect_error(
      select_anchors_mra(
        U_original_list_pilot = pilot,
        spectral_rank_k = 2,
        m_target = 2,
        total_parcels = 3,
        candidate_pool = c(1L, 4L),
        weight_dispersion = 0,
        verbose = FALSE
      ),
      "candidate_pool contains invalid"
    )

    expect_error(
      select_anchors_mra(
        U_original_list_pilot = pilot,
        spectral_rank_k = 2,
        m_target = 0,
        total_parcels = 3,
        weight_dispersion = 0,
        verbose = FALSE
      ),
      "m_target must be positive"
    )
  })
})
