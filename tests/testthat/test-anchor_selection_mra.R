library(testthat)

describe("select_anchors_mra basic functionality", {
  make_pilot <- function(n_subj, V, k) {
    replicate(n_subj, matrix(seq_len(V * k), nrow = V, ncol = k), simplify = FALSE)
  }

  test_that("returns sorted anchors with expected count", {
    pilot <- make_pilot(3, 4, 2)
    res <- select_anchors_mra(
      U_original_list_pilot = pilot,
      k_spectral_rank = 2,
      m_target = 2,
      total_parcels = 4,
      weight_dispersion = 0,
      min_anchors_for_metrics = 3,
      verbose = FALSE
    )
    expect_length(res, 2)
    expect_type(res, "integer")
    expect_equal(res, c(1L, 2L))
  })

  test_that("respects initial_selection and candidate_pool ordering", {
    pilot <- make_pilot(2, 4, 2)
    res <- select_anchors_mra(
      U_original_list_pilot = pilot,
      k_spectral_rank = 2,
      m_target = 3,
      total_parcels = 4,
      initial_selection = c(2L),
      candidate_pool = c(4L, 3L, 1L),
      weight_dispersion = 0,
      min_anchors_for_metrics = 4,
      verbose = FALSE
    )
    expect_length(res, 3)
    expect_equal(res, c(2L, 3L, 4L))
  })

  test_that("warns when candidate pool cannot satisfy target", {
    pilot <- make_pilot(2, 3, 2)
    expect_warning(
      res <- select_anchors_mra(
        U_original_list_pilot = pilot,
        k_spectral_rank = 2,
        m_target = 2,
        total_parcels = 3,
        initial_selection = c(1L),
        candidate_pool = integer(0),
        weight_dispersion = 0,
        min_anchors_for_metrics = 3,
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
        k_spectral_rank = 2,
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
        k_spectral_rank = 2,
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
        k_spectral_rank = 2,
        m_target = 0,
        total_parcels = 3,
        weight_dispersion = 0,
        verbose = FALSE
      ),
      "m_target must be positive"
    )
  })
})
