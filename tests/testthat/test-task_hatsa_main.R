# test-task_hatsa_main.R

library(testthat)
skip_on_cran()

describe("run_task_hatsa end-to-end pipeline", {

make_toy_subject_data <- function(N, T, V) {
  lapply(1:N, function(i) matrix(rnorm(T*V), nrow=T, ncol=V))
}
make_toy_task_data <- function(N, C, V) {
  lapply(1:N, function(i) matrix(rnorm(C*V), nrow=C, ncol=V))
}

N <- 3; T <- 10; V <- 5; C <- 4
subject_data_list <- make_toy_subject_data(N, T, V)
task_data_list <- make_toy_task_data(N, C, V)
anchor_indices <- 1:2
spectral_rank_k <- 2
parcel_names <- paste0("P", 1:V)

# Minimal working example for core_hatsa

test_that("run_task_hatsa works for core_hatsa", {
  res <- run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    parcel_names = parcel_names,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
  expect_equal(length(res$U_aligned_list), N)
  expect_true(all(sapply(res$U_aligned_list, function(x) is.null(x) || (is.matrix(x) && ncol(x) == spectral_rank_k && nrow(x) == V))))
})

# Minimal working example for lambda_blend

test_that("run_task_hatsa works for lambda_blend", {
  res <- run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = task_data_list,
    task_method = "lambda_blend",
    lambda_blend_value = 0.2,
    parcel_names = parcel_names,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
  expect_equal(length(res$U_aligned_list), N)
})

# Minimal working example for gev_patch (if supported)
test_that("run_task_hatsa works for gev_patch", {
  res <- run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = task_data_list,
    task_method = "gev_patch",
    k_gev_dims = 2,
    parcel_names = parcel_names,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
  expect_equal(length(res$U_aligned_list), N)
})

# Handles missing/NULL task_data_list for core_hatsa
test_that("run_task_hatsa allows NULL task_data_list for core_hatsa", {
  expect_warning(run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    parcel_names = parcel_names,
    verbose = FALSE
  ), "backward compatibility")
})

# Handles errors for missing required arguments
test_that("run_task_hatsa errors for missing required arguments", {
  expect_error(run_task_hatsa(
    subject_data_list = NULL,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    parcel_names = parcel_names,
    verbose = FALSE
  ))
  expect_error(run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = NULL,
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    parcel_names = parcel_names,
    verbose = FALSE
  ))
})

# Handles row_augmentation and reliability_scores_list
test_that("run_task_hatsa works with row_augmentation and reliability_scores_list", {
  reliability_scores_list <- lapply(1:N, function(i) runif(C))
  res <- run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = task_data_list,
    task_method = "lambda_blend",
    row_augmentation = TRUE,
    omega_mode = "adaptive",
    reliability_scores_list = reliability_scores_list,
    parcel_names = parcel_names,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
})

# Handles custom W_task_helper_func
test_that("run_task_hatsa works with custom W_task_helper_func", {
  custom_fun <- function(mat, parcel_names, k_conn_task_pos, k_conn_task_neg, ...) {
    Matrix::Diagonal(n = ncol(mat), x = 1)
  }
  res <- run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = task_data_list,
    task_method = "lambda_blend",
    W_task_helper_func = custom_fun,
    parcel_names = parcel_names,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
})

# Edge cases: empty subject_data_list, anchor_indices out of bounds
test_that("run_task_hatsa errors for empty subject_data_list or out-of-bounds anchors", {
  expect_error(run_task_hatsa(
    subject_data_list = list(),
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    parcel_names = parcel_names,
    verbose = FALSE
  ))
  expect_error(run_task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = c(100),
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    parcel_names = parcel_names,
    verbose = FALSE
  ))
}) })
