# test-task_hatsa.R

library(testthat)
skip_on_cran()

# Helper functions to generate test data
make_toy_subject_data <- function(N, T, V) {
  lapply(1:N, function(i) matrix(rnorm(T*V), nrow=T, ncol=V))
}
make_toy_task_data <- function(N, C, V) {
  lapply(1:N, function(i) matrix(rnorm(C*V), nrow=C, ncol=V))
}

# Test data
N <- 3; T <- 10; V <- 5; C <- 4
subject_data_list <- make_toy_subject_data(N, T, V)
task_data_list <- make_toy_task_data(N, C, V)
anchor_indices <- 1:2
spectral_rank_k <- 2
parcel_names <- paste0("P", 1:V)

test_that("task_hatsa basic usage works", {
  res <- task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
  expect_equal(length(res$U_aligned_list), N)
})

test_that("task_hatsa_opts works correctly", {
  # Create custom options
  opts <- task_hatsa_opts(
    lambda_blend_value = 0.2,
    k_conn_pos = 5,
    k_conn_neg = 5,
    n_refine = 3
  )
  
  # Test with the options
  res <- task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = task_data_list,
    task_method = "lambda_blend",
    opts = opts,
    verbose = FALSE
  )
  
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
})

test_that("task_hatsa works with ... arguments", {
  # Test passing options directly via ...
  res <- task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = task_data_list,
    task_method = "lambda_blend",
    lambda_blend_value = 0.3,  # via ...
    k_conn_pos = 6,            # via ...
    verbose = FALSE
  )
  
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
})

test_that("task_hatsa handles graph_mode parameter", {
  res <- task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_method = "core_hatsa",
    graph_mode = "schur_complement",  # Explicit graph_mode
    verbose = FALSE
  )
  
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
})

test_that("task_hatsa works with gev_patch method", {
  res <- task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = task_data_list,
    task_method = "gev_patch",
    k_gev_dims = 2,  # via ...
    verbose = FALSE
  )
  
  expect_type(res, "list")
  expect_true("U_aligned_list" %in% names(res))
}) 