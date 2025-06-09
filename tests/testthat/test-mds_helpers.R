library(testthat)

# Tests for functions in mds_helpers.R


test_that("run_cmdscale_safe returns cmdscale output and handles errors", {
  dist_mat <- matrix(c(
    0, 1, 2,
    1, 0, 1,
    2, 1, 0
  ), nrow = 3, byrow = TRUE)
  res_safe <- run_cmdscale_safe(dist_mat, k_mds = 2, add = TRUE)
  res_base <- cmdscale(as.dist(dist_mat), k = 2, eig = TRUE, add = TRUE)
  expect_equal(res_safe$points, res_base$points, tolerance = 1e-8)
  expect_equal(res_safe$eig, res_base$eig)

  bad_mat <- matrix(c(0, 1, 2, 3), nrow = 2)
  expect_warning(res_err <- run_cmdscale_safe(bad_mat), "MDS computation")
  expect_null(res_err)
})


test_that("compute_mds_distance_matrix delegates to generic and handles errors", {
  dummy_obj <- structure(list(), class = "dummy_proj")
  riemannian_distance_matrix_spd.dummy_proj <- function(object, type, verbose = FALSE) {
    matrix(c(0, 1, 1, 0), nrow = 2)
  }
  res <- compute_mds_distance_matrix(dummy_obj, "cov_coeffs", verbose = FALSE)
  expect_true(is.matrix(res))
  expect_equal(res[1, 2], 1)

  riemannian_distance_matrix_spd.dummy_proj <- function(object, type, verbose = FALSE) {
    stop("fail")
  }
  expect_warning(res_err <- compute_mds_distance_matrix(dummy_obj, "cov_coeffs"),
                 "distance matrix computation failed")
  expect_null(res_err)
  rm(riemannian_distance_matrix_spd.dummy_proj, envir = .GlobalEnv)
})


test_that("prepare_mds_plot_df merges subject info and labels", {
  mds_fit <- list(points = matrix(c(0, 1, 2, 3), ncol = 2))
  valid_idx <- c(2, 4)
  mask <- c(FALSE, TRUE, FALSE, TRUE)
  subj_info <- data.frame(
    age = c(25, 30, 35, 40),
    group = c("A", "B", "A", "B"),
    subject_label = paste0("S", 1:4)
  )

  res <- prepare_mds_plot_df(mds_fit, valid_idx, mask,
                             n_total_subjects = 4,
                             subject_info = subj_info,
                             plot_labels = TRUE)
  expect_equal(nrow(res), 2)
  expect_equal(res$subject_idx_original, valid_idx)
  expect_equal(res$subject_label_plot, subj_info$subject_label[valid_idx])
  expect_true(all(c("age", "group", "subject_label") %in% colnames(res)))

  subj_info_bad <- subj_info[1:3, ]
  expect_warning(res_bad <- prepare_mds_plot_df(mds_fit, valid_idx, mask,
                          n_total_subjects = 4, subject_info = subj_info_bad),
                 "row count")
  expect_equal(colnames(res_bad), c("Dim1", "Dim2", "subject_idx_original", "subject_label_plot"))
})
