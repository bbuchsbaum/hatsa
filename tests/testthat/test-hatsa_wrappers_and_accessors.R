library(testthat)
skip_on_cran()

.gen_subject_data <- function(N, V, T_len = 20) {
  lapply(seq_len(N), function(i) matrix(rnorm(T_len * V), ncol = V))
}

.gen_task_data <- function(N, T_len, C, V) {
  # Task data should be C x V (conditions x parcels) for activation matrices
  lapply(seq_len(N), function(i) matrix(rnorm(C * V), nrow = C, ncol = V))
}

# Test hatsa wrapper with preset parameters

test_that("hatsa wrapper applies preset values", {
  set.seed(1)
  subj <- .gen_subject_data(2, 6, 20)
  res <- suppressMessages(hatsa(subj, anchors = 1:3, components = 2, preset = "fast"))
  expect_s3_class(res, "hatsa_projector")
  expect_equal(res$parameters$k_conn_pos, 5)
  expect_equal(res$parameters$n_refine, 1)
  expect_equal(res$parameters$V_p, 6)
  expect_length(res$parameters$anchor_indices, 3)
})

# Test task_hatsa automatic method selection

test_that("task_hatsa selects method automatically", {
  set.seed(2)
  subj <- .gen_subject_data(2, 5, 15)
  task <- .gen_task_data(2, 15, 12, 5)  # >10 conditions -> GEV
  
  # Create opts with fast preset
  opts <- task_hatsa_opts(
    k_conn_pos = 5,
    k_conn_neg = 2,
    n_refine = 1
  )
  
  res <- suppressMessages(task_hatsa(
    subject_data_list = subj,
    task_data_list = task,
    anchor_indices = 1:2,
    spectral_rank_k = 2,
    task_method = "gev_patch",  # Since we know it should be GEV
    opts = opts,
    verbose = FALSE
  ))
  expect_s3_class(res, "task_hatsa_projector")
  expect_equal(res$parameters$task_method, "gev_patch")
})

# Test accessor functions on hatsa output

test_that("hatsa accessor helpers return expected structures", {
  set.seed(3)
  subj <- .gen_subject_data(3, 4, 18)
  obj <- suppressMessages(hatsa(subj, anchors = 1:2, components = 2, preset = "fast"))

  aligned <- get_aligned_data(obj)
  expect_equal(length(aligned), 3)
  expect_true(all(vapply(aligned, is.matrix, logical(1))))
  expect_equal(dim(aligned[[1]]), c(4, 2))

  template <- get_template(obj)
  expect_equal(dim(template), c(length(obj$parameters$anchor_indices), 2))

  rots <- get_rotations(obj)
  expect_equal(length(rots), 3)
  expect_true(all(vapply(rots, function(m) is.matrix(m) && all(dim(m)==c(2,2)), logical(1))))

  anchors <- get_anchor_indices(obj)
  expect_equal(anchors, obj$parameters$anchor_indices)

  qmet <- get_quality_metrics(obj)
  expect_named(qmet, c("mean_anchor_error", "rotation_dispersion", "eigenvalue_gaps",
                       "n_subjects", "n_components", "n_voxels"))

  out <- capture.output(met <- hatsa_summary(obj))
  expect_true(is.list(met))
  expect_equal(met$n_subjects, 3)
})

