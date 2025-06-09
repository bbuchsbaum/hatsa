library(testthat)

create_quality_projector <- function(N = 3, V = 5, k = 2, anchor_indices = 1:2) {
  U_list <- lapply(seq_len(N), function(i) {
    matrix(seq((i - 1) * V * k + 1, i * V * k), nrow = V, ncol = k, byrow = TRUE)
  })
  R_list <- replicate(N, diag(k), simplify = FALSE)
  U_aligned <- Map(function(U, R) U %*% R, U_list, R_list)
  v <- Reduce("+", U_aligned) / N
  s <- do.call(rbind, U_aligned)
  block_idx <- split(seq_len(N * V), rep(seq_len(N), each = V))
  Lambda_list <- replicate(N, rep(1, k), simplify = FALSE)
  params <- list(k = k, N_subjects = N, V_p = V, anchor_indices = anchor_indices)
  obj <- list(
    v = v,
    s = s,
    sdev = rep(1, k),
    preproc = NULL,
    block_indices = block_idx,
    R_final_list = R_list,
    U_original_list = U_list,
    Lambda_original_list = Lambda_list,
    Lambda_original_gaps_list = replicate(N, rep(1, k - 1), simplify = FALSE),
    T_anchor_final = v[anchor_indices, , drop = FALSE],
    parameters = params,
    method = "hatsa_core",
    ._cache = list()
  )
  class(obj) <- c("hatsa_projector", "multiblock_biprojector", "projector", "list")
  obj
}

test_that("get_anchor_sketches returns anchor rows for each subject", {
  proj <- create_quality_projector()
  a_orig <- get_anchor_sketches(proj, type = "original")
  expect_length(a_orig, proj$parameters$N_subjects)
  expect_equal(a_orig[[1]], proj$U_original_list[[1]][proj$parameters$anchor_indices, , drop = FALSE])
  a_aligned <- get_anchor_sketches(proj, type = "aligned")
  expect_equal(a_aligned, a_orig)
  single <- get_anchor_sketches(proj, type = "original", subject_idx = 2)
  expect_equal(single, a_orig[[2]])
})

test_that("parcel_quality_metrics computes reconstruction and variance", {
  proj <- create_quality_projector()
  recon <- parcel_quality_metrics(proj, type = "reconstruction")
  expect_equal(dim(recon), c(proj$parameters$V_p, proj$parameters$N_subjects))
  expect_true(all(recon == 0))
  varmet <- parcel_quality_metrics(proj, type = "variance_captured")
  expected_var <- cbind(
    rowSums(proj$U_original_list[[1]]^2),
    rowSums(proj$U_original_list[[2]]^2),
    rowSums(proj$U_original_list[[3]]^2)
  )
  expect_equal(varmet, expected_var)
  expect_error(parcel_quality_metrics(proj, type = "connectivity_preservation"),
               "not yet implemented")
})

test_that("get_condition_numbers and summary_enhanced compute expected values", {
  proj <- create_quality_projector()
  cond <- get_condition_numbers(proj, what = "all")
  anchor_cond <- sapply(get_anchor_sketches(proj, type = "original"),
                        function(A) kappa(A, exact = TRUE))
  expect_equal(cond$kappa_T_anchor_final,
               kappa(proj$T_anchor_final, exact = TRUE))
  expect_equal(cond$kappa_anchor_sketches_per_subject, anchor_cond)
  expect_equal(cond$kappa_anchor_sketches_mean, mean(anchor_cond))
  expect_equal(cond$kappa_anchor_sketches_max, max(anchor_cond))

  sum_enh <- summary_enhanced(proj)
  expect_s3_class(sum_enh, "summary_enhanced_hatsa")
  expect_true(is.list(sum_enh$condition_numbers))
  expect_equal(sum_enh$condition_numbers$kappa_T_anchor_final,
               cond$kappa_T_anchor_final)
})
