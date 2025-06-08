library(testthat)

# Helper to create a minimal hatsa_projector for reconstruction_error
create_mock_recon_projector <- function(N_subjects = 2, V_p = 6, k = 2, anchor_indices = c(1,2)) {
  U_orig_list <- lapply(seq_len(N_subjects), function(i) matrix(rnorm(V_p * k), V_p, k))
  R_list <- lapply(seq_len(N_subjects), function(i) {
    qr_Q(qr(matrix(rnorm(k * k), k, k)))
  })
  U_aligned <- mapply(function(U, R) U %*% R, U_orig_list, R_list, SIMPLIFY = FALSE)
  s <- do.call(rbind, U_aligned)
  block_indices <- split(seq_len(N_subjects * V_p), rep(seq_len(N_subjects), each = V_p))
  T_anchor_final <- Reduce("+", lapply(seq_len(N_subjects), function(i) U_orig_list[[i]][anchor_indices, ] %*% R_list[[i]])) / N_subjects
  obj <- list(
    s = s,
    block_indices = block_indices,
    R_final_list = R_list,
    U_original_list = U_orig_list,
    T_anchor_final = T_anchor_final,
    parameters = list(k = k, N_subjects = N_subjects, V_p = V_p, anchor_indices = anchor_indices)
  )
  class(obj) <- c("hatsa_projector", "multiblock_biprojector", "projector", "list")
  obj
}


test_that("all_parcels reconstruction yields near-zero error", {
  set.seed(123)
  proj <- create_mock_recon_projector()
  res <- reconstruction_error(proj, type = "all_parcels")
  expect_type(res, "list")
  expect_length(res$per_subject_error, proj$parameters$N_subjects)
  expect_true(all(res$per_subject_error < 1e-6))
  expect_lt(res$mean_error, 1e-6)
})


test_that("non_anchors reconstruction returns per-parcel errors", {
  set.seed(456)
  proj <- create_mock_recon_projector()
  res <- reconstruction_error(proj, type = "non_anchors")
  non_anchor_indices <- setdiff(seq_len(proj$parameters$V_p), proj$parameters$anchor_indices)
  expect_type(res, "list")
  expect_length(res$per_subject_error, proj$parameters$N_subjects)
  expect_length(res$per_parcel_error, length(non_anchor_indices))
  expect_true(all(res$per_subject_error >= 0))
  expect_true(all(res$per_parcel_error >= 0))
})
