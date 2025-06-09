describe("Voxel Projection Functionality", {

# Helper function to generate mock coordinates
.generate_mock_coords <- function(N_points, N_dim = 3) {
  matrix(rnorm(N_points * N_dim), ncol = N_dim)
}

# Helper function to generate mock spectral components (U_orig_parcel, Lambda_orig_parcel)
.generate_mock_parcel_components <- function(V_p, k) {
  list(
    U_orig_parcel = matrix(rnorm(V_p * k), ncol = k),
    Lambda_orig_parcel = sort(runif(k, 0.1, 5), decreasing = TRUE) # Ensure positive and sorted
  )
}

# Helper function to generate mock voxel time-series list
.generate_mock_voxel_ts_list <- function(N_subjects, V_v, T_i_mean = 50, T_i_sd = 5) {
  lapply(1:N_subjects, function(i) {
    matrix(rnorm(round(rnorm(1, T_i_mean, T_i_sd)) * V_v), ncol = V_v)
  })
}

# Simplified helper from the other test file (or could be shared if test helpers are centralized)
.generate_mock_subject_data_parcels <- function(N_subjects, V_p, T_i_mean = 100, T_i_sd = 10) {
  lapply(1:N_subjects, function(i) {
    matrix(rnorm(round(rnorm(1, T_i_mean, T_i_sd)) * V_p), ncol = V_p)
  })
}

.get_default_hatsa_params_for_voxel_test <- function(V_p, k) {
  list(
    anchor_indices = sample(1:V_p, min(V_p, k + 1)), # Ensure enough for anchors
    spectral_rank_k = k,
    k_conn_pos = min(5, V_p -1),
    k_conn_neg = min(5, V_p-1),
    n_refine = 1 # Keep low for speed
  )
}


test_that("compute_voxel_basis_nystrom: basic functionality and dimensions", {
  V_p <- 30  # Number of parcels
  V_v <- 100 # Number of voxels
  k <- 5    # Number of components

  parcel_coords <- .generate_mock_coords(V_p)
  voxel_coords <- .generate_mock_coords(V_v)
  parcel_comps <- .generate_mock_parcel_components(V_p, k)

  # Test with default row_normalize_W = FALSE
  phi_voxel <- compute_voxel_basis_nystrom(
    voxel_coords = voxel_coords,
    parcel_coords = parcel_coords,
    U_orig_parcel = parcel_comps$U_orig_parcel,
    Lambda_orig_parcel = parcel_comps$Lambda_orig_parcel,
    n_nearest_parcels = 5,
    kernel_sigma = 5.0
  )

  expect_true(is.matrix(phi_voxel))
  expect_equal(nrow(phi_voxel), V_v)
  expect_equal(ncol(phi_voxel), k)
  expect_true(all(is.finite(phi_voxel)))

  # Test with row_normalize_W = TRUE
  phi_voxel_norm <- compute_voxel_basis_nystrom(
    voxel_coords = voxel_coords,
    parcel_coords = parcel_coords,
    U_orig_parcel = parcel_comps$U_orig_parcel,
    Lambda_orig_parcel = parcel_comps$Lambda_orig_parcel,
    n_nearest_parcels = 5,
    kernel_sigma = 5.0,
    row_normalize_W = TRUE
  )
  expect_true(is.matrix(phi_voxel_norm))
  expect_equal(nrow(phi_voxel_norm), V_v)
  expect_equal(ncol(phi_voxel_norm), k)
  expect_true(all(is.finite(phi_voxel_norm)))
  
  # Test kernel_sigma = "auto"
  phi_voxel_auto_sigma <- compute_voxel_basis_nystrom(
    voxel_coords = voxel_coords,
    parcel_coords = parcel_coords,
    U_orig_parcel = parcel_comps$U_orig_parcel,
    Lambda_orig_parcel = parcel_comps$Lambda_orig_parcel,
    n_nearest_parcels = 5,
    kernel_sigma = "auto"
  )
  expect_true(is.matrix(phi_voxel_auto_sigma))
  expect_equal(nrow(phi_voxel_auto_sigma), V_v)
  expect_equal(ncol(phi_voxel_auto_sigma), k)

  # Edge case: k=0
  parcel_comps_k0 <- .generate_mock_parcel_components(V_p, 0)
  phi_voxel_k0 <- compute_voxel_basis_nystrom(
    voxel_coords = voxel_coords,
    parcel_coords = parcel_coords,
    U_orig_parcel = parcel_comps_k0$U_orig_parcel, # V_p x 0 matrix
    Lambda_orig_parcel = parcel_comps_k0$Lambda_orig_parcel # empty numeric
  )
  expect_true(is.matrix(phi_voxel_k0))
  expect_equal(nrow(phi_voxel_k0), V_v)
  expect_equal(ncol(phi_voxel_k0), 0)

  # Edge case: V_v = 0
  phi_voxel_Vv0 <- compute_voxel_basis_nystrom(
    voxel_coords = .generate_mock_coords(0),
    parcel_coords = parcel_coords,
    U_orig_parcel = parcel_comps$U_orig_parcel,
    Lambda_orig_parcel = parcel_comps$Lambda_orig_parcel
  )
  expect_true(is.matrix(phi_voxel_Vv0))
  expect_equal(nrow(phi_voxel_Vv0), 0)
  expect_equal(ncol(phi_voxel_Vv0), k)
  
  # Edge case: V_p = 0 (should return V_v x k matrix of zeros)
   # This case is tricky because nn2 would fail if data has 0 rows.
   # compute_voxel_basis_nystrom has a check: `if (V_p == 0 || k == 0) return(matrix(0, nrow = V_v, ncol = k))`
   # So we need to ensure U_orig_parcel also matches this V_p=0 scenario.
  U_orig_parcel_Vp0 = matrix(0, nrow=0, ncol=k)
  phi_voxel_Vp0 <- compute_voxel_basis_nystrom(
    voxel_coords = voxel_coords,
    parcel_coords = .generate_mock_coords(0), # V_p = 0
    U_orig_parcel = U_orig_parcel_Vp0,      # 0 x k matrix
    Lambda_orig_parcel = parcel_comps$Lambda_orig_parcel # k length vector
  )
  expect_true(is.matrix(phi_voxel_Vp0))
  expect_equal(nrow(phi_voxel_Vp0), V_v)
  expect_equal(ncol(phi_voxel_Vp0), k)
  if (k > 0 && V_v > 0) expect_true(all(phi_voxel_Vp0 == 0))
  
  # Error for invalid n_nearest_parcels
  expect_error(compute_voxel_basis_nystrom(
    voxel_coords, parcel_coords, parcel_comps$U_orig_parcel, parcel_comps$Lambda_orig_parcel,
    n_nearest_parcels = 0
  ))
  # Error when n_nearest_parcels exceeds number of parcels
  expect_error(compute_voxel_basis_nystrom(
    voxel_coords, parcel_coords, parcel_comps$U_orig_parcel, parcel_comps$Lambda_orig_parcel,
    n_nearest_parcels = V_p + 1
  ))
  
  # Error for invalid kernel_sigma
  expect_error(compute_voxel_basis_nystrom(
    voxel_coords, parcel_coords, parcel_comps$U_orig_parcel, parcel_comps$Lambda_orig_parcel,
    kernel_sigma = -1
  ))
  expect_error(compute_voxel_basis_nystrom(
    voxel_coords, parcel_coords, parcel_comps$U_orig_parcel, parcel_comps$Lambda_orig_parcel,
    kernel_sigma = "invalid_string"
  ))
  
})

test_that("project_voxels.hatsa_projector: basic functionality and dimensions", {
  V_p_fit <- 20
  N_subjects_fit <- 2
  k_fit <- 3
  
  V_v_proj <- 50 

  # 1. Create a mock hatsa_projector object
  # This typically involves running run_hatsa_core or manually constructing one.
  # For simplicity, we use run_hatsa_core with small data.
  fit_parcel_data <- .generate_mock_subject_data_parcels(N_subjects_fit, V_p_fit)
  fit_params <- .get_default_hatsa_params_for_voxel_test(V_p_fit, k_fit)

  hatsa_fitted_obj <- suppressMessages(try(run_hatsa_core(
    subject_data_list = fit_parcel_data,
    anchor_indices = fit_params$anchor_indices,
    spectral_rank_k = fit_params$spectral_rank_k,
    k_conn_pos = fit_params$k_conn_pos,
    k_conn_neg = fit_params$k_conn_neg,
    n_refine = fit_params$n_refine
  ), silent = TRUE))

  if (inherits(hatsa_fitted_obj, "try-error")) {
    skip(paste0("run_hatsa_core failed during setup for project_voxels tests: ", 
                attr(hatsa_fitted_obj, "condition")$message))
    return()
  }

  # 2. Prepare inputs for project_voxels
  # Use a different number of subjects for projection to test generalizability
  N_subjects_proj <- N_subjects_fit # Must match for now due to U_orig, Lambda_orig, R retrieval
  # N_subjects_proj <- 1 # Simpler case for initial test
  
  voxel_ts_list <- .generate_mock_voxel_ts_list(N_subjects_proj, V_v_proj)
  voxel_coords_proj <- .generate_mock_coords(V_v_proj)
  # Parcel coords must match those used for the hatsa_fitted_obj. 
  # This is a bit tricky as run_hatsa_core doesn't directly store/return them.
  # For testing, we assume parcel_coords are externally known and consistent.
  # Ideally, hatsa_projector might store parcel_coords if always needed for voxel projection.
  # For now, just generate some parcel_coords with the correct V_p_fit.
  parcel_coords_fit <- .generate_mock_coords(V_p_fit)
  
  # Test .validate_coordinate_inputs directly (internal helper)
  # Using expect_message or expect_warning if applicable for its interactive messages
  # For non-interactive, it returns NULL. We can test it doesn't error.
  expect_null(.validate_coordinate_inputs(voxel_coords_proj, parcel_coords_fit))
  # Test with deliberately mismatched scales to try and trigger a message (if interactive testing was set up)
  # e.g., .validate_coordinate_inputs(voxel_coords_proj * 100, parcel_coords_fit)

  projected_voxel_data <- suppressMessages(try(project_voxels(
    object = hatsa_fitted_obj,
    voxel_timeseries_list = voxel_ts_list,
    voxel_coords = voxel_coords_proj,
    parcel_coords = parcel_coords_fit, # Must match V_p of hatsa_fitted_obj
    n_nearest_parcels = 3,
    kernel_sigma = "auto"
  ), silent = TRUE))

  if (inherits(projected_voxel_data, "try-error")) {
    fail(paste0("project_voxels.hatsa_projector failed: ", 
                attr(projected_voxel_data, "condition")$message))
    return()
  }
  
  expect_true(is.list(projected_voxel_data))
  expect_length(projected_voxel_data, N_subjects_proj)

  for (i in 1:N_subjects_proj) {
    expect_true(is.matrix(projected_voxel_data[[i]]))
    T_i_subject <- nrow(voxel_ts_list[[i]]) # Get the actual T_i for this subject
    expect_equal(nrow(projected_voxel_data[[i]]), T_i_subject)
    expect_equal(ncol(projected_voxel_data[[i]]), k_fit)
    expect_true(all(is.finite(projected_voxel_data[[i]])))
  }
  
  # Test error conditions for project_voxels
  # - voxel_timeseries_list not a list
  expect_error(project_voxels(hatsa_fitted_obj, voxel_timeseries_list = voxel_ts_list[[1]], voxel_coords_proj, parcel_coords_fit))
  # - voxel_coords not a matrix or wrong ncol
  expect_error(project_voxels(hatsa_fitted_obj, voxel_ts_list, voxel_coords = as.data.frame(voxel_coords_proj), parcel_coords_fit))
  expect_error(project_voxels(hatsa_fitted_obj, voxel_ts_list, voxel_coords = voxel_coords_proj[,1:2], parcel_coords_fit))
  # - parcel_coords not a matrix or wrong ncol or V_p mismatch
  expect_error(project_voxels(hatsa_fitted_obj, voxel_ts_list, voxel_coords_proj, parcel_coords = .generate_mock_coords(V_p_fit + 1)))
  # - voxel_timeseries_list[[i]] ncol mismatch with voxel_coords nrow
  bad_voxel_ts_list <- voxel_ts_list
  bad_voxel_ts_list[[1]] <- matrix(rnorm(10 * (V_v_proj+1)), ncol=V_v_proj+1)
  # This should ideally produce a warning and NA matrix for that subject, not a hard error for the whole call
  # The current code in project_voxels.hatsa_projector does this.
  res_bad_ts <- suppressMessages(project_voxels(hatsa_fitted_obj, bad_voxel_ts_list, voxel_coords_proj, parcel_coords_fit))
  expect_true(is.na(res_bad_ts[[1]][1,1])) # Check if it was NA-ed out
  if (length(res_bad_ts) > 1) expect_false(is.na(res_bad_ts[[2]][1,1])) # Check other subjects are fine

})

test_that("project_voxels works with precomputed W_vox_parc", {
  V_p_fit <- 15
  N_subjects_fit <- 2
  k_fit <- 3
  V_v_proj <- 40

  fit_parcel_data <- .generate_mock_subject_data_parcels(N_subjects_fit, V_p_fit)
  fit_params <- .get_default_hatsa_params_for_voxel_test(V_p_fit, k_fit)
  hatsa_fitted_obj <- suppressMessages(run_hatsa_core(
    subject_data_list = fit_parcel_data,
    anchor_indices = fit_params$anchor_indices,
    spectral_rank_k = fit_params$spectral_rank_k,
    k_conn_pos = fit_params$k_conn_pos,
    k_conn_neg = fit_params$k_conn_neg,
    n_refine = fit_params$n_refine
  ))

  voxel_ts_list <- .generate_mock_voxel_ts_list(N_subjects_fit, V_v_proj)
  voxel_coords_proj <- .generate_mock_coords(V_v_proj)
  parcel_coords_fit <- .generate_mock_coords(V_p_fit)

  # Precompute W_vox_parc similarly to compute_voxel_basis_nystrom
  nn_res <- RANN::nn2(data = parcel_coords_fit, query = voxel_coords_proj, k = 3, treetype = "kd")
  first_dists <- sqrt(nn_res$nn.dists[, 1])
  med_dist <- median(first_dists, na.rm = TRUE)
  sigma_eff <- if (is.finite(med_dist) && med_dist > 1e-6) med_dist / sqrt(2) else 1.0
  sims <- exp(-nn_res$nn.dists / (2 * sigma_eff^2))
  W_pre <- as(Matrix::sparseMatrix(
    i = rep(1:nrow(voxel_coords_proj), each = 3),
    j = as.vector(t(nn_res$nn.idx)),
    x = as.vector(t(sims)),
    dims = c(nrow(voxel_coords_proj), V_p_fit),
    repr = "T"
  ), "CsparseMatrix")

  res_default <- suppressMessages(project_voxels(
    object = hatsa_fitted_obj,
    voxel_timeseries_list = voxel_ts_list,
    voxel_coords = voxel_coords_proj,
    parcel_coords = parcel_coords_fit,
    n_nearest_parcels = 3,
    kernel_sigma = "auto"
  ))

  res_precomp <- suppressMessages(project_voxels(
    object = hatsa_fitted_obj,
    voxel_timeseries_list = voxel_ts_list,
    voxel_coords = voxel_coords_proj,
    parcel_coords = parcel_coords_fit,
    n_nearest_parcels = 3,
    kernel_sigma = "auto",
    W_vox_parc = W_pre
  ))

  expect_equal(res_precomp, res_default)
})

# Further tests for project_voxels could include:
# - Consistency check: If two subjects have identical voxel time-series and identical
#   Phi_voxel_i (which implies identical U_orig_i, Lambda_orig_i from the model, and identical coords),
#   then C_voxel_coeffs_i should be identical. C_voxel_aligned_i should then only differ by R_i.
#   This is complex to set up perfectly due to U_orig_i differing per subject.
# - A simpler consistency: if R_i is identity for all subjects (e.g. if T_anchor_final was based on subject 1 and U_orig_1)
#   then aligned and unaligned coefficients should be similar for subject 1.

})

