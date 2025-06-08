# Test for end-to-end functionality of run_task_hatsa
# THFIX-007

context("Integration Test for run_task_hatsa")

test_that("run_task_hatsa with core_hatsa method runs and produces valid output (sequential)", {
  # Ensure future plan is sequential for this test block for deterministic behavior
  # If future is used internally by run_task_hatsa, it should respect this.
  # If run_task_hatsa itself sets a plan, this might be overridden.
  # For now, assume we control the top-level plan or run_task_hatsa accepts it.
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  # 1. Synthetic Data Generation
  N_subjects_test <- 3
  V_p_test <- 10 # Number of parcels
  T_subj_test <- 50 # Time points per subject
  k_test <- 3 # Spectral rank
  num_anchors_test <- 5

  subject_data_list_test <- vector("list", N_subjects_test)
  parcel_names_test <- paste0("P", 1:V_p_test)
  for (i in 1:N_subjects_test) {
    mat <- matrix(rnorm(T_subj_test * V_p_test), nrow = T_subj_test, ncol = V_p_test)
    colnames(mat) <- parcel_names_test
    subject_data_list_test[[i]] <- mat
  }

  # 2. Parameters for run_task_hatsa
  anchor_indices_test <- 1:num_anchors_test
  
  # Ensure all necessary parameters for run_task_hatsa are provided.
  # It's assumed run_task_hatsa is available in the environment (e.g., via devtools::load_all())
  
  # We'll call expect_no_error to catch any errors during the run
  result_core_hatsa <- NULL
  expect_no_error({
    result_core_hatsa <- run_task_hatsa(
      subject_data_list = subject_data_list_test,
      task_data_list = NULL, # For core_hatsa method
      anchor_indices = anchor_indices_test,
      spectral_rank_k = k_test,
      k_conn_pos = 3,
      k_conn_neg = 3,
      n_refine = 2,
      task_method = "core_hatsa", # Single string (will be matched against c("lambda_blend", "gev_patch", "core_hatsa"))
      lambda_blend_value = 0.5, # Default, not used by core_hatsa
      row_augmentation = FALSE, # Explicitly disable row augmentation when task_data_list is NULL
      omega_mode = "fixed", # Single string (will be matched against c("fixed", "adaptive"))
      verbose = FALSE # Keep tests quiet
    )
  })

  # 3. Assertions
  expect_true(!is.null(result_core_hatsa), "run_task_hatsa should return a non-NULL object.")
  expect_s3_class(result_core_hatsa, "task_hatsa_projector")
  expect_s3_class(result_core_hatsa, "hatsa_projector") # Inherits

  # Check parameters stored in the object
  expect_equal(result_core_hatsa$parameters$k, k_test)
  expect_equal(result_core_hatsa$parameters$N_subjects, N_subjects_test)
  expect_equal(result_core_hatsa$parameters$V_p, V_p_test)
  expect_equal(length(result_core_hatsa$parameters$anchor_indices), num_anchors_test)
  expect_equal(result_core_hatsa$parameters$task_method, "core_hatsa")

  # Check dimensions of key components
  # $v: group-level template (V_p x k)
  expect_true(is.matrix(result_core_hatsa$v), "result$v should be a matrix")
  expect_equal(dim(result_core_hatsa$v), c(V_p_test, k_test))

  # $s: stacked aligned sketches ((N_subjects * V_p) x k)
  expect_true(is.matrix(result_core_hatsa$s), "result$s should be a matrix")
  expect_equal(dim(result_core_hatsa$s), c(N_subjects_test * V_p_test, k_test))
  expect_false(any(is.na(result_core_hatsa$s)), "Stacked scores 's' should not contain NAs for core_hatsa with valid inputs.")


  # $R_final_list: list of rotation matrices (k x k)
  expect_true(is.list(result_core_hatsa$R_final_list), "result$R_final_list should be a list")
  expect_equal(length(result_core_hatsa$R_final_list), N_subjects_test)
  for (i in 1:N_subjects_test) {
    expect_true(is.matrix(result_core_hatsa$R_final_list[[i]]), paste("R_final_list[[ D ]]", i, "is not a matrix"))
    expect_equal(dim(result_core_hatsa$R_final_list[[i]]), c(k_test, k_test))
  }

  # $U_original_list: list of original sketches (V_p x k)
  expect_true(is.list(result_core_hatsa$U_original_list), "result$U_original_list should be a list")
  expect_equal(length(result_core_hatsa$U_original_list), N_subjects_test)
  for (i in 1:N_subjects_test) {
    expect_true(is.matrix(result_core_hatsa$U_original_list[[i]]), paste("U_original_list[[ D ]]",i ,"is not a matrix"))
    expect_equal(dim(result_core_hatsa$U_original_list[[i]]), c(V_p_test, k_test))
  }
  
  # $U_aligned_list (stored in task_hatsa_projector from task_hatsa_helpers)
  expect_true(is.list(result_core_hatsa$U_aligned_list), "result$U_aligned_list should be a list")
  expect_equal(length(result_core_hatsa$U_aligned_list), N_subjects_test)
  for (i in 1:N_subjects_test) {
    expect_true(is.matrix(result_core_hatsa$U_aligned_list[[i]]), paste("U_aligned_list[[ D ]]",i ,"is not a matrix"))
    expect_equal(dim(result_core_hatsa$U_aligned_list[[i]]), c(V_p_test, k_test))
    # Check consistency with stacked scores 's'
    s_block <- result_core_hatsa$s[((i-1)*V_p_test + 1):(i*V_p_test), , drop = FALSE]
    expect_equal(result_core_hatsa$U_aligned_list[[i]], s_block, 
                 info = paste("U_aligned_list[[ D ]]", i, "should match corresponding block in s_stacked"))

  }


  # $T_anchor_final: group anchor template (num_anchors x k)
  expect_true(is.matrix(result_core_hatsa$T_anchor_final), "result$T_anchor_final should be a matrix")
  expect_equal(dim(result_core_hatsa$T_anchor_final), c(num_anchors_test, k_test))
  
  # $Lambda_original_list
  expect_true(is.list(result_core_hatsa$Lambda_original_list), "result$Lambda_original_list should be a list")
  expect_equal(length(result_core_hatsa$Lambda_original_list), N_subjects_test)
  for (i in 1:N_subjects_test) {
    expect_true(is.numeric(result_core_hatsa$Lambda_original_list[[i]]), paste("Lambda_original_list[[ D ]]",i ,"is not numeric"))
    expect_equal(length(result_core_hatsa$Lambda_original_list[[i]]), k_test)
  }

  # Check that task-specific components are NULL or empty for core_hatsa
  expect_true(is.null(result_core_hatsa$W_task_list) || 
              (is.list(result_core_hatsa$W_task_list) && length(result_core_hatsa$W_task_list) == 0) ||
              all(sapply(result_core_hatsa$W_task_list, is.null)))
              
  expect_true(is.null(result_core_hatsa$L_task_list) || 
              (is.list(result_core_hatsa$L_task_list) && length(result_core_hatsa$L_task_list) == 0) ||
              all(sapply(result_core_hatsa$L_task_list, is.null)))
              
  expect_true(is.null(result_core_hatsa$W_hybrid_list) || 
              (is.list(result_core_hatsa$W_hybrid_list) && length(result_core_hatsa$W_hybrid_list) == 0) ||
              all(sapply(result_core_hatsa$W_hybrid_list, is.null)))
              
  expect_true(is.null(result_core_hatsa$L_hybrid_list) || 
              (is.list(result_core_hatsa$L_hybrid_list) && length(result_core_hatsa$L_hybrid_list) == 0) ||
              all(sapply(result_core_hatsa$L_hybrid_list, is.null)))
              
  expect_true(is.null(result_core_hatsa$U_task_list) || 
              (is.list(result_core_hatsa$U_task_list) && length(result_core_hatsa$U_task_list) == 0) ||
              all(sapply(result_core_hatsa$U_task_list, is.null)))
              
  expect_true(is.null(result_core_hatsa$Lambda_task_list) || 
              (is.list(result_core_hatsa$Lambda_task_list) && length(result_core_hatsa$Lambda_task_list) == 0) ||
              all(sapply(result_core_hatsa$Lambda_task_list, is.null)))

})

test_that("run_task_hatsa with lambda_blend method runs and produces valid output (sequential)", {
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    skip("Matrix package not available, skipping lambda_blend test that uses sparse matrices.")
  }

  # 1. Synthetic Data Generation (similar to core_hatsa test)
  N_subjects_test <- 2 # Reduced for quicker test
  V_p_test <- 8
  T_subj_test <- 40
  k_test <- 2
  num_anchors_test <- 4

  subject_data_list_test <- vector("list", N_subjects_test)
  parcel_names_test <- paste0("P", 1:V_p_test)
  for (i in 1:N_subjects_test) {
    mat <- matrix(rnorm(T_subj_test * V_p_test), nrow = T_subj_test, ncol = V_p_test)
    colnames(mat) <- parcel_names_test
    subject_data_list_test[[i]] <- mat
  }

  # Mock task_data_list: Each element is a list containing W_task_i for that subject
  # W_task_i should be a V_p x V_p sparse matrix
  task_data_list_test <- vector("list", N_subjects_test)
  for (i in 1:N_subjects_test) {
    # Create a simple sparse positive definite-like matrix for W_task_i
    diag_vals <- runif(V_p_test, 0.5, 1.5)
    off_diag_count <- floor(V_p_test * V_p_test * 0.1) # 10% sparsity
    row_indices <- sample(1:V_p_test, off_diag_count, replace = TRUE)
    col_indices <- sample(1:V_p_test, off_diag_count, replace = TRUE)
    # Ensure symmetry for off-diagonals
    all_rows <- c(row_indices, col_indices, 1:V_p_test)
    all_cols <- c(col_indices, row_indices, 1:V_p_test) # Symmetrize
    all_vals <- c(runif(off_diag_count, -0.2, 0.2), runif(off_diag_count, -0.2, 0.2), diag_vals)
    
    W_task_mat_sparse <- Matrix::sparseMatrix(i = all_rows, j = all_cols, x = all_vals, 
                                        dims = c(V_p_test, V_p_test),
                                        symmetric = FALSE) # Build then forceSymmetric
    W_task_mat_sparse <- Matrix::forceSymmetric(W_task_mat_sparse, uplo="L") # Ensure it's symmetric after creation
    # Ensure positive definiteness by adding to diagonal if needed (simplified approach)
    # A more robust way is to ensure it's a valid covariance or regularized Laplacian.
    # For testing, we hope the construction of W_hybrid results in something reasonable.
    # What shape_basis expects for lambda_blend is W_task_i (z-scored). 
    # The z-scoring happens inside compute_task_matrices based on task_data_type.
    # For simplicity, we provide a W_task_i directly if task_matrix_method="from_precomputed_Wtask"
    # Let's assume task_matrix_method = "from_features" and provide some mock task features.
    # Each task_data_list[[i]] should be a list of matrices, where each matrix is features_c x V_p
    num_task_conditions <- 2
    num_task_features <- 5 # e.g. 5 features per condition
    task_features_subj_i <- vector("list", num_task_conditions)
    for(cond in 1:num_task_conditions){
        task_features_subj_i[[cond]] <- matrix(rnorm(num_task_features * V_p_test), 
                                                 nrow=num_task_features, ncol=V_p_test)
    }
    task_data_list_test[[i]] <- task_features_subj_i
  }

  anchor_indices_test <- 1:num_anchors_test
  result_lambda_blend <- NULL

  expect_no_error({
    result_lambda_blend <- run_task_hatsa(
      subject_data_list = subject_data_list_test,
      task_data_list = task_data_list_test, 
      anchor_indices = anchor_indices_test,
      spectral_rank_k = k_test,
      k_conn_pos = 3,
      k_conn_neg = 3,
      n_refine = 1, # Reduced for quicker test
      task_method = "lambda_blend", # Single string
      lambda_blend_value = 0.5,
      omega_mode = "fixed", # Single string
      row_augmentation = TRUE, 
      scale_omega_trace = TRUE, 
      verbose = FALSE
    )
  })

  expect_true(!is.null(result_lambda_blend), "run_task_hatsa (lambda_blend) should return non-NULL.")
  expect_s3_class(result_lambda_blend, "task_hatsa_projector")

  # Basic parameter checks
  expect_equal(result_lambda_blend$parameters$k, k_test)
  expect_equal(result_lambda_blend$parameters$N_subjects, N_subjects_test)
  expect_equal(result_lambda_blend$parameters$V_p, V_p_test)
  expect_equal(result_lambda_blend$parameters$task_method, "lambda_blend")
  expect_equal(result_lambda_blend$parameters$lambda_blend_value, 0.5)

  # Check dimensions of core components (v, s, R_final_list, U_original_list, T_anchor_final)
  # (Similar checks as in core_hatsa test)
  expect_true(is.matrix(result_lambda_blend$v) && all(dim(result_lambda_blend$v) == c(V_p_test, k_test)))
  expect_true(is.matrix(result_lambda_blend$s) && all(dim(result_lambda_blend$s) == c(N_subjects_test * V_p_test, k_test)))
  expect_false(any(is.na(result_lambda_blend$s)), "Stacked scores 's' should not contain NAs for lambda_blend with valid inputs.")

  expect_true(is.list(result_lambda_blend$R_final_list) && length(result_lambda_blend$R_final_list) == N_subjects_test)
  expect_true(all(sapply(result_lambda_blend$R_final_list, function(m) is.matrix(m) && all(dim(m) == c(k_test, k_test)))))
  
  expect_true(is.list(result_lambda_blend$U_original_list) && length(result_lambda_blend$U_original_list) == N_subjects_test)
  expect_true(all(sapply(result_lambda_blend$U_original_list, function(m) is.matrix(m) && all(dim(m) == c(V_p_test, k_test)))))

  expect_true(is.matrix(result_lambda_blend$T_anchor_final) && all(dim(result_lambda_blend$T_anchor_final) == c(num_anchors_test, k_test)))
  
  # Check presence and dimensions of task-specific or hybrid components
  # For lambda_blend, we expect W_hybrid_list, L_hybrid_list, U_hybrid_list, Lambda_hybrid_list
  # Or, the U_aligned_list would be the U_hybrid_list essentially.
  # The constructor task_hatsa_projector stores U_aligned_list which contains the final sketches (hybrid in this case)

  # W_task_list check - allow for NULL, empty list, or list of NULLs
  expect_true(
    is.null(result_lambda_blend$W_task_list) ||
    (is.list(result_lambda_blend$W_task_list) && 
     (length(result_lambda_blend$W_task_list) == N_subjects_test || 
      (length(result_lambda_blend$W_task_list) > 0 && 
       all(sapply(result_lambda_blend$W_task_list, function(w) is.null(w) || inherits(w, "Matrix"))))))
  )

  # W_hybrid_list check
  expect_true(
    is.null(result_lambda_blend$W_hybrid_list) ||
    (is.list(result_lambda_blend$W_hybrid_list) && 
     (length(result_lambda_blend$W_hybrid_list) == 0 || 
      all(sapply(result_lambda_blend$W_hybrid_list, function(w) is.null(w) || inherits(w, "Matrix")))))
  )

  # L_hybrid_list check
  expect_true(
    is.null(result_lambda_blend$L_hybrid_list) ||
    (is.list(result_lambda_blend$L_hybrid_list) && 
     (length(result_lambda_blend$L_hybrid_list) == 0 || 
      all(sapply(result_lambda_blend$L_hybrid_list, function(l) is.null(l) || inherits(l, "Matrix")))))
  )

  # U_aligned_list check
  expect_true(
    is.null(result_lambda_blend$U_aligned_list) ||
    (is.list(result_lambda_blend$U_aligned_list) && 
     length(result_lambda_blend$U_aligned_list) == N_subjects_test)
  )

  # Check original connectivity-based sketches
  expect_true(
    is.null(result_lambda_blend$U_conn_list) ||
    (is.list(result_lambda_blend$U_conn_list) && 
     (length(result_lambda_blend$U_conn_list) == 0 || 
      all(sapply(result_lambda_blend$U_conn_list, function(u) is.null(u) || is.matrix(u)))))
  )
  
  expect_true(
    is.null(result_lambda_blend$Lambda_conn_list) ||
    (is.list(result_lambda_blend$Lambda_conn_list) && 
     (length(result_lambda_blend$Lambda_conn_list) == 0 || 
      all(sapply(result_lambda_blend$Lambda_conn_list, function(l) is.null(l) || is.numeric(l)))))
  )

})

# TODO THFIX-007 Task 3: Test with future_plan = "multisession" if appropriate and stable

test_that("run_task_hatsa with core_hatsa method runs with future_plan = 'multisession'", {
  if (!requireNamespace("future", quietly = TRUE)) {
    skip("future package not available, skipping multisession test.")
  }
  if (future::availableCores() < 2) {
    skip("Less than 2 cores available, skipping multisession test as it might behave like sequential.")
  }

  old_plan <- future::plan(future::multisession) # Set to multisession for this test
  on.exit(future::plan(old_plan), add = TRUE)

  # 1. Synthetic Data Generation (same as the first core_hatsa sequential test)
  N_subjects_test <- 3
  V_p_test <- 10 
  T_subj_test <- 50 
  k_test <- 3 
  num_anchors_test <- 5

  subject_data_list_test <- vector("list", N_subjects_test)
  parcel_names_test <- paste0("P", 1:V_p_test)
  for (i in 1:N_subjects_test) {
    mat <- matrix(rnorm(T_subj_test * V_p_test), nrow = T_subj_test, ncol = V_p_test)
    colnames(mat) <- parcel_names_test
    subject_data_list_test[[i]] <- mat
  }
  anchor_indices_test <- 1:num_anchors_test

  result_multisession <- NULL
  expect_no_error({
    result_multisession <- run_task_hatsa(
      subject_data_list = subject_data_list_test,
      task_data_list = NULL, 
      anchor_indices = anchor_indices_test,
      spectral_rank_k = k_test,
      k_conn_pos = 3,
      k_conn_neg = 3,
      n_refine = 2,
      task_method = "core_hatsa", # Single string
      lambda_blend_value = 0.5,
      omega_mode = "fixed", # Single string
      row_augmentation = FALSE, # Explicitly disable row augmentation when task_data_list is NULL 
      scale_omega_trace = TRUE, 
      verbose = FALSE
    )
  })

  # 3. Assertions (identical structure to the sequential core_hatsa test)
  expect_true(!is.null(result_multisession), "run_task_hatsa (multisession) should return non-NULL.")
  expect_s3_class(result_multisession, "task_hatsa_projector")
  expect_s3_class(result_multisession, "hatsa_projector")

  expect_equal(result_multisession$parameters$k, k_test)
  expect_equal(result_multisession$parameters$N_subjects, N_subjects_test)
  expect_equal(result_multisession$parameters$V_p, V_p_test)
  expect_equal(length(result_multisession$parameters$anchor_indices), num_anchors_test)
  expect_equal(result_multisession$parameters$task_method, "core_hatsa")

  expect_true(is.matrix(result_multisession$v) && all(dim(result_multisession$v) == c(V_p_test, k_test)))
  expect_true(is.matrix(result_multisession$s) && all(dim(result_multisession$s) == c(N_subjects_test * V_p_test, k_test)))
  # For multisession, NA check might be too strict if minor numerical differences lead to issues previously unseen.
  # However, with core_hatsa and simple data, it should ideally still be NA-free if sequential was.
  # If this fails, it might indicate an issue with how data is handled/aggregated in parallel.
  expect_false(any(is.na(result_multisession$s)), "Stacked scores 's' (multisession) should not contain NAs for core_hatsa with valid inputs.")

  expect_true(is.list(result_multisession$R_final_list) && length(result_multisession$R_final_list) == N_subjects_test)
  expect_true(all(sapply(result_multisession$R_final_list, function(m) is.matrix(m) && all(dim(m) == c(k_test, k_test)))))
  
  expect_true(is.list(result_multisession$U_original_list) && length(result_multisession$U_original_list) == N_subjects_test)
  expect_true(all(sapply(result_multisession$U_original_list, function(m) is.matrix(m) && all(dim(m) == c(V_p_test, k_test)))))

  expect_true(is.list(result_multisession$U_aligned_list) && length(result_multisession$U_aligned_list) == N_subjects_test)
  expect_true(all(sapply(result_multisession$U_aligned_list, function(u) is.matrix(u) && all(dim(u) == c(V_p_test, k_test)))))
  for (i in 1:N_subjects_test) {
    s_block <- result_multisession$s[((i-1)*V_p_test + 1):(i*V_p_test), , drop = FALSE]
    expect_equal(result_multisession$U_aligned_list[[i]], s_block, 
                 info = paste("U_aligned_list[[ D ]] (multisession)", i, "should match corresponding block in s_stacked"))
  }

  expect_true(is.matrix(result_multisession$T_anchor_final) && all(dim(result_multisession$T_anchor_final) == c(num_anchors_test, k_test)))
  
  expect_true(is.list(result_multisession$Lambda_original_list) && length(result_multisession$Lambda_original_list) == N_subjects_test)
  expect_true(all(sapply(result_multisession$Lambda_original_list, function(l) is.numeric(l) && length(l) == k_test)))

  expect_true(is.null(result_multisession$W_task_list) || 
              (is.list(result_multisession$W_task_list) && length(result_multisession$W_task_list) == 0) ||
              all(sapply(result_multisession$W_task_list, is.null)))
              
  expect_true(is.null(result_multisession$L_task_list) || 
              (is.list(result_multisession$L_task_list) && length(result_multisession$L_task_list) == 0) ||
              all(sapply(result_multisession$L_task_list, is.null)))
              
  expect_true(is.null(result_multisession$W_hybrid_list) || 
              (is.list(result_multisession$W_hybrid_list) && length(result_multisession$W_hybrid_list) == 0) ||
              all(sapply(result_multisession$W_hybrid_list, is.null)))
              
  expect_true(is.null(result_multisession$L_hybrid_list) || 
              (is.list(result_multisession$L_hybrid_list) && length(result_multisession$L_hybrid_list) == 0) ||
              all(sapply(result_multisession$L_hybrid_list, is.null)))
              
  expect_true(is.null(result_multisession$U_task_list) || 
              (is.list(result_multisession$U_task_list) && length(result_multisession$U_task_list) == 0) ||
              all(sapply(result_multisession$U_task_list, is.null)))
              
  expect_true(is.null(result_multisession$Lambda_task_list) || 
              (is.list(result_multisession$Lambda_task_list) && length(result_multisession$Lambda_task_list) == 0) ||
              all(sapply(result_multisession$Lambda_task_list, is.null)))

})

