# R/task_hatsa_helpers.R

#' Validate arguments and initialize variables for task_hatsa
#'
#' @param subject_data_list List of subject data matrices
#' @param anchor_indices Indices of anchor parcels
#' @param spectral_rank_k Spectral rank k
#' @param task_data_list List of task data
#' @param task_method Method for task incorporation
#' @param lambda_blend_value Blend value for lambda_blend task method
#' @param k_gev_dims Number of dimensions for gev_patch method
#' @param row_augmentation Whether to augment with task rows
#' @param residualize_condition_anchors Whether to residualize condition anchors
#' @param omega_weights Fixed weights for omega_mode
#' @param omega_mode Mode for omega calculation
#' @param reliability_scores_list List of reliability scores for adaptive weighting
#' @param scale_omega_trace Whether to scale omega trace
#' @param alpha_laplacian Alpha parameter for graph Laplacian
#' @param degree_type_laplacian Degree type for graph Laplacian
#' @param k_conn_pos Number of positive connections for connectivity graph
#' @param k_conn_neg Number of negative connections for connectivity graph
#' @param k_conn_task_pos Number of positive connections for task graph
#' @param k_conn_task_neg Number of negative connections for task graph
#' @param similarity_method_task Similarity method for task graph
#' @param W_task_helper_func Function to compute task graph from task data
#' @param n_refine Number of GPA refinement iterations
#' @param check_redundancy Whether to check for redundancy
#' @param redundancy_threshold Threshold for redundancy check
#' @param residualize_k_conn_proj Number of eigenvectors to remove for residualization
#' @param residualize_k_conn_labels Number of nearest neighbors to resparsify for residualization
#' @param gev_lambda_max Maximum eigenvalue for gev_patch method
#' @param gev_epsilon_reg Regularization parameter for gev_patch method
#' @param parcel_names Names of parcels
#' @param verbose Whether to print messages
#'
#' @return List of validated arguments and initialized variables
#' @noRd
validate_and_initialize_args <- function(
    subject_data_list,
    anchor_indices,
    spectral_rank_k,
    task_data_list,
    task_method,
    lambda_blend_value,
    k_gev_dims,
    row_augmentation,
    residualize_condition_anchors,
    omega_weights,
    omega_mode,
    reliability_scores_list,
    scale_omega_trace,
    alpha_laplacian,
    degree_type_laplacian,
    k_conn_pos,
    k_conn_neg,
    k_conn_task_pos,
    k_conn_task_neg,
    similarity_method_task,
    W_task_helper_func,
    n_refine,
    check_redundancy,
    redundancy_threshold,
    residualize_k_conn_proj,
    residualize_k_conn_labels,
    gev_lambda_max,
    gev_epsilon_reg,
    parcel_names,
    verbose
) {
    if (verbose) message_stage(sprintf("Starting task_hatsa run (method: %s)...", task_method), interactive_only = TRUE)

    if (!is.list(subject_data_list) || length(subject_data_list) == 0) {
        stop("subject_data_list must be a non-empty list.")
    }
    N_subjects <- length(subject_data_list)

    first_valid_subj_idx <- which(sapply(subject_data_list, function(x) !is.null(x) && is.matrix(x) && ncol(x) > 0))[1]
    if (is.na(first_valid_subj_idx)) stop("No valid subject data found in subject_data_list.")
    V_p <- ncol(subject_data_list[[first_valid_subj_idx]])
    if (V_p == 0) stop("Subject data has 0 columns (V_p = 0).")

    if (is.null(parcel_names)) {
        parcel_names <- paste0("P", 1:V_p)
    } else if (length(parcel_names) != V_p) {
        stop(sprintf("parcel_names length (%d) does not match inferred V_p (%d).", length(parcel_names), V_p))
    }

    if (!is.numeric(anchor_indices) || any(anchor_indices < 1) || any(anchor_indices > V_p)) {
        stop("anchor_indices must be numeric indices between 1 and V_p.")
    }
    m_parcel_rows <- length(anchor_indices)
    if (m_parcel_rows == 0) warning("No parcel anchor indices provided.")
    if (spectral_rank_k <= 0) stop("spectral_rank_k must be a positive integer.")
    if (spectral_rank_k >= V_p) warning("spectral_rank_k should generally be less than V_p.")

    # First check if this is core_hatsa and handle it separately (allow NULL task_data_list)
    if (task_method == "core_hatsa") {
        # For core_hatsa, we allow NULL task_data_list but disable row_augmentation if it's NULL
        if (row_augmentation && is.null(task_data_list)) {
            if (verbose) message("row_augmentation is TRUE but task_data_list is NULL. Disabling row augmentation for core_hatsa.")
            row_augmentation <- FALSE
        }
    } else {
        # For non-core_hatsa methods, task_data_list is required
        needs_task_data <- task_method != "core_hatsa" || row_augmentation
        if (needs_task_data && (is.null(task_data_list) || length(task_data_list) != N_subjects)) {
            stop("task_data_list is required and must have the same length as subject_data_list when task_method is not 'core_hatsa' or row_augmentation is TRUE.")
        }
    }

    # Handle W_task_helper_func
    current_W_task_helper_func <- W_task_helper_func
    if (task_method != "core_hatsa") {
        if (is.null(current_W_task_helper_func)) {
            warning("W_task_helper_func is NULL. Assuming task_data_list contains activation matrices (C x Vp) for compute_W_task_from_activations.")
            current_W_task_helper_func <- hatsa::compute_W_task_from_activations
        } else if (!is.function(current_W_task_helper_func)) {
            stop("W_task_helper_func must be a function.")
        }
    }

    # Construct the output list with all necessary parameters explicitly included
    out_list <- list(
        N_subjects = N_subjects,
        V_p = V_p,
        m_parcel_rows = m_parcel_rows,
        parcel_names = parcel_names,
        task_method = task_method,
        row_augmentation = row_augmentation,
        lambda_blend_value = lambda_blend_value,
        k_gev_dims = k_gev_dims,
        residualize_condition_anchors = residualize_condition_anchors,
        omega_weights = omega_weights,
        omega_mode = omega_mode,
        reliability_scores_list = reliability_scores_list,
        scale_omega_trace = scale_omega_trace,
        alpha_laplacian = alpha_laplacian,
        degree_type_laplacian = degree_type_laplacian,
        k_conn_pos = k_conn_pos,
        k_conn_neg = k_conn_neg,
        k_conn_task_pos = k_conn_task_pos,
        k_conn_task_neg = k_conn_task_neg,
        similarity_method_task = similarity_method_task,
        W_task_helper_func = current_W_task_helper_func,
        n_refine = n_refine,
        check_redundancy = check_redundancy,
        redundancy_threshold = redundancy_threshold,
        residualize_k_conn_proj = residualize_k_conn_proj,
        residualize_k_conn_labels = residualize_k_conn_labels,
        gev_lambda_max = gev_lambda_max,
        gev_epsilon_reg = gev_epsilon_reg,
        verbose = verbose,
        subject_data_list = subject_data_list,
        anchor_indices = anchor_indices,
        spectral_rank_k = spectral_rank_k,
        task_data_list = task_data_list
    )

    return(out_list)
}

#' Process each subject for basis shaping
#'
#' @param args Validated arguments from validate_and_initialize_args
#' @param subject_data_list List of subject data matrices
#' @param task_data_list List of task data
#'
#' @return List of processing results for each subject
#' @importFrom future.apply future_lapply
#' @noRd
process_subjects <- function(
    args,
    subject_data_list,
    task_data_list
) {
    N_subjects <- args$N_subjects
    verbose <- args$verbose
    
    if (verbose) message_stage(sprintf("Processing %d subjects (Graph construction, Basis shaping) using future_lapply...", N_subjects), interactive_only = TRUE)

    subject_indices <- 1:N_subjects

    # The core logic for a single subject, to be applied in parallel.
    # This function will be serialized and sent to workers if future.globals = FALSE is not specific enough,
    # or if it captures a large environment. Here, it's defined locally so its environment is small.
    # However, process_single_subject and its callees are in the package namespace / helper file scope.
    single_subject_processor_for_future <- function(idx, s_data_list, t_data_list, common_args) {
        # Note: `verbose` from `process_subjects` scope won't be available here if future.globals=FALSE
        # unless explicitly passed or part of common_args and common_args$verbose is used.
        # process_single_subject uses common_args$verbose, which is fine.

        current_subj_data <- s_data_list[[idx]]
        current_task_data <- if (!is.null(t_data_list) && length(t_data_list) >= idx) t_data_list[[idx]] else NULL

        tryCatch({
            result_val <- process_single_subject(
                subject_idx = idx, 
                subj_data_i = current_subj_data, 
                task_data_i = current_task_data, 
                args = common_args
            )
            list(ok = TRUE, subject_idx = idx, result = result_val)
        }, error = function(e) {
            # Optionally log the full error: common_args$verbose might not be reliable here depending on future setup.
            # For now, just capture message.
            warning(sprintf("Error in worker for subject %d: %s", idx, e$message)) # Warning will appear on worker console
            list(ok = FALSE, subject_idx = idx, error_message = e$message, error_object = e)
        })
    }

    all_subject_results_structured <- future.apply::future_lapply(
        X = subject_indices,
        FUN = single_subject_processor_for_future,
        # Pass data and args explicitly to FUN for each iteration
        s_data_list = subject_data_list,
        t_data_list = task_data_list,
        common_args = args, # args now contains all validated & necessary parameters
        future.seed = TRUE,
        future.globals = FALSE # Prevent shipping unnecessary parts of process_subjects environment
                               # Requires FUN to get all its needs from its arguments.
    )

    if (verbose) message_stage(sprintf("Finished parallel processing for %d subjects. Restructuring results...", N_subjects), interactive_only = TRUE)

    # Initialize lists to store results and collect errors
    W_conn_list          <- vector("list", N_subjects)
    L_conn_list          <- vector("list", N_subjects)
    W_task_list          <- vector("list", N_subjects)
    L_task_list          <- vector("list", N_subjects)
    U_original_list      <- vector("list", N_subjects)
    Lambda_original_list <- vector("list", N_subjects)
    U_patch_list         <- vector("list", N_subjects)
    Lambda_patch_list    <- vector("list", N_subjects)
    gev_diagnostics_list <- vector("list", N_subjects)
    qc_metrics_list      <- vector("list", N_subjects)
    
    failed_subjects_info <- list()

    for (i in seq_along(all_subject_results_structured)) {
        res_wrapper <- all_subject_results_structured[[i]]
        # If the result itself is an error (e.g., from future framework), treat as failure
        if (inherits(res_wrapper, "error")) {
             # This case might happen if the error is from the future framework itself, not caught by our tryCatch
            subj_idx <- subject_indices[i] # Best guess for subject index
            err_msg <- conditionMessage(res_wrapper)
            warning(sprintf("Future framework error for subject index %d (approx): %s", subj_idx, err_msg))
            failed_subjects_info[[length(failed_subjects_info) + 1]] <- list(subject_idx = subj_idx, error_message = err_msg)
            # Fill with NULLs for this subject
            W_conn_list[[subj_idx]] <- NULL
            L_conn_list[[subj_idx]] <- NULL
            W_task_list[[subj_idx]] <- NULL
            L_task_list[[subj_idx]] <- NULL
            U_original_list[[subj_idx]] <- NULL
            Lambda_original_list[[subj_idx]] <- NULL
            U_patch_list[[subj_idx]] <- NULL
            Lambda_patch_list[[subj_idx]] <- NULL
            gev_diagnostics_list[[subj_idx]] <- NULL
            qc_metrics_list[[subj_idx]] <- list(rho_redundancy = NA, was_residualized = FALSE) # Default QC
            next
        }


        subj_idx <- res_wrapper$subject_idx # Guarantees correct subject index

        if (res_wrapper$ok) {
            res <- res_wrapper$result
            W_conn_list[[subj_idx]]          <- if (!is.null(res)) res$W_conn else NULL
            L_conn_list[[subj_idx]]          <- if (!is.null(res)) res$L_conn else NULL
            W_task_list[[subj_idx]]          <- if (!is.null(res)) res$W_task else NULL
            L_task_list[[subj_idx]]          <- if (!is.null(res)) res$L_task else NULL
            U_original_list[[subj_idx]]      <- if (!is.null(res)) res$U_original else NULL
            Lambda_original_list[[subj_idx]] <- if (!is.null(res)) res$Lambda_original else NULL
            U_patch_list[[subj_idx]]         <- if (!is.null(res)) res$U_patch else NULL
            Lambda_patch_list[[subj_idx]]    <- if (!is.null(res)) res$Lambda_patch else NULL
            gev_diagnostics_list[[subj_idx]] <- if (!is.null(res)) res$gev_diagnostics else NULL
            qc_metrics_list[[subj_idx]]      <- if (!is.null(res) && !is.null(res$qc_metrics)) {
                                                  res$qc_metrics
                                              } else {
                                                  list(rho_redundancy = NA, was_residualized = FALSE)
                                              }
        } else {
            failed_subjects_info[[length(failed_subjects_info) + 1]] <- list(subject_idx = subj_idx, error_message = res_wrapper$error_message)
            # Fill with NULLs for this subject, qc_metrics with defaults
            W_conn_list[[subj_idx]] <- NULL
            L_conn_list[[subj_idx]] <- NULL
            W_task_list[[subj_idx]] <- NULL
            L_task_list[[subj_idx]] <- NULL
            U_original_list[[subj_idx]] <- NULL
            Lambda_original_list[[subj_idx]] <- NULL
            U_patch_list[[subj_idx]] <- NULL
            Lambda_patch_list[[subj_idx]] <- NULL
            gev_diagnostics_list[[subj_idx]] <- NULL
            qc_metrics_list[[subj_idx]] <- list(rho_redundancy = NA, was_residualized = FALSE)
        }
    }
    
    if (length(failed_subjects_info) > 0) {
        warning_message <- sprintf("Processing failed for %d subject(s):
", length(failed_subjects_info))
        for (fail_info in failed_subjects_info) {
            warning_message <- paste0(warning_message, sprintf("  - Subject %d: %s
", fail_info$subject_idx, fail_info$error_message))
        }
        warning(warning_message, call. = FALSE)
        # Potentially add failed_subjects_info to the returned list or a QC object if more detail is needed upstream.
    }
    
    # The previous loop for restructuring is now replaced by the lapply calls above.
    # The old loop for handling NULL result_i has been incorporated into the lapply logic.

    return(list(
        W_conn_list = W_conn_list,
        L_conn_list = L_conn_list,
        W_task_list = W_task_list,
        L_task_list = L_task_list,
        U_original_list = U_original_list,
        Lambda_original_list = Lambda_original_list,
        U_patch_list = U_patch_list,
        Lambda_patch_list = Lambda_patch_list,
        gev_diagnostics_list = gev_diagnostics_list,
        qc_metrics_list = qc_metrics_list
    ))
}

#' Process a single subject
#'
#' @param subject_idx Subject index.
#' @param subj_data_i Subject data matrix
#' @param task_data_i Task data for the subject
#' @param args Validated arguments
#'
#' @return List of processing results for this subject
#' @noRd
process_single_subject <- function(subject_idx, subj_data_i, task_data_i, args) {
    # Initialize results
    result <- list(
        W_conn = NULL,
        L_conn = NULL,
        W_task = NULL,
        L_task = NULL,
        U_original = NULL, 
        Lambda_original = NULL,
        U_patch = NULL,
        Lambda_patch = NULL,
        gev_diagnostics = NULL,
        qc_metrics = list(rho_redundancy = NA, was_residualized = FALSE)
    )
    
    # Extract arguments from args list for local use if needed, or pass args directly
    V_p <- args$V_p
    verbose <- args$verbose
    parcel_names <- args$parcel_names
    k_conn_pos <- args$k_conn_pos
    k_conn_neg <- args$k_conn_neg
    spectral_rank_k <- args$spectral_rank_k # Used in shape_basis
    task_method <- args$task_method # Used in shape_basis and here
    
    # Basic validation of subject data
    if (is.null(subj_data_i) || !is.matrix(subj_data_i) || ncol(subj_data_i) != V_p) {
        warning(sprintf("Skipping subject %d due to invalid data (expected %d columns).", subject_idx, V_p))
        return(result)
    }
    
    # --- Compute W_conn and L_conn ---
    W_conn_i <- tryCatch({
        compute_subject_connectivity_graph_sparse(subj_data_i, parcel_names, k_conn_pos, k_conn_neg)
    }, error = function(e) {
        warning(sprintf("Error computing W_conn for subject %d: %s. Skipping subject.", subject_idx, e$message)); NULL
    })
    if (is.null(W_conn_i)) return(result) # Return initialized result list (all NULLs essentially)
    result$W_conn <- W_conn_i
    
    L_conn_i <- tryCatch({
        compute_graph_laplacian_sparse(W_conn_i, alpha = args$alpha_laplacian, degree_type = args$degree_type_laplacian)
    }, error = function(e) {
        warning(sprintf("Error computing L_conn for subject %d: %s. Skipping subject.", subject_idx, e$message)); NULL
    })
    if (is.null(L_conn_i)) return(result)
    result$L_conn <- L_conn_i
    
    # --- Compute W_task and L_task if needed ---
    W_task_i_for_basis <- NULL # Initialize to NULL, will be set if task_method requires it
    L_task_i_for_basis <- NULL # Initialize to NULL

    if (task_method != "core_hatsa") {
        task_result_data <- compute_task_matrices(subject_idx, task_data_i, args, W_conn_i, L_conn_i)
        result$W_task <- task_result_data$W_task_i # Store final W_task (raw or residualized)
        result$L_task <- task_result_data$L_task_i # Store final L_task
        result$qc_metrics <- task_result_data$qc_metrics # Update QC from task processing

        # These are used by shape_basis
        W_task_i_for_basis <- result$W_task 
        L_task_i_for_basis <- result$L_task
    }
    
    # --- Perform basis shaping based on task_method ---
    # Pass W_conn_i and potentially W_task_i_for_basis to shape_basis for lambda_blend
    basis_result <- shape_basis(subject_idx, L_conn_i, L_task_i_for_basis, args, W_conn_i, W_task_i_for_basis)
    result$U_original <- basis_result$U_original
    result$Lambda_original <- basis_result$Lambda_original
    result$U_patch <- basis_result$U_patch
    result$Lambda_patch <- basis_result$Lambda_patch
    result$gev_diagnostics <- basis_result$gev_diagnostics
    
    return(result)
}

#' Compute task matrices (W_task and L_task)
#'
#' @param i Subject index. Renamed to subject_idx for clarity.
#' @param task_data_i Task data for subject i
#' @param args Validated arguments
#' @param W_conn_i Connectivity matrix for subject i
#' @param L_conn_i Laplacian of W_conn_i, needed for residualization projection
#'
#' @return List with W_task_i, L_task_i, and qc_metrics
#' @noRd
compute_task_matrices <- function(subject_idx, task_data_i, args, W_conn_i, L_conn_i) {
    result <- list(
        W_task_i = NULL,
        L_task_i = NULL,
        qc_metrics = list(rho_redundancy = NA, was_residualized = FALSE)
    )
    
    # Extract arguments
    verbose <- args$verbose
    parcel_names <- args$parcel_names
    k_conn_task_pos <- args$k_conn_task_pos
    k_conn_task_neg <- args$k_conn_task_neg
    similarity_method_task <- args$similarity_method_task
    W_task_helper_func <- args$W_task_helper_func
    check_redundancy <- args$check_redundancy
    redundancy_threshold <- args$redundancy_threshold
    alpha_laplacian <- args$alpha_laplacian
    degree_type_laplacian <- args$degree_type_laplacian
    
    # Basic validation for task data
    if (is.null(task_data_i) || !is.matrix(task_data_i)) {
        warning(sprintf("Task data for subject %d is missing or not a matrix. Cannot compute W_task.", subject_idx))
        return(result)
    }
    
    # Compute W_task_raw
    W_task_i_raw <- tryCatch({
        W_task_helper_func(task_data_i, parcel_names = parcel_names,
                          k_conn_task_pos = k_conn_task_pos, k_conn_task_neg = k_conn_task_neg,
                          similarity_method = similarity_method_task)
    }, error = function(e) {
        warning(sprintf("Error computing W_task_raw for subject %d: %s.", subject_idx, e$message)); NULL
    })
    
    if (is.null(W_task_i_raw)) return(result)
    
    W_task_i <- W_task_i_raw # Start with raw
    
    # Redundancy Check
    if (check_redundancy) {
        rho_i <- tryCatch({
            compute_graph_correlation(W_conn_i, W_task_i_raw)
        }, error = function(e) {
            warning(sprintf("Error computing graph correlation for subject %d: %s", subject_idx, e$message)); NA
        })
        
        result$qc_metrics$rho_redundancy <- rho_i
        
        if (!is.na(rho_i) && rho_i >= redundancy_threshold) {
            if (verbose) message(sprintf("  - Subject %d: W_task/W_conn redundancy rho=%.3f >= %.3f. Residualizing W_task.", subject_idx, rho_i, redundancy_threshold))
            W_task_i <- tryCatch({
                residualize_graph_on_subspace(
                    W_graph_to_residualize = W_task_i_raw,
                    L_graph_for_projection = L_conn_i,
                    k_eigenvectors_to_remove = args$residualize_k_conn_proj,
                    k_nn_resparsify = args$residualize_k_conn_labels
                )
            }, error = function(e) {
                warning(sprintf("Error residualizing W_task for subject %d: %s. Using W_task_raw.", subject_idx, e$message)); W_task_i_raw
            })
            result$qc_metrics$was_residualized <- TRUE
        }
    }
    
    # Compute L_task from the final W_task_i (raw or residualized)
    L_task_i <- tryCatch({
        compute_graph_laplacian_sparse(W_task_i, alpha = alpha_laplacian, degree_type = degree_type_laplacian)
    }, error = function(e) {
        warning(sprintf("Error computing L_task for subject %d: %s.", subject_idx, e$message)); NULL
    })
    
    result$W_task_i <- W_task_i
    result$L_task_i <- L_task_i
    
    return(result)
}

#' Shape basis based on task method
#'
#' @param subject_idx Subject index.
#' @param L_conn_i Laplacian matrix for connectivity
#' @param L_task_i Laplacian matrix for task
#' @param args Validated arguments
#' @param W_conn_i Connectivity matrix for subject i
#' @param W_task_i Task matrix for subject i
#'
#' @return List with basis matrices and eigenvalues
#' @noRd
shape_basis <- function(subject_idx, L_conn_i, L_task_i, args, W_conn_i, W_task_i) {
    result <- list(
        U_original = NULL,
        Lambda_original = NULL,
        U_patch = NULL,
        Lambda_patch = NULL,
        gev_diagnostics = NULL
    )
    
    task_method <- args$task_method
    spectral_rank_k <- args$spectral_rank_k
    
    if (task_method == "core_hatsa") {
        sketch <- tryCatch({
            compute_spectral_sketch_sparse(L_conn_i, spectral_rank_k, eigenvalue_tol = 1e-8)
        }, error = function(e) {
            warning(sprintf("Error computing spectral sketch (core) for subject %d: %s.", subject_idx, e$message)); NULL
        })
        if (is.null(sketch)) return(result)
        result$U_original <- sketch$vectors
        result$Lambda_original <- sketch$values
        
    } else if (task_method == "lambda_blend") {
        if (is.null(W_conn_i)) { # W_conn_i is essential
            warning(sprintf("W_conn_i is NULL for subject %d (task_method=%s). Cannot perform blend. Skipping basis for this subject.", subject_idx, task_method))
            return(result) # Return empty result
        }
        if (is.null(W_task_i) && args$lambda_blend_value > 0) { # W_task_i is essential if lambda > 0
            warning(sprintf("W_task_i is NULL for subject %d (task_method=%s) but lambda_blend_value > 0. Cannot perform blend. Using core W_conn for basis.", subject_idx, task_method))
            # Fallback to core L_conn logic if W_task is missing but was expected
            if (is.null(L_conn_i)) { # Should not happen if W_conn_i was present
                 warning(sprintf("L_conn_i also NULL for subject %d. Skipping basis.", subject_idx))
                 return(result)
            }
            sketch <- tryCatch({
                compute_spectral_sketch_sparse(L_conn_i, spectral_rank_k, eigenvalue_tol = 1e-8)
            }, error = function(e) {
                warning(sprintf("Error computing spectral sketch (lambda_blend fallback to core L_conn) for subject %d: %s.", subject_idx, e$message)); NULL
            })
        } else {
            # Blend W matrices first
            if (args$lambda_blend_value == 0) { # Effectively core_hatsa
                W_hybrid_i <- W_conn_i
            } else if (args$lambda_blend_value == 1 && !is.null(W_task_i)) { # Effectively task-only graph
                W_hybrid_i <- W_task_i
            } else if (is.null(W_task_i)) { # lambda > 0 but W_task is null, already warned, use W_conn
                 W_hybrid_i <- W_conn_i
            } else { # Actual blend
                W_hybrid_i <- (1 - args$lambda_blend_value) * W_conn_i + args$lambda_blend_value * W_task_i
            }
            
            # Compute Laplacian from the blended W_hybrid_i
            L_hybrid_i <- tryCatch({
                compute_graph_laplacian_sparse(W_hybrid_i, alpha = args$alpha_laplacian, degree_type = args$degree_type_laplacian)
            }, error = function(e) {
                warning(sprintf("Error computing L_hybrid_i from blended W for subject %d: %s. Skipping basis.", subject_idx, e$message)); NULL
            })

            if(is.null(L_hybrid_i)) return(result) # Stop if Laplacian computation failed

            sketch <- tryCatch({
                compute_spectral_sketch_sparse(L_hybrid_i, spectral_rank_k, eigenvalue_tol = 1e-8)
            }, error = function(e) {
                warning(sprintf("Error computing spectral sketch (lambda_blend) for subject %d: %s.", subject_idx, e$message)); NULL
            })
        }
        if (is.null(sketch)) return(result) # U_original and Lambda_original will be NULL
        result$U_original <- sketch$vectors
        result$Lambda_original <- sketch$values
        
    } else if (task_method == "gev_patch") {
        # Core sketch
        core_sketch <- tryCatch({
            compute_spectral_sketch_sparse(L_conn_i, spectral_rank_k, eigenvalue_tol = 1e-8)
        }, error = function(e) {
            warning(sprintf("Error computing spectral sketch (core for GEV) for subject %d: %s.", subject_idx, e$message)); NULL
        })
        if (is.null(core_sketch)) return(result)
        result$U_original <- core_sketch$vectors
        result$Lambda_original <- core_sketch$values
        
        # GEV Patch
        if (is.null(L_task_i)) {
            warning(sprintf("L_task is NULL for subject %d (task_method=%s). Cannot compute GEV patch.", subject_idx, task_method))
        } else {
            gev_results <- tryCatch({
                solve_gev_laplacian_primme(
                    L_task_i, L_conn_i,
                    k_request = args$k_gev_dims * 2, # Request more, filter later
                    lambda_max_thresh = args$gev_lambda_max,
                    epsilon_reg_B = args$gev_epsilon_reg
                )
            }, error = function(e) {
                warning(sprintf("Error solving GEV for subject %d: %s.", subject_idx, e$message)); NULL
            })
            
            if (!is.null(gev_results)) {
                result$U_patch <- gev_results$vectors
                result$Lambda_patch <- gev_results$values
                result$gev_diagnostics <- compute_gev_spectrum_diagnostics(gev_results$values, args$gev_lambda_max)
            }
        }
    }
    
    return(result)
}

#' Perform anchor augmentation
#'
#' @param args Validated arguments
#' @param processing_results Results from process_subjects
#' @param anchor_indices Indices of anchor parcels
#'
#' @return List with augmentation results
#' @noRd
perform_anchor_augmentation <- function(args, processing_results, anchor_indices) {
    verbose <- args$verbose
    N_subjects <- args$N_subjects
    row_augmentation <- args$row_augmentation
    task_data_list <- args$task_data_list
    m_parcel_rows <- args$m_parcel_rows
    spectral_rank_k <- args$spectral_rank_k
    
    U_original_list <- processing_results$U_original_list
    
    if (verbose) message_stage("Performing Anchor Augmentation (if enabled)...", interactive_only = TRUE)
    
    valid_subject_indices_for_gpa <- which(!sapply(U_original_list, is.null))
    if (length(valid_subject_indices_for_gpa) < 2) {
        stop(sprintf("Need at least 2 subjects with valid computed bases for GPA, found only %d.", length(valid_subject_indices_for_gpa)))
    }
    
    # Initialize storage
    A_originals_list_for_gpa <- vector("list", N_subjects)
    m_task_rows_effective <- 0
    condition_labels_for_anchors <- NULL
    condition_labels_set <- FALSE
    
    for (i in valid_subject_indices_for_gpa) {
        augmentation_result <- augment_anchors_for_subject(
            i, 
            U_original_list[[i]], 
            row_augmentation,
            if (!is.null(task_data_list) && length(task_data_list) >= i) task_data_list[[i]] else NULL,
            anchor_indices,
            args
        )
        
        A_originals_list_for_gpa[[i]] <- augmentation_result$A_augmented
        
        if (m_task_rows_effective == 0 && augmentation_result$C_subj > 0) {
            m_task_rows_effective <- augmentation_result$C_subj
        }
        
        if (!condition_labels_set && !is.null(augmentation_result$condition_labels)) {
            condition_labels_for_anchors <- augmentation_result$condition_labels
            condition_labels_set <- TRUE
        }
    }
    
    # Check consistency of augmented matrices before GPA
    final_m_total_rows <- m_parcel_rows + m_task_rows_effective
    valid_A_indices_for_gpa_final <- sapply(A_originals_list_for_gpa, function(A) {
        !is.null(A) && is.matrix(A) && nrow(A) == final_m_total_rows && ncol(A) == spectral_rank_k
    })
    num_valid_for_gpa <- sum(valid_A_indices_for_gpa_final)
    
    if (num_valid_for_gpa < 2) {
        stop(sprintf("Need at least 2 subjects with valid anchor matrices for GPA (found %d with %d rows). Check basis computation and anchor augmentation steps.", num_valid_for_gpa, final_m_total_rows))
    }
    
    if (verbose) message(sprintf("Proceeding to GPA with %d subjects having valid anchor matrices (%d parcel + %d task rows).", num_valid_for_gpa, m_parcel_rows, m_task_rows_effective))
    
    return(list(
        A_originals_list_for_gpa = A_originals_list_for_gpa,
        valid_A_indices_for_gpa_final = valid_A_indices_for_gpa_final,
        m_task_rows_effective = m_task_rows_effective,
        condition_labels_for_anchors = condition_labels_for_anchors
    ))
}

#' Augment anchors for a single subject
#'
#' @param i Subject index
#' @param U_basis_i Basis matrix for subject i
#' @param row_augmentation Whether to augment with task rows
#' @param task_data_i Task data for subject i
#' @param anchor_indices Indices of anchor parcels
#' @param args Validated arguments
#'
#' @return List with augmented anchor matrix and related information
#' @noRd
augment_anchors_for_subject <- function(i, U_basis_i, row_augmentation, task_data_i, anchor_indices, args) {
    result <- list(
        A_augmented = NULL,
        C_subj = 0,
        condition_labels = NULL
    )
    
    verbose <- args$verbose
    V_p <- args$V_p
    
    # Extract parcel anchors
    if (any(anchor_indices > nrow(U_basis_i))) {
        stop(sprintf("Subject %d: anchor_indices out of bounds for basis matrix.", i))
    }
    
    A_parc_i <- U_basis_i[anchor_indices, , drop = FALSE]
    result$A_augmented <- A_parc_i # Default if no augmentation
    
    if (row_augmentation && !is.null(task_data_i)) {
        # Check if task_data_i is suitable (e.g., matrix C x Vp)
        is_suitable_task_data <- is.matrix(task_data_i) && nrow(task_data_i) > 0 && ncol(task_data_i) == V_p
        
        if (is_suitable_task_data) {
            result$C_subj <- nrow(task_data_i)
            
            # Store condition labels if available
            if (!is.null(rownames(task_data_i))) {
                result$condition_labels <- rownames(task_data_i)
            }
            
            # Project features
            Act_i <- t(task_data_i) # Transpose for projection
            Z_i_projected <- tryCatch({
                project_features_to_spectral_space(feature_matrix = Act_i, U_basis = U_basis_i)
            }, error = function(e) {
                warning(sprintf("Error projecting task features for subject %d: %s.", i, e$message)); NULL
            })
            
            if (!is.null(Z_i_projected)) {
                # Output is k x C, need C x k for rbind
                Z_i <- t(Z_i_projected)
                
                # Optional Residualization
                if (args$residualize_condition_anchors) {
                    if (verbose) message(sprintf("  - Subject %d: Residualizing condition anchors.", i))
                    Z_i <- tryCatch({
                        residualize_matrix_on_subspace(matrix_to_residualize = Z_i, subspace_basis_matrix = A_parc_i)
                    }, error = function(e) {
                        warning(sprintf("Error residualizing task anchors for subject %d: %s.", i, e$message)); Z_i
                    })
                }
                
                # Build augmented matrix
                result$A_augmented <- build_augmented_anchor_matrix(A_parc_i, Z_i)
            } else {
                if (verbose) message(sprintf("  - Subject %d: Task feature projection failed. Using parcel anchors only.", i))
            }
        } else {
            if (verbose) message(sprintf("  - Subject %d: Task data not suitable for row augmentation (expected C x Vp matrix). Using parcel anchors only.", i))
        }
    }
    
    return(result)
}

#' Prepare and run GPA refinement
#'
#' @param args Validated arguments
#' @param augmentation_results Results from perform_anchor_augmentation
#' @param reliability_scores_list List of reliability scores for adaptive weighting
#'
#' @return Results from GPA refinement
#' @noRd
prepare_and_run_gpa <- function(args, augmentation_results, reliability_scores_list) {
    verbose <- args$verbose
    N_subjects <- args$N_subjects
    original_omega_mode <- args$omega_mode # Store the original mode
    n_refine <- args$n_refine
    m_parcel_rows <- args$m_parcel_rows
    spectral_rank_k <- args$spectral_rank_k
    scale_omega_trace <- args$scale_omega_trace
    
    m_task_rows_effective <- augmentation_results$m_task_rows_effective
    A_originals_list_for_gpa <- augmentation_results$A_originals_list_for_gpa
    valid_A_indices_for_gpa_final <- augmentation_results$valid_A_indices_for_gpa_final
    
    if (verbose) message_stage(sprintf("Performing GPA (%d iterations)...", n_refine), interactive_only = TRUE)
    
    # Prepare reliability scores list for GPA
    reliability_scores_for_gpa <- prepare_reliability_scores(
        N_subjects, 
        m_task_rows_effective, 
        reliability_scores_list, 
        valid_A_indices_for_gpa_final, 
        original_omega_mode, # Pass original mode here
        verbose
    )
    
    # Determine effective omega_mode for the solver call
    effective_omega_mode <- original_omega_mode
    if (original_omega_mode == "adaptive" && is.null(reliability_scores_for_gpa)) {
        if (verbose && m_task_rows_effective > 0) { # Only message if adaptive was relevant
            message("  - GPA: omega_mode was 'adaptive' but no valid reliability scores were available. Defaulting to 'fixed' omega_mode for GPA solver.")
        }
        effective_omega_mode <- "fixed"
    }

    # Run GPA refinement
    gpa_results <- tryCatch({
        perform_gpa_refinement(
            A_originals_list = A_originals_list_for_gpa,
            n_refine = n_refine,
            k = spectral_rank_k,
            m_parcel_rows = m_parcel_rows,
            m_task_rows = m_task_rows_effective,
            omega_mode = effective_omega_mode, # Use the potentially adjusted mode
            fixed_omega_weights = args$omega_weights,
            reliability_scores_list = reliability_scores_for_gpa, # This might be NULL
            scale_omega_trace = scale_omega_trace
        )
    }, error = function(e) {
        stop(sprintf("Error during GPA refinement: %s", e$message)); NULL
    })
    
    if (is.null(gpa_results)) stop("GPA refinement failed.")
    
    return(gpa_results)
}

#' Prepare reliability scores for GPA
#'
#' @param N_subjects Number of subjects
#' @param m_task_rows_effective Effective number of task rows
#' @param reliability_scores_list List of reliability scores
#' @param valid_A_indices Valid indices for GPA
#' @param omega_mode Mode for omega calculation
#' @param verbose Whether to print messages
#'
#' @return Prepared reliability scores list
#' @noRd
prepare_reliability_scores <- function(N_subjects, m_task_rows_effective, reliability_scores_list, valid_A_indices, omega_mode, verbose) {
    if (omega_mode == "adaptive" && m_task_rows_effective > 0 && !is.null(reliability_scores_list)) {
        reliability_scores_for_gpa <- vector("list", N_subjects)
        valid_reliability_provided <- FALSE
        
        for (i in 1:N_subjects) {
            if (valid_A_indices[i] && length(reliability_scores_list) >= i && !is.null(reliability_scores_list[[i]])) {
                if (length(reliability_scores_list[[i]]) == m_task_rows_effective) {
                    reliability_scores_for_gpa[[i]] <- reliability_scores_list[[i]]
                    valid_reliability_provided <- TRUE
                } else {
                    warning(sprintf("Subject %d: reliability_scores length (%d) doesn't match effective task rows (%d). Ignoring for this subject.", 
                                   i, length(reliability_scores_list[[i]]), m_task_rows_effective))
                }
            }
        }
        
        if (!valid_reliability_provided) {
            # Warning is now more informational, as prepare_and_run_gpa will handle the mode switch.
            if (verbose) {
                message("  - GPA Setup: omega_mode='adaptive' but no valid reliability_scores_list provided or matched subjects. GPA will use fixed weights.")
            }
            reliability_scores_for_gpa <- NULL
        }
        
        return(reliability_scores_for_gpa)
    }
    
    return(NULL)
}

#' Perform patch alignment if GEV method is used
#'
#' @param args Validated arguments
#' @param processing_results Results from process_subjects
#' @param gpa_results Results from GPA
#'
#' @return Patch alignment results or NULL
#' @noRd
perform_patch_alignment <- function(args, processing_results, gpa_results) {
    task_method <- args$task_method
    verbose <- args$verbose
    N_subjects <- args$N_subjects
    m_parcel_rows <- args$m_parcel_rows
    n_refine <- args$n_refine
    anchor_indices <- args$anchor_indices
    
    U_patch_list <- processing_results$U_patch_list
    
    R_patch_list <- NULL
    
    if (task_method == "gev_patch") {
        if (verbose) message_stage("Performing GEV Patch Alignment...", interactive_only = TRUE)
        
        # Check if any patches were actually computed
        valid_patch_indices <- which(!sapply(U_patch_list, function(p) is.null(p) || ncol(p) == 0))
        
        if (length(valid_patch_indices) >= 2) {
            patch_results <- align_gev_patches(valid_patch_indices, U_patch_list, anchor_indices, m_parcel_rows, n_refine, verbose)
            R_patch_list <- patch_results$R_patch_list
        } else {
            if (verbose) message("Not enough valid GEV patches computed (need >= 2) for patch alignment.")
        }
    }
    
    return(list(R_patch_list = R_patch_list))
}

#' Align GEV patches
#'
#' @param valid_patch_indices Indices of valid patches
#' @param U_patch_list List of patch matrices
#' @param anchor_indices Indices of anchor parcels
#' @param m_parcel_rows Number of parcel rows
#' @param n_refine Number of GPA refinement iterations
#' @param verbose Whether to print messages
#'
#' @return List with patch rotation matrices
#' @noRd
align_gev_patches <- function(valid_patch_indices, U_patch_list, anchor_indices, m_parcel_rows, n_refine, verbose) {
    k_patch_dims <- ncol(U_patch_list[[valid_patch_indices[1]]])
    N_subjects <- length(U_patch_list)
    
    # Create anchor list for patches
    A_patch_originals_list <- vector("list", N_subjects)
    valid_indices_for_patch_gpa <- logical(N_subjects)
    
    for (i in valid_patch_indices) {
        U_patch_i <- U_patch_list[[i]]
        
        if (ncol(U_patch_i) == k_patch_dims) {
            if (any(anchor_indices > nrow(U_patch_i))) {
                warning(sprintf("Cannot extract anchors for GEV patch for subject %d (indices out of bounds). Skipping patch alignment for this subject.", i))
                next
            }
            
            A_patch_originals_list[[i]] <- U_patch_i[anchor_indices, , drop = FALSE]
            valid_indices_for_patch_gpa[i] <- TRUE
        } else {
            warning(sprintf("Subject %d GEV patch dimension (%d) differs from expected (%d). Skipping patch alignment for this subject.", 
                           i, ncol(U_patch_i), k_patch_dims))
        }
    }
    
    if (sum(valid_indices_for_patch_gpa) >= 2) {
        gpa_patch_results <- tryCatch({
            # Use unweighted GPA for patches by default (m_task_rows = 0)
            perform_gpa_refinement(
                A_originals_list = A_patch_originals_list,
                n_refine = n_refine,
                k = k_patch_dims,
                m_parcel_rows = m_parcel_rows,
                m_task_rows = 0 # No task augmentation for patches
            )
        }, error = function(e) {
            warning(sprintf("Error during GEV Patch GPA refinement: %s.", e$message)); NULL
        })
        
        if (!is.null(gpa_patch_results)) {
            return(list(R_patch_list = gpa_patch_results$R_final_list))
        } else {
            warning("GEV Patch GPA failed. Patch rotations not computed.")
            return(list(R_patch_list = NULL))
        }
    } else {
        warning("Fewer than 2 subjects have valid, consistent GEV patches for patch alignment.")
        return(list(R_patch_list = NULL))
    }
}

#' Construct aligned U matrices and final output
#'
#' @param args Validated arguments
#' @param processing_results Results from process_subjects
#' @param augmentation_results Results from perform_anchor_augmentation
#' @param gpa_results Results from GPA
#' @param patch_results Results from patch alignment
#'
#' @return Final task_hatsa output
#' @noRd
construct_output <- function(args, processing_results, augmentation_results, gpa_results, patch_results) {
    verbose <- args$verbose
    N_subjects <- args$N_subjects
    
    U_original_list <- processing_results$U_original_list
    Lambda_original_list <- processing_results$Lambda_original_list
    U_patch_list <- processing_results$U_patch_list
    Lambda_patch_list <- processing_results$Lambda_patch_list
    gev_diagnostics_list <- processing_results$gev_diagnostics_list
    qc_metrics_list <- processing_results$qc_metrics_list
    
    # Extract task-related outputs if available
    W_task_list <- processing_results$W_task_list
    L_task_list <- processing_results$L_task_list
    W_hybrid_list <- processing_results$W_hybrid_list  
    L_hybrid_list <- processing_results$L_hybrid_list
    U_task_list <- processing_results$U_task_list
    Lambda_task_list <- processing_results$Lambda_task_list
    
    R_final_list <- gpa_results$R_final_list
    T_anchor_final <- gpa_results$T_anchor_final
    R_patch_list <- patch_results$R_patch_list
    
    if (verbose) message_stage("Constructing output object...", interactive_only = TRUE)
    
    # Compute aligned U matrices (core basis)
    U_aligned_list <- compute_aligned_U_matrices(N_subjects, U_original_list, R_final_list)
    
    # Build the final result list
    result_list <- list(
        method = "task_hatsa", # This will be set/overridden by the constructor based on parameters
        parameters = create_parameters_list(args),
        # Core alignment results
        R_final_list = R_final_list,
        T_anchor_final = T_anchor_final,
        U_original_list = U_original_list,
        Lambda_original_list = Lambda_original_list,
        U_aligned_list = U_aligned_list,
        
        # Task-specific outputs
        W_task_list = W_task_list,
        L_task_list = L_task_list,
        W_hybrid_list = W_hybrid_list,
        L_hybrid_list = L_hybrid_list,
        U_task_list = U_task_list,
        Lambda_task_list = Lambda_task_list,
        
        # QC and augmentation info
        qc_metrics = qc_metrics_list,
        anchor_augmentation_info = list(
            m_parcel_rows = args$m_parcel_rows,
            m_task_rows_effective = augmentation_results$m_task_rows_effective,
            condition_labels = augmentation_results$condition_labels_for_anchors,
            was_residualized = args$residualize_condition_anchors,
            omega_mode_used = args$omega_mode,
            omega_weights_params = args$omega_weights,
            trace_scaled = args$scale_omega_trace
        ),
        gev_patch_data = if (args$task_method == "gev_patch") list(
            U_patch_list = U_patch_list,
            Lambda_patch_list = Lambda_patch_list,
            R_patch_list = R_patch_list,
            diagnostics = gev_diagnostics_list
        ) else NULL
    )
    
    if (verbose) message_stage("task_hatsa run completed.", interactive_only = TRUE)
    # Call the constructor to create the S3 object
    return(task_hatsa_projector(result_list))
}

#' Compute aligned U matrices
#'
#' @param N_subjects Number of subjects
#' @param U_original_list List of original U matrices
#' @param R_final_list List of rotation matrices
#'
#' @return List of aligned U matrices
#' @noRd
compute_aligned_U_matrices <- function(N_subjects, U_original_list, R_final_list) {
    U_aligned_list <- vector("list", N_subjects)
    
    for (i in 1:N_subjects) {
        if (!is.null(U_original_list[[i]]) && !is.null(R_final_list[[i]])) {
            U_aligned_list[[i]] <- tryCatch({
                U_original_list[[i]] %*% R_final_list[[i]]
            }, error = function(e) {
                warning(sprintf("Error aligning U for subject %d: %s", i, e$message)); NULL
            })
        }
    }
    
    return(U_aligned_list)
}

#' Create parameters list for output
#'
#' @param args Validated arguments
#'
#' @return List of parameters for output
#' @noRd
create_parameters_list <- function(args) {
    return(list(
        k = args$spectral_rank_k,
        spectral_rank_k = args$spectral_rank_k,
        task_method = args$task_method,
        lambda_blend_value = args$lambda_blend_value,
        k_gev_dims = args$k_gev_dims,
        row_augmentation = args$row_augmentation,
        residualize_condition_anchors = args$residualize_condition_anchors,
        omega_mode = args$omega_mode,
        fixed_omega_weights = args$omega_weights,
        scale_omega_trace = args$scale_omega_trace,
        alpha_laplacian = args$alpha_laplacian,
        degree_type_laplacian = args$degree_type_laplacian,
        k_conn_pos = args$k_conn_pos,
        k_conn_neg = args$k_conn_neg,
        k_conn_task_pos = args$k_conn_task_pos,
        k_conn_task_neg = args$k_conn_task_neg,
        similarity_method_task = args$similarity_method_task,
        n_refine = args$n_refine,
        check_redundancy = args$check_redundancy,
        redundancy_threshold = args$redundancy_threshold,
        residualize_k_conn_proj = args$residualize_k_conn_proj,
        residualize_k_conn_labels = args$residualize_k_conn_labels,
        gev_lambda_max = args$gev_lambda_max,
        gev_epsilon_reg = args$gev_epsilon_reg,
        V_p = args$V_p,
        N_subjects = args$N_subjects,
        anchor_indices = args$anchor_indices
    ))
}



