# R/task_hatsa_main.R

#' Run Task-Informed HATSA (task_hatsa)
#'
#' Performs Hyperalignment via Task-Informed Shared Analysis (task_hatsa) on a list of subject data,
#' incorporating task information during basis shaping and/or alignment refinement.
#'
#' @param subject_data_list List where each element is a numeric matrix (`T_i x V_p`)
#'   of time-series data for one subject.
#' @param anchor_indices Integer vector, indices of the canonical anchor parcels (`1` to `V_p`).
#' @param spectral_rank_k Integer, the desired dimensionality of the primary spectral sketch.
#' @param task_data_list List (parallel to `subject_data_list`). Each element provides
#'   task-related data for a subject. Structure depends on `task_method` and subsequent
#'   processing (e.g., a `C x V_p` matrix of activations/betas for `compute_W_task_from_activations`
#'   or anchor augmentation). Required if `task_method != "core_hatsa"` or `row_augmentation=TRUE`.
#'   Assumed cross-validated if necessary before input.
#' @param task_method Character string: `"core_hatsa"`, `"lambda_blend"`, or `"gev_patch"`. Default: `"lambda_blend"`.
#' @param lambda_blend_value Numeric `lambda` in `[0,1]`. Weight for `L_task` in blend. Default 0.15.
#' @param k_gev_dims Integer, requested dimension for GEV patches. Default 10. Used if `task_method == "gev_patch"`.
#' @param row_augmentation Logical. If `TRUE`, add projected task features to anchor matrices
#'   for GPA refinement. Requires suitable `task_data_list`. Default `TRUE` if suitable data provided.
#' @param residualize_condition_anchors Logical. If `TRUE` and `row_augmentation` is `TRUE`,
#'   residualize projected task anchors against parcel anchors. Default `FALSE`.
#' @param omega_weights List specifying fixed weights for weighted Procrustes (e.g.,
#'   `list(parcel = 1.0, condition = 0.5)`). Used if `row_augmentation=TRUE` and `omega_mode == "fixed"`.
#'   Defaults handled by `solve_procrustes_rotation_weighted`.
#' @param omega_mode Character string: `"fixed"` or `"adaptive"`. Controls weighting in GPA. Default `"fixed"`.
#' @param reliability_scores_list List (parallel to `subject_data_list`), each element a numeric
#'   vector of reliability scores (e.g., R^2) for task data (length `C`). Used if `omega_mode == "adaptive"`.
#' @param scale_omega_trace Logical. Whether to rescale weights in weighted GPA so trace equals total anchors. Default `TRUE`.
#' @param alpha_laplacian Numeric, laziness parameter for graph Laplacians (`L = I - alpha D^{-1} W`). Default 0.93.
#' @param degree_type_laplacian Character string (`"abs"`, `"positive"`, `"signed"`). Type of degree calculation for Laplacian. Default `"abs"`.
#' @param k_conn_pos Integer >= 0. k-NN sparsification for positive edges in `W_conn`.
#' @param k_conn_neg Integer >= 0. k-NN sparsification for negative edges in `W_conn`.
#' @param k_conn_task_pos Integer >= 0. k-NN sparsification for positive edges in `W_task`.
#' @param k_conn_task_neg Integer >= 0. k-NN sparsification for negative edges in `W_task`.
#' @param similarity_method_task Character string or function. Method to compute similarity for `W_task`
#'   (e.g., "pearson", "spearman"). Default "pearson".
#' @param W_task_helper_func Function. The specific function to compute `W_task` (e.g.,
#'   `compute_W_task_from_activations`, `compute_W_task_from_encoding`). If `NULL`, attempts
#'   to infer based on `task_data_list` structure (currently assumes activations `C x Vp`). Default `NULL`.
#' @param n_refine Integer >= 0. Number of GPA refinement iterations.
#' @param check_redundancy Logical. If `TRUE`, check correlation between `W_conn` and `W_task`. Default `TRUE`.
#' @param redundancy_threshold Numeric. Spearman rho threshold for triggering `W_task` residualization. Default 0.45.
#' @param residualize_k_conn_proj Integer. Number of `L_conn` eigenvectors to project `W_task` out of. Default 64.
#' @param residualize_k_conn_labels Integer. k-NN value for re-sparsifying `W_task_res` after residualization. Default 10.
#' @param gev_lambda_max Numeric. Max GEV eigenvalue `lambda` to retain for patches. Default 0.8.
#' @param gev_epsilon_reg Numeric. Small regularization for `L_conn` in GEV. Default 1e-6.
#' @param parcel_names Optional character vector of parcel names. If `NULL`, names like "P1", "P2"... are generated.
#' @param verbose Logical. Print progress messages? Default `TRUE`.
#'
#' @return A list representing the `task_hatsa_projector` object (structure TBD, needs final class definition).
#'         Contains aligned bases, rotations, parameters, QC metrics, etc.
#'
#' @export
#' @importFrom Matrix Matrix Diagonal t crossprod sparseMatrix drop0 forceSymmetric is
#' @importFrom stats cor sd median runif
#' @importFrom utils head
#' @importFrom methods as
run_task_hatsa <- function(
    subject_data_list,
    anchor_indices,
    spectral_rank_k,
    task_data_list = NULL,
    task_method = c("lambda_blend", "gev_patch", "core_hatsa"),
    lambda_blend_value = 0.15,
    k_gev_dims = 10,
    row_augmentation = TRUE, # Default TRUE, but logic checks if task_data allows it
    residualize_condition_anchors = FALSE,
    omega_weights = NULL, # Defaults handled in weighted solver
    omega_mode = c("fixed", "adaptive"),
    reliability_scores_list = NULL,
    scale_omega_trace = TRUE,
    alpha_laplacian = 0.93,
    degree_type_laplacian = c("abs", "positive", "signed"),
    k_conn_pos = 10, # Default based on typical usage
    k_conn_neg = 10, # Default based on typical usage
    k_conn_task_pos = 10, # Default, align with k_conn_pos/neg
    k_conn_task_neg = 10, # Default, align with k_conn_pos/neg
    similarity_method_task = "pearson",
    W_task_helper_func = NULL, # Will default based on task_data_list structure
    n_refine = 5, # Default based on typical usage
    check_redundancy = TRUE,
    redundancy_threshold = 0.45,
    residualize_k_conn_proj = 64,
    residualize_k_conn_labels = 10,
    gev_lambda_max = 0.8,
    gev_epsilon_reg = 1e-6,
    parcel_names = NULL,
    verbose = TRUE
) {
    # Match arguments from vectors of choices
    task_method <- match.arg(task_method)
    omega_mode <- match.arg(omega_mode)
    degree_type_laplacian <- match.arg(degree_type_laplacian)
    
    # --- Step 1: Validate arguments and initialize variables ---
    args <- validate_and_initialize_args(
        subject_data_list = subject_data_list,
        anchor_indices = anchor_indices,
            spectral_rank_k = spectral_rank_k,
        task_data_list = task_data_list,
            task_method = task_method,
            lambda_blend_value = lambda_blend_value,
            k_gev_dims = k_gev_dims,
            row_augmentation = row_augmentation,
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
        W_task_helper_func = W_task_helper_func,
            n_refine = n_refine,
            check_redundancy = check_redundancy,
            redundancy_threshold = redundancy_threshold,
            residualize_k_conn_proj = residualize_k_conn_proj,
            residualize_k_conn_labels = residualize_k_conn_labels,
            gev_lambda_max = gev_lambda_max,
            gev_epsilon_reg = gev_epsilon_reg,
        parcel_names = parcel_names,
        verbose = verbose
    )

    # --- Step 2: Process each subject for graph construction and basis shaping ---
    processing_results <- process_subjects(
        args = args,
        subject_data_list = subject_data_list,
        task_data_list = task_data_list,
        anchor_indices = anchor_indices
    )

    # --- Step 3: Anchor Augmentation ---
    augmentation_results <- perform_anchor_augmentation(
        args = args,
        processing_results = processing_results, 
        anchor_indices = anchor_indices
    )

    # --- Step 4: Iterative Refinement (GPA) ---
    gpa_results <- prepare_and_run_gpa(
        args = args,
        augmentation_results = augmentation_results,
        reliability_scores_list = reliability_scores_list
    )

    # --- Step 5: Patch Alignment (If GEV) ---
    patch_results <- perform_patch_alignment(
        args = args,
        processing_results = processing_results,
        gpa_results = gpa_results
    )

    # --- Step 6: Construct Output ---
    result <- construct_output(
        args = args,
        processing_results = processing_results,
        augmentation_results = augmentation_results,
        gpa_results = gpa_results,
        patch_results = patch_results
    )

    return(result)
}

# # Helper for status messages (can be moved to utils)
# message_stage <- function(message_text, verbose = TRUE) {
#     if (verbose) {
#         message(rep("-", nchar(message_text)))
#         message(message_text)
#         message(rep("-", nchar(message_text)))
#     }
# }

#' Advanced Options for Task-Informed HATSA
#'
#' Creates a list of advanced parameters for the task_hatsa function.
#'
#' @param lambda_blend_value Numeric `lambda` in `[0,1]`. Weight for `L_task` in blend. Default 0.15.
#' @param k_gev_dims Integer, requested dimension for GEV patches. Default 10. Used if `task_method == "gev_patch"`.
#' @param row_augmentation Logical. If `TRUE`, add projected task features to anchor matrices
#'   for GPA refinement. Requires suitable `task_data_list`. Default `TRUE` if suitable data provided.
#' @param residualize_condition_anchors Logical. If `TRUE` and `row_augmentation` is `TRUE`,
#'   residualize projected task anchors against parcel anchors. Default `FALSE`.
#' @param omega_weights List specifying fixed weights for weighted Procrustes (e.g.,
#'   `list(parcel = 1.0, condition = 0.5)`). Used if `row_augmentation=TRUE` and `omega_mode == "fixed"`.
#'   Defaults handled by `solve_procrustes_rotation_weighted`.
#' @param omega_mode Character string: `"fixed"` or `"adaptive"`. Controls weighting in GPA. Default `"fixed"`.
#' @param reliability_scores_list List (parallel to `subject_data_list`), each element a numeric
#'   vector of reliability scores (e.g., R^2) for task data (length `C`). Used if `omega_mode == "adaptive"`.
#' @param scale_omega_trace Logical. Whether to rescale weights in weighted GPA so trace equals total anchors. Default `TRUE`.
#' @param alpha_laplacian Numeric, laziness parameter for graph Laplacians (`L = I - alpha D^{-1} W`). Default 0.93.
#' @param degree_type_laplacian Character string (`"abs"`, `"positive"`, `"signed"`). Type of degree calculation for Laplacian. Default `"abs"`.
#' @param k_conn_pos Integer >= 0. k-NN sparsification for positive edges in `W_conn`.
#' @param k_conn_neg Integer >= 0. k-NN sparsification for negative edges in `W_conn`.
#' @param k_conn_task_pos Integer >= 0. k-NN sparsification for positive edges in `W_task`.
#' @param k_conn_task_neg Integer >= 0. k-NN sparsification for negative edges in `W_task`.
#' @param similarity_method_task Character string or function. Method to compute similarity for `W_task`
#'   (e.g., "pearson", "spearman"). Default "pearson".
#' @param W_task_helper_func Function. The specific function to compute `W_task` (e.g.,
#'   `compute_W_task_from_activations`, `compute_W_task_from_encoding`). If `NULL`, attempts
#'   to infer based on `task_data_list` structure (currently assumes activations `C x Vp`). Default `NULL`.
#' @param n_refine Integer >= 0. Number of GPA refinement iterations.
#' @param check_redundancy Logical. If `TRUE`, check correlation between `W_conn` and `W_task`. Default `TRUE`.
#' @param redundancy_threshold Numeric. Spearman rho threshold for triggering `W_task` residualization. Default 0.45.
#' @param residualize_k_conn_proj Integer. Number of `L_conn` eigenvectors to project `W_task` out of. Default 64.
#' @param residualize_k_conn_labels Integer. k-NN value for re-sparsifying `W_task_res` after residualization. Default 10.
#' @param gev_lambda_max Numeric. Max GEV eigenvalue `lambda` to retain for patches. Default 0.8.
#' @param gev_epsilon_reg Numeric. Small regularization for `L_conn` in GEV. Default 1e-6.
#' @param parcel_names Optional character vector of parcel names. If `NULL`, names like "P1", "P2"... are generated.
#'
#' @return A list of options to pass to the task_hatsa function.
#'
#' @export
task_hatsa_opts <- function(
    lambda_blend_value = 0.15,
    k_gev_dims = 10,
    row_augmentation = TRUE,
    residualize_condition_anchors = FALSE,
    omega_weights = NULL,
    omega_mode = c("fixed", "adaptive"),
    reliability_scores_list = NULL,
    scale_omega_trace = TRUE,
    alpha_laplacian = 0.93,
    degree_type_laplacian = c("abs", "positive", "signed"),
    k_conn_pos = 10,
    k_conn_neg = 10,
    k_conn_task_pos = 10,
    k_conn_task_neg = 10,
    similarity_method_task = "pearson",
    W_task_helper_func = NULL,
    n_refine = 5,
    check_redundancy = TRUE,
    redundancy_threshold = 0.45,
    residualize_k_conn_proj = 64,
    residualize_k_conn_labels = 10,
    gev_lambda_max = 0.8,
    gev_epsilon_reg = 1e-6,
    parcel_names = NULL
) {
    omega_mode <- match.arg(omega_mode)
    degree_type_laplacian <- match.arg(degree_type_laplacian)
    
    list(
        lambda_blend_value = lambda_blend_value,
        k_gev_dims = k_gev_dims,
        row_augmentation = row_augmentation,
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
        W_task_helper_func = W_task_helper_func,
        n_refine = n_refine,
        check_redundancy = check_redundancy,
        redundancy_threshold = redundancy_threshold,
        residualize_k_conn_proj = residualize_k_conn_proj,
        residualize_k_conn_labels = residualize_k_conn_labels,
        gev_lambda_max = gev_lambda_max,
        gev_epsilon_reg = gev_epsilon_reg,
        parcel_names = parcel_names
    )
}

#' Hyperalignment via Task-Informed Shared Analysis
#' 
#' A user-friendly interface to run Task-Informed HATSA (task_hatsa), which incorporates 
#' task information during basis shaping and/or alignment refinement.
#'
#' @param subject_data_list List where each element is a numeric matrix (`T_i x V_p`)
#'   of time-series data for one subject.
#' @param anchor_indices Integer vector, indices of the canonical anchor parcels (`1` to `V_p`).
#' @param spectral_rank_k Integer, the desired dimensionality of the primary spectral sketch.
#'   Default is 40.
#' @param task_data_list List (parallel to `subject_data_list`). Each element provides
#'   task-related data for a subject. Structure depends on `task_method` and subsequent
#'   processing. Default is NULL.
#' @param task_method Character string: `"lambda_blend"`, `"gev_patch"`, or `"core_hatsa"`. 
#'   Default: `"lambda_blend"`.
#' @param opts List of advanced options created by `task_hatsa_opts()`.
#' @param verbose Logical. Print progress messages? Default `TRUE`.
#' @param ... Additional arguments passed to `task_hatsa_opts()` if not provided directly via `opts`.
#'
#' @return A list representing the `task_hatsa_projector` object containing aligned bases, 
#'         rotations, parameters, QC metrics, etc.
#'
#' @export
task_hatsa <- function(
    subject_data_list,
    anchor_indices,
    spectral_rank_k = 40,
    task_data_list = NULL,
    task_method = c("lambda_blend", "gev_patch", "core_hatsa"),
    opts = task_hatsa_opts(),
    verbose = TRUE,
    ...
) {
    task_method <- match.arg(task_method)
    
    # Handle additional options passed via ...
    if (...length() > 0) {
        additional_opts <- list(...)
        # Override existing opts with additional options
        for (opt_name in names(additional_opts)) {
            opts[[opt_name]] <- additional_opts[[opt_name]]
        }
    }
    
    # Call the existing implementation with all parameters
    run_task_hatsa(
        subject_data_list = subject_data_list,
        anchor_indices = anchor_indices,
        spectral_rank_k = spectral_rank_k,
        task_data_list = task_data_list,
        task_method = task_method,
        lambda_blend_value = opts$lambda_blend_value,
        k_gev_dims = opts$k_gev_dims,
        row_augmentation = opts$row_augmentation,
        residualize_condition_anchors = opts$residualize_condition_anchors,
        omega_weights = opts$omega_weights,
        omega_mode = opts$omega_mode,
        reliability_scores_list = opts$reliability_scores_list,
        scale_omega_trace = opts$scale_omega_trace,
        alpha_laplacian = opts$alpha_laplacian,
        degree_type_laplacian = opts$degree_type_laplacian,
        k_conn_pos = opts$k_conn_pos,
        k_conn_neg = opts$k_conn_neg,
        k_conn_task_pos = opts$k_conn_task_pos,
        k_conn_task_neg = opts$k_conn_task_neg,
        similarity_method_task = opts$similarity_method_task,
        W_task_helper_func = opts$W_task_helper_func,
        n_refine = opts$n_refine,
        check_redundancy = opts$check_redundancy,
        redundancy_threshold = opts$redundancy_threshold,
        residualize_k_conn_proj = opts$residualize_k_conn_proj,
        residualize_k_conn_labels = opts$residualize_k_conn_labels,
        gev_lambda_max = opts$gev_lambda_max,
        gev_epsilon_reg = opts$gev_epsilon_reg,
        parcel_names = opts$parcel_names,
        verbose = verbose
    )
}

# Now we need to save the original run_task_hatsa implementation
# and redefine it to call task_hatsa, maintaining backward compatibility

# First, rename the original function to .task_hatsa_engine 
.task_hatsa_engine <- run_task_hatsa

# Now redefine run_task_hatsa to call task_hatsa while preserving the original signature
#' @rdname task_hatsa
#' @export
run_task_hatsa <- function(
    subject_data_list,
    anchor_indices,
    spectral_rank_k,
    task_data_list = NULL,
    task_method = c("lambda_blend", "gev_patch", "core_hatsa"),
    lambda_blend_value = 0.15,
    k_gev_dims = 10,
    row_augmentation = TRUE,
    residualize_condition_anchors = FALSE,
    omega_weights = NULL,
    omega_mode = c("fixed", "adaptive"),
    reliability_scores_list = NULL,
    scale_omega_trace = TRUE,
    alpha_laplacian = 0.93,
    degree_type_laplacian = c("abs", "positive", "signed"),
    k_conn_pos = 10,
    k_conn_neg = 10,
    k_conn_task_pos = 10,
    k_conn_task_neg = 10,
    similarity_method_task = "pearson",
    W_task_helper_func = NULL,
    n_refine = 5,
    check_redundancy = TRUE,
    redundancy_threshold = 0.45,
    residualize_k_conn_proj = 64,
    residualize_k_conn_labels = 10,
    gev_lambda_max = 0.8,
    gev_epsilon_reg = 1e-6,
    parcel_names = NULL,
    verbose = TRUE
) {
    # Print a deprecation message but don't break anything
    .Deprecated("task_hatsa", msg = "run_task_hatsa() is kept for backward compatibility; new code should call task_hatsa()")
    
    # Call the original implementation
    .task_hatsa_engine(
        subject_data_list = subject_data_list,
        anchor_indices = anchor_indices,
        spectral_rank_k = spectral_rank_k,
        task_data_list = task_data_list,
        task_method = task_method,
        lambda_blend_value = lambda_blend_value,
        k_gev_dims = k_gev_dims,
        row_augmentation = row_augmentation,
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
        W_task_helper_func = W_task_helper_func,
        n_refine = n_refine,
        check_redundancy = check_redundancy,
        redundancy_threshold = redundancy_threshold,
        residualize_k_conn_proj = residualize_k_conn_proj,
        residualize_k_conn_labels = residualize_k_conn_labels,
        gev_lambda_max = gev_lambda_max,
        gev_epsilon_reg = gev_epsilon_reg,
        parcel_names = parcel_names,
        verbose = verbose
    )
}
