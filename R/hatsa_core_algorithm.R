#' Run the Core HATSA Algorithm
#'
#' Implements the Core HATSA algorithm to align functional connectivity patterns
#' across subjects. This version uses sparse matrices for graph representations,
#' efficient eigendecomposition via `PRIMME`, and incorporates robustness
#' improvements based on detailed audits.
#'
#' @param subject_data_list A list of dense numeric matrices. Each matrix `X_i`
#'   corresponds to a subject, with dimensions `T_i` (time points) x `V_p` (parcels).
#' @param anchor_indices A numeric vector of 1-based indices for the selected
#'   anchor parcels. Duplicate indices will be removed.
#' @param spectral_rank_k An integer specifying the dimensionality (`k`) of the
#'   low-dimensional spectral sketch. Must be non-negative. `k=0` yields empty sketches.
#' @param k_conn_pos An integer. For graph sparsification, number of strongest
#'   positive connections to retain per parcel.
#' @param k_conn_neg An integer. For graph sparsification, number of strongest
#'   negative connections to retain per parcel.
#' @param n_refine An integer, number of GPA refinement iterations.
#' @param use_dtw Logical, defaults to `FALSE`. If `TRUE` (not yet fully implemented),
#'   Dynamic Time Warping would be considered in graph construction similarity.
#' @param n_cores Integer number of CPU cores to use. If `> 1` and the platform
#'   supports forking (i.e., non-Windows), per-subject computations are
#'   parallelized via \code{parallel::mclapply}. Defaults to 1 (sequential).
#'
#' @return A `hatsa_projector` object. This S3 object inherits from
#'   `multiblock_biprojector` (from the `multivarious` package) and contains
#'   the results of the HATSA analysis. Key components include:
#'   \itemize{
#'     \item{\code{v}: The mean aligned sketch (group-level template, V_p x k matrix).}
#'     \item{\code{s}: Stacked aligned sketches for all subjects ((N*V_p) x k matrix).}
#'     \item{\code{sdev}: Component standard deviations (vector of length k, currently defaults to 1s).}
#'     \item{\code{preproc}: Preprocessing object (currently `prep(pass())`).}
#'     \item{\code{block_indices}: List defining subject blocks in the `s` matrix.}
#'     \item{\code{R_final_list}: List of subject-specific rotation matrices (k x k).}
#'     \item{\code{U_original_list}: List of subject-specific original (unaligned) sketch matrices (V_p x k).}
#'     \item{\code{Lambda_original_list}: List of subject-specific original eigenvalues (vector of length k) from the parcel-level decomposition. Essential for Nyström voxel projection.}
#'     \item{\code{T_anchor_final}: The final group anchor template used for alignment (V_a x k matrix, where V_a is number of anchors).}
#'     \item{\code{parameters}: List of input parameters used for the HATSA run (e.g., `k`, `V_p`, `N_subjects`, anchor details, sparsification parameters).}
#'     \item{\code{method}: Character string, "hatsa_core".}
#'     \item{\code{U_aligned_list}: (Note: while used to compute `v` and `s`, direct aligned sketches per subject are also stored if `project_block` needs them or for direct inspection, typically same as `object$s` reshaped per block)}
#'   }
#'   This object can be used with S3 methods like `print`, `summary`, `coef`,
#'   `scores`, `predict` (for new parcel data), and `project_voxels` (for new
#'   voxel data).
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' N_subjects <- 3 # Small N for quick example
#' V_p_parcels <- 25
#' T_times_avg <- 50
#'
#' # Create a list of matrices (time x parcels)
#' subject_data <- lapply(1:N_subjects, function(i) {
#'   T_i <- T_times_avg + sample(-5:5, 1)
#'   matrix(stats::rnorm(T_i * V_p_parcels), nrow = T_i, ncol = V_p_parcels)
#' })
#'
#' # Assign parcel names (optional but good practice)
#' parcel_names_vec <- paste0("Parcel_", 1:V_p_parcels)
#' subject_data <- lapply(subject_data, function(mat) {
#'   colnames(mat) <- parcel_names_vec
#'   mat
#' })
#'
#' # Define HATSA parameters
#' # Ensure number of anchors >= k for stable Procrustes
#' anchor_idx <- sample(1:V_p_parcels, 7)
#' k_spectral <- 5 # k=5, num_anchors=7 is valid
#' k_pos <- 4
#' k_neg <- 2
#' n_iter_refine <- 2
#'
#' # Run Core HATSA (requires Matrix and PRIMME packages)
#' hatsa_results <- NULL
#' if (requireNamespace("Matrix", quietly = TRUE) &&
#'     requireNamespace("PRIMME", quietly = TRUE)) {
#'   hatsa_results <- tryCatch(
#'     run_hatsa_core(
#'       subject_data_list = subject_data,
#'       anchor_indices = anchor_idx,
#'       spectral_rank_k = k_spectral,
#'       k_conn_pos = k_pos,
#'       k_conn_neg = k_neg,
#'       n_refine = n_iter_refine
#'     ),
#'     error = function(e) {
#'       message("Example run failed: ", e$message)
#'       NULL
#'     }
#'   )
#'
#'   # Inspect the results object
#'   if (!is.null(hatsa_results)) {
#'     print(hatsa_results)
#'     summary_info <- summary(hatsa_results)
#'     print(summary_info)
#'
#'     # Get coefficients (mean aligned sketch)
#'     group_template <- coef(hatsa_results)
#'     # print(dim(group_template)) # Should be V_p x k
#'
#'     # Get stacked scores (aligned sketches for all subjects)
#'     all_scores <- scores(hatsa_results)
#'     # print(dim(all_scores)) # Should be (N*V_p) x k
#'
#'     # Get block indices to map scores to subjects
#'     indices <- block_indices(hatsa_results)
#'     # subject1_scores <- all_scores[indices[[1]], ]
#'     # print(dim(subject1_scores)) # Should be V_p x k
#'   }
#' } else {
#'   if (interactive()) message("Matrix and PRIMME packages needed for this example.")
#' }
#'
#' @seealso \code{\link{hatsa_projector}}, \code{\link{project_voxels.hatsa_projector}}
#' @author Expert R Developer (GPT)
#' @export
#' @importFrom stats setNames rnorm runif sd
#' @importFrom multivarious pass prep
#' @importFrom parallel mclapply detectCores
run_hatsa_core <- function(subject_data_list, anchor_indices, spectral_rank_k,
                           k_conn_pos, k_conn_neg, n_refine, use_dtw = FALSE,
                           n_cores = 1L) {

  if (length(subject_data_list) > 0 && !is.null(subject_data_list[[1]])) {
    stopifnot(is.matrix(subject_data_list[[1]]), !inherits(subject_data_list[[1]], "Matrix"))
    V_p <- ncol(subject_data_list[[1]])
    pnames <- colnames(subject_data_list[[1]])
    if (is.null(pnames) || length(pnames) != V_p) pnames <- paste0("P", 1:V_p)
  } else { 
    V_p <- 0
    pnames <- character(0)
  }

  val_results <- validate_hatsa_inputs(subject_data_list, anchor_indices, spectral_rank_k,
                                       k_conn_pos, k_conn_neg, n_refine, V_p)
  unique_anchor_indices <- val_results$unique_anchor_indices
  
  num_subjects <- length(subject_data_list)
  n_cores <- min(as.integer(n_cores), parallel::detectCores())
  if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
  U_original_list <- vector("list", num_subjects)
  Lambda_original_list <- vector("list", num_subjects)
  Lambda_original_gaps_list <- vector("list", num_subjects)
  
  message_stage("Stage 1: Computing initial spectral sketches...", interactive_only = TRUE)
  if (num_subjects > 0) {
    process_one <- function(i) {
      X_i <- subject_data_list[[i]]
      current_pnames <- colnames(X_i)
      if (is.null(current_pnames) || length(current_pnames) != V_p) current_pnames <- pnames

      W_conn_i_sparse <- compute_subject_connectivity_graph_sparse(X_i, current_pnames,
                                                                  k_conn_pos, k_conn_neg, use_dtw)
      L_conn_i_sparse <- compute_graph_laplacian_sparse(W_conn_i_sparse)

      sketch_result <- compute_spectral_sketch_sparse(L_conn_i_sparse, spectral_rank_k, eigenvalue_tol = 1e-8)
      lambdas_i <- sketch_result$values
      gaps_i <- if (!is.null(lambdas_i) && length(lambdas_i) > 1) {
        denominators <- lambdas_i[-length(lambdas_i)]
        denominators[abs(denominators) < .Machine$double.eps^0.5] <- NA
        (lambdas_i[-1] - lambdas_i[-length(lambdas_i)]) / denominators
      } else {
        numeric(0)
      }
      list(U = sketch_result$vectors, Lambda = lambdas_i, gaps = gaps_i)
    }

    res_list <- if (.Platform$OS.type == "windows" || n_cores == 1L) {
      lapply(seq_len(num_subjects), process_one)
    } else {
      parallel::mclapply(seq_len(num_subjects), process_one, mc.cores = n_cores)
    }

    for (i in seq_len(num_subjects)) {
      U_original_list[[i]] <- res_list[[i]]$U
      Lambda_original_list[[i]] <- res_list[[i]]$Lambda
      Lambda_original_gaps_list[[i]] <- res_list[[i]]$gaps
    }
  }
  message_stage("Stage 1 complete.", interactive_only = TRUE)

  message_stage("Stage 2: Performing iterative refinement (GPA)...", interactive_only = TRUE)
  A_originals_list <- lapply(U_original_list, function(U_orig_subj) {
    if (is.null(U_orig_subj) || nrow(U_orig_subj) == 0 || ncol(U_orig_subj) == 0 || length(unique_anchor_indices) == 0) {
      return(matrix(0, nrow = length(unique_anchor_indices), ncol = spectral_rank_k))
    }
    U_orig_subj[unique_anchor_indices, , drop = FALSE]
  })
  
  gpa_results <- perform_gpa_refinement(A_originals_list, n_refine, spectral_rank_k)
  R_final_list <- gpa_results$R_final_list
  T_anchor_final <- gpa_results$T_anchor_final
  message_stage("Stage 2 complete.", interactive_only = TRUE)

  message_stage("Stage 3: Applying final rotations...", interactive_only = TRUE)
  U_aligned_list <- vector("list", num_subjects)
  if (num_subjects > 0) {
    rotate_one <- function(i) {
      U_orig_i <- U_original_list[[i]]
      R_final_i <- R_final_list[[i]]
      if (!is.null(U_orig_i) && nrow(U_orig_i) == V_p && ncol(U_orig_i) == spectral_rank_k &&
          !is.null(R_final_i) && nrow(R_final_i) == spectral_rank_k && ncol(R_final_i) == spectral_rank_k) {
        U_orig_i %*% R_final_i
      } else {
        matrix(0, nrow = V_p, ncol = spectral_rank_k)
      }
    }

    U_aligned_list <- if (.Platform$OS.type == "windows" || n_cores == 1L) {
      lapply(seq_len(num_subjects), rotate_one)
    } else {
      parallel::mclapply(seq_len(num_subjects), rotate_one, mc.cores = n_cores)
    }
  }
  message_stage("Stage 3 complete. Core HATSA finished.", interactive_only = TRUE)
  
  hatsa_core_results <- list(
    U_aligned_list = U_aligned_list,
    R_final_list = R_final_list,
    U_original_list = U_original_list,
    Lambda_original_list = Lambda_original_list, 
    Lambda_original_gaps_list = Lambda_original_gaps_list,
    T_anchor_final = T_anchor_final
  )

  parameters <- list(
    k = spectral_rank_k,
    N_subjects = num_subjects,
    V_p = V_p, 
    method = "hatsa_core",
    anchor_indices = unique_anchor_indices,
    k_conn_pos = k_conn_pos,
    k_conn_neg = k_conn_neg,
    n_refine = n_refine,
    use_dtw = use_dtw,
    n_cores = n_cores
  )
  
  projector_object <- hatsa_projector(hatsa_core_results = hatsa_core_results, 
                                      parameters = parameters)
  
  return(projector_object)
}

