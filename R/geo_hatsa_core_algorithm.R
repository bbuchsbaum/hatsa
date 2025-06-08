#' Run Geometric Harmonized Tensors SVD Alignment (Geo-HATSA) Core Algorithm
#'
#' This function implements a variant of the HATSA core algorithm where the
#' Generalized Procrustes Analysis (GPA) step for refining rotations and the
#' group anchor template uses a geometric approach on SO(k) via
#' `perform_geometric_gpa_refinement`.
#'
#' @param subject_data_list A list of subject-level parcel time-series matrices.
#'   Each element is a numeric matrix (T_i time points x V_p parcels).
#' @param anchor_indices A numeric vector of indices for the anchor parcels
#'   (1-based, referring to columns of matrices in `subject_data_list`).
#' @param spectral_rank_k Integer, the number of spectral components (k) to retain.
#' @param k_conn_pos Integer, number of positive degree neighbors for k-NN graph construction.
#' @param k_conn_neg Integer, number of negative degree neighbors (if applicable, often same as pos).
#' @param n_refine Integer, number of GPA refinement iterations.
#' @param V_p Integer, number of parcels (vertices per subject). If NULL, inferred.
#' @param use_dtw Logical, whether to use Dynamic Time Warping for FC graph (passed to
#'   `compute_subject_connectivity_graph_sparse`). Default FALSE.
#' @param graph_mode Character string, specifies the graph construction method.
#'   Options include `"anchor_block"` (default), `"schur"` (Schur complement),
#'   `"full"` (standard kNN graph). Passed to
#'   `compute_subject_connectivity_graph_sparse`.
#' @param schur_eps Numeric, epsilon for Schur complement graph construction if
#'   `graph_mode = "schur"`. Default: 0.01.
#' @param eigengap_tol Numeric, tolerance for checking the eigengap after spectral
#'   sketch computation. If `lambda_k - lambda_{k+1} < eigengap_tol`, a warning
#'   is issued for that subject. Default: 1e-9.
#' @param rotation_mode Character string, passed to `perform_geometric_gpa_refinement`.
#'   One of `"svd"` (default) or `"riemannian"`. See that function for details.
#' @param verbose Logical, if TRUE, prints progress messages. Default TRUE.
#' @param frechet_mean_options A list of options to pass to `hatsa::frechet_mean_so_k`
#'   within `perform_geometric_gpa_refinement`. Default: `list()`.
#' @param gpa_tol Numeric, tolerance for convergence in `perform_geometric_gpa_refinement`.
#'   Default: 1e-7.
#'
#' @return A list containing the core Geo-HATSA results and QC metrics:
#'   \itemize{
#'     \item{\code{U_aligned_list}: List of subject-specific aligned sketch matrices (V_p x k).}
#'     \item{\code{R_final_list}: List of subject-specific rotation matrices (k x k).}
#'     \item{\code{U_original_list}: List of subject-specific original sketch matrices (V_p x k).}
#'     \item{\code{Lambda_original_list}: List of subject-specific original eigenvalues (length k).}
#'     \item{\code{Lambda_original_gaps_list}: List of subject-specific eigengap ratios.}
#'     \item{\code{T_anchor_final}: The final group anchor template matrix (N_anchors x k).}
#'     \item{\code{R_bar_final}: The final Fréchet mean of rotation matrices from Geo-GPA.}
#'     \item{\code{qc_metrics}: A data frame with per-subject QC flags (laplacian_computed_ok, 
#'           sketch_computed_ok, eigengap_sufficient).}
#'   }
#'   This list is intended to be used as input to the `hatsa_projector` constructor,
#'   though `R_bar_final` is an additional output specific to Geo-HATSA.
#'
#' @export
#' @seealso \code{\link{perform_geometric_gpa_refinement}}, \code{\link{hatsa_projector}}
#' @importFrom stats cor
run_geo_hatsa_core <- function(subject_data_list, 
                               anchor_indices, 
                               spectral_rank_k, 
                               k_conn_pos = 10, k_conn_neg = 10, 
                               n_refine = 10, 
                               V_p = NULL,
                               use_dtw = FALSE,
                               graph_mode = "anchor_block",
                               schur_eps = 0.01,
                               eigengap_tol = 1e-9,
                               rotation_mode = "svd",
                               verbose = TRUE,
                               frechet_mean_options = list(),
                               gpa_tol = 1e-7) {

  if (verbose) message_stage("Starting Geometric HATSA (Geo-HATSA) Core Algorithm...", interactive_only = TRUE)

  # --- Input Validation (simplified from run_hatsa_core, focusing on key aspects) ---
  if (!is.list(subject_data_list) || length(subject_data_list) == 0) {
    stop("subject_data_list must be a non-empty list of matrices.")
  }
  if (is.null(V_p)) {
    first_subj_data <- Filter(Negate(is.null), subject_data_list)
    if (length(first_subj_data) > 0 && is.matrix(first_subj_data[[1]])) {
      V_p <- ncol(first_subj_data[[1]])
    } else {
      stop("V_p not specified and cannot infer from first subject data.")
    }
  }
  if (V_p == 0) stop("V_p (number of parcels) cannot be zero.")
  
  pnames <- paste0("P", 1:V_p) # Default parcel names
  val_results <- validate_hatsa_inputs(subject_data_list, anchor_indices, spectral_rank_k,
                                       k_conn_pos, k_conn_neg, n_refine, V_p) # Reuse existing validator
  unique_anchor_indices <- val_results$unique_anchor_indices
  num_anchors <- length(unique_anchor_indices)

  num_subjects <- length(subject_data_list)
  U_original_list <- vector("list", num_subjects)
  Lambda_original_list <- vector("list", num_subjects)
  Lambda_original_gaps_list <- vector("list", num_subjects)

  # Initialize QC metrics dataframe
  qc_metrics <- data.frame(
      subject_id = seq_len(num_subjects),
      laplacian_computed_ok = TRUE,
      sketch_computed_ok = TRUE,
      eigengap_sufficient = TRUE,
      stringsAsFactors = FALSE
  )

  # --- Stage 1: Compute initial spectral sketches (same as run_hatsa_core) ---
  if (verbose) message_stage("Geo-HATSA Stage 1: Computing initial spectral sketches...", interactive_only = TRUE)
  if (num_subjects > 0) {
    for (i in 1:num_subjects) {
      X_i <- subject_data_list[[i]]
      if (is.null(X_i) || !is.matrix(X_i) || ncol(X_i) != V_p) {
          message(sprintf("  Geo-HATSA Stage 1: Subject %d has invalid data. Skipping.", i))
          qc_metrics$laplacian_computed_ok[i] <- FALSE
          qc_metrics$sketch_computed_ok[i] <- FALSE
          U_original_list[[i]] <- matrix(NA, nrow=V_p, ncol=spectral_rank_k)
          Lambda_original_list[[i]] <- rep(NA, spectral_rank_k)
          next
      }
      current_pnames <- colnames(X_i)
      if(is.null(current_pnames) || length(current_pnames) != V_p) current_pnames <- pnames

      W_conn_i_sparse <- tryCatch({
          compute_subject_connectivity_graph_sparse(X_i, current_pnames, k_conn_pos, k_conn_neg, use_dtw, 
                                                    graph_mode = graph_mode, schur_options = list(eps = schur_eps))
      }, error = function(e) {
          message(sprintf("  Geo-HATSA Stage 1: Error computing connectivity graph for subject %d: %s", i, e$message))
          qc_metrics$laplacian_computed_ok[i] <- FALSE
          NULL
      })
      
      if (is.null(W_conn_i_sparse)) {
          qc_metrics$sketch_computed_ok[i] <- FALSE
          U_original_list[[i]] <- matrix(NA, nrow=V_p, ncol=spectral_rank_k)
          Lambda_original_list[[i]] <- rep(NA, spectral_rank_k)
          next
      }
      L_conn_i_sparse <- compute_graph_laplacian_sparse(W_conn_i_sparse)
      
      sketch_result <- tryCatch({
          compute_spectral_sketch_sparse(L_conn_i_sparse, spectral_rank_k, eigenvalue_tol = 1e-8)
      }, error = function(e) {
          message(sprintf("  Geo-HATSA Stage 1: Error computing spectral sketch for subject %d: %s", i, e$message))
          qc_metrics$sketch_computed_ok[i] <- FALSE
          NULL
      })
      
      if (is.null(sketch_result) || is.null(sketch_result$vectors) || is.null(sketch_result$values)){
          U_original_list[[i]] <- matrix(NA, nrow=V_p, ncol=spectral_rank_k)
          Lambda_original_list[[i]] <- rep(NA, spectral_rank_k)
          qc_metrics$sketch_computed_ok[i] <- FALSE # Already set but good to be explicit
          next
      }
      U_original_list[[i]] <- sketch_result$vectors
      Lambda_original_list[[i]] <- sketch_result$values

      lambdas_i <- Lambda_original_list[[i]]
      if (!is.null(lambdas_i) && length(lambdas_i) > 1) {
        # Eigengap check
        if (spectral_rank_k < length(lambdas_i)) {
            eigengap <- lambdas_i[spectral_rank_k] - lambdas_i[spectral_rank_k + 1]
            if (eigengap < eigengap_tol) {
                qc_metrics$eigengap_sufficient[i] <- FALSE
                if(verbose) message(sprintf("  Geo-HATSA Stage 1: Warning for subject %d: Eigengap (%.2e) at k=%d is below tolerance (%.2e). Rotations might be unstable.", 
                                    i, eigengap, spectral_rank_k, eigengap_tol))
            }
        } else if (length(lambdas_i) == spectral_rank_k && spectral_rank_k > 0) {
            # Cannot compute k+1th eigenvalue, assume gap is sufficient if sketch dim matches request and is >0
            # Or consider it insufficient if strict k vs k+1 gap is required. For now, mark as sufficient.
             qc_metrics$eigengap_sufficient[i] <- TRUE 
        } else {
            qc_metrics$eigengap_sufficient[i] <- FALSE # Not enough eigenvalues for rank k
        }
        
        denominators <- lambdas_i[-length(lambdas_i)]
        denominators[abs(denominators) < .Machine$double.eps^0.5] <- NA # Avoid division by zero
        gaps_i <- (lambdas_i[-1] - lambdas_i[-length(lambdas_i)]) / denominators
        Lambda_original_gaps_list[[i]] <- gaps_i
      } else {
        Lambda_original_gaps_list[[i]] <- numeric(0)
        qc_metrics$eigengap_sufficient[i] <- FALSE # Not enough eigenvalues
      }
      
      if (num_subjects > 10 && i %% floor(num_subjects/5) == 0 && interactive() && verbose) {
          message(sprintf("  Geo-HATSA Stage 1: Processed subject %d/%d", i, num_subjects))
      }
    }
  }
  if (verbose) message_stage("Geo-HATSA Stage 1 complete.", interactive_only = TRUE)

  # --- Stage 2: Performing iterative GEOMETRIC refinement (Geo-GPA) ---
  if (verbose) message_stage("Geo-HATSA Stage 2: Performing iterative geometric refinement (Geo-GPA)...", interactive_only = TRUE)
  
  # Prepare A_originals_list (anchor sketches)
  A_originals_list <- lapply(U_original_list, function(U_orig_subj) {
    if (is.null(U_orig_subj) || !is.matrix(U_orig_subj) || nrow(U_orig_subj) != V_p || ncol(U_orig_subj) != spectral_rank_k || num_anchors == 0) {
      return(matrix(NA, nrow = num_anchors, ncol = spectral_rank_k))
    }
    U_orig_subj[unique_anchor_indices, , drop = FALSE]
  })
  
  # Filter out any NA matrices from A_originals_list that resulted from failed U_original_list entries
  valid_A_indices_for_gpa <- !sapply(A_originals_list, function(A) any(is.na(A)) || nrow(A) != num_anchors || ncol(A) != spectral_rank_k)
  A_originals_for_gpa <- A_originals_list[valid_A_indices_for_gpa]
  
  initial_R_for_gpa <- replicate(sum(valid_A_indices_for_gpa), diag(spectral_rank_k), simplify=FALSE) # Default

  if(length(A_originals_for_gpa) < 2 && verbose){
      message("Geo-HATSA Stage 2: Fewer than 2 subjects have valid anchor sketches for GPA. Skipping refinement.")
      # Create dummy geo_gpa_results if GPA is skipped
      R_final_list_placeholder <- vector("list", num_subjects)
      for(idx in 1:num_subjects) R_final_list_placeholder[[idx]] <- if(qc_metrics$sketch_computed_ok[idx]) diag(spectral_rank_k) else matrix(NA, spectral_rank_k, spectral_rank_k)
      geo_gpa_results <- list(
          R_final_list = R_final_list_placeholder,
          T_anchor_final = matrix(NA, nrow = num_anchors, ncol = spectral_rank_k),
          R_bar_final = diag(spectral_rank_k)
      )
  } else if (length(A_originals_for_gpa) >= 2) {
      geo_gpa_results <- perform_geometric_gpa_refinement(
        A_originals_list = A_originals_for_gpa, 
        n_refine = n_refine,
        k = spectral_rank_k,
        m_rows = num_anchors, 
        tol = gpa_tol,
        rotation_mode = rotation_mode,
        frechet_mean_options = frechet_mean_options,
        verbose = verbose,
        initial_R_list = initial_R_for_gpa
      )
      # Map results back to the full list of subjects
      R_final_list_full <- vector("list", num_subjects)
      gpa_result_idx <- 1
      for(subj_idx in 1:num_subjects){
          if(valid_A_indices_for_gpa[subj_idx]){
              R_final_list_full[[subj_idx]] <- geo_gpa_results$R_final_list[[gpa_result_idx]]
              gpa_result_idx <- gpa_result_idx + 1
          } else {
              R_final_list_full[[subj_idx]] <- matrix(NA, spectral_rank_k, spectral_rank_k)
          }
      }
      geo_gpa_results$R_final_list <- R_final_list_full
  } else { # Should not happen if length(A_originals_for_gpa) < 2 handled above
       geo_gpa_results <- list(R_final_list = replicate(num_subjects, matrix(NA, spectral_rank_k, spectral_rank_k), simplify=FALSE),
                               T_anchor_final = matrix(NA, num_anchors, spectral_rank_k), R_bar_final = diag(spectral_rank_k))
  }
  
  R_final_list <- geo_gpa_results$R_final_list
  T_anchor_final <- geo_gpa_results$T_anchor_final # This is T_anchor_geo
  R_bar_final <- geo_gpa_results$R_bar_final
  
  if (verbose) message_stage("Geo-HATSA Stage 2 complete.", interactive_only = TRUE)

  # --- Stage 3: Applying final rotations (same as run_hatsa_core) ---
  if (verbose) message_stage("Geo-HATSA Stage 3: Applying final rotations...", interactive_only = TRUE)
  U_aligned_list <- vector("list", num_subjects)
  if (num_subjects > 0) {
    for (i in 1:num_subjects) {
      U_orig_i <- U_original_list[[i]]
      R_final_i <- R_final_list[[i]]
      
      if (!is.null(U_orig_i) && is.matrix(U_orig_i) && !any(is.na(U_orig_i)) && 
          nrow(U_orig_i) == V_p && ncol(U_orig_i) == spectral_rank_k &&
          !is.null(R_final_i) && is.matrix(R_final_i) && !any(is.na(R_final_i)) &&
          nrow(R_final_i) == spectral_rank_k && ncol(R_final_i) == spectral_rank_k) {
        U_aligned_list[[i]] <- U_orig_i %*% R_final_i
      } else { 
        U_aligned_list[[i]] <- matrix(NA, nrow = V_p, ncol = spectral_rank_k) 
        if(verbose && (is.null(U_orig_i) || any(is.na(U_orig_i)) || is.null(R_final_i) || any(is.na(R_final_i)))){
            message(sprintf("  Geo-HATSA Stage 3: Subject %d has NA/NULL U_original or R_final. Aligned sketch set to NA.", i))
        }
      }
    }
  }
  if (verbose) message_stage("Geo-HATSA Stage 3 complete. Core Geo-HATSA finished.", interactive_only = TRUE)
  
  # Consolidate results
  geo_hatsa_core_results <- list(
    U_aligned_list = U_aligned_list,
    R_final_list = R_final_list,
    U_original_list = U_original_list,
    Lambda_original_list = Lambda_original_list, 
    Lambda_original_gaps_list = Lambda_original_gaps_list,
    T_anchor_final = T_anchor_final,
    R_bar_final = R_bar_final, 
    qc_metrics = qc_metrics
  )

  # Parameters list (similar to run_hatsa_core, maybe add a geo_hatsa_specific field)
  # parameters <- list(
  #   k = spectral_rank_k,
  #   N_subjects = num_subjects,
  #   V_p = V_p, 
  #   method = "geo_hatsa_core", # Indicate Geo-HATSA was used
  #   anchor_indices = unique_anchor_indices,
  #   k_conn_pos = k_conn_pos,
  #   k_conn_neg = k_conn_neg,
  #   n_refine = n_refine,
  #   use_dtw = use_dtw,
  #   gpa_tol = gpa_tol 
  #   # frechet_mean_options could also be stored if needed for reproducibility
  # )
  
  # The projector expects results from `hatsa_core_results` and parameters.
  # If we use the same `hatsa_projector`, `R_bar_final` will be an extra item.
  # Or, one could define a `geo_hatsa_projector` if R_bar_final needs special handling.
  # For now, just returning the list. The projector constructor would need to know about it.

  return(geo_hatsa_core_results)
} 