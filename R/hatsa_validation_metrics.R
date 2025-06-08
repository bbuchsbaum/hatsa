#' @file hatsa_validation_metrics.R
#' @title Validation Metrics for HATSA
#' @description Functions to compute various validation metrics for assessing
#' HATSA performance, typically using toy data with known ground truth.
#'
#' @importFrom stats cor
#' @importFrom vegan procrustes
NULL

#' Compute Group Template (v) Recovery Metrics
#'
#' Evaluates how well the estimated group-level spectral template (`object$v`
#' from a `hatsa_projector` object) matches a true ground-truth spectral basis.
#'
#' @param hatsa_object A fitted `hatsa_projector` object.
#' @param U_true The ground-truth spectral basis (matrix, e.g., Vp x k).
#' @param ... Additional arguments passed to `vegan::procrustes`.
#'
#' @return A list containing:
#'   \item{correlation}{Pearson correlation between the vectorized aligned estimated
#'     template and the vectorized true template.}
#'   \item{frobenius_norm_diff}{Frobenius norm of the difference between the
#'     aligned estimated template and the true template.}
#'   \item{procrustes_result}{The result of the `vegan::procrustes` call.}
#'   \item{v_aligned}{The Procrustes-aligned estimated group template.}
#'
#' @export
#' @examples
#' # See core-hatsa-toy-example.Rmd for usage.
compute_v_recovery <- function(hatsa_object, U_true, ...) {
  if (!inherits(hatsa_object, "hatsa_projector")) {
    stop("`hatsa_object` must be of class 'hatsa_projector'.")
  }
  if (!is.matrix(U_true) || !is.numeric(U_true)) {
    stop("`U_true` must be a numeric matrix.")
  }
  if (is.null(hatsa_object$v)) {
      stop("`hatsa_object$v` is NULL. Cannot compute recovery.")
  }
  if (nrow(hatsa_object$v) != nrow(U_true) || ncol(hatsa_object$v) != ncol(U_true)) {
    stop("Dimensions of `hatsa_object$v` and `U_true` must match.")
  }

  v_estimated <- hatsa_object$v

  # Align estimated v to U_true
  procr_fit <- vegan::procrustes(X = U_true, Y = v_estimated, symmetric = FALSE, ...)
  v_aligned <- procr_fit$Yrot

  # Metrics
  correlation <- stats::cor(as.vector(v_aligned), as.vector(U_true))
  frobenius_norm_diff <- norm(v_aligned - U_true, type = "F")

  return(list(
    correlation = correlation,
    frobenius_norm_diff = frobenius_norm_diff,
    procrustes_result = procr_fit,
    v_aligned = v_aligned
  ))
}

#' Compute Anchor Template (T_anchor_final) Recovery Metrics
#'
#' Evaluates how well the estimated group anchor template (`object$T_anchor_final`
#' from a `hatsa_projector` object) matches the true spectral basis at the
#' anchor locations.
#'
#' @param hatsa_object A fitted `hatsa_projector` object.
#' @param U_true The ground-truth spectral basis (matrix, Vp x k) for all parcels.
#' @param anchor_indices_true A numeric vector of indices specifying which rows of
#'   `U_true` correspond to the anchor parcels.
#' @param ... Additional arguments passed to `vegan::procrustes`.
#'
#' @return A list containing:
#'   \item{correlation}{Pearson correlation between the vectorized aligned estimated
#'     anchor template and the vectorized true anchor template.}
#'   \item{frobenius_norm_diff}{Frobenius norm of the difference between the
#'     aligned estimated anchor template and the true anchor template.}
#'   \item{procrustes_result}{The result of the `vegan::procrustes` call.}
#'   \item{T_anchor_aligned}{The Procrustes-aligned estimated anchor template.}
#'
#' @export
#' @examples
#' # See core-hatsa-toy-example.Rmd for usage.
compute_anchor_template_recovery <- function(hatsa_object, U_true, anchor_indices_true, ...) {
  if (!inherits(hatsa_object, "hatsa_projector")) {
    stop("`hatsa_object` must be of class 'hatsa_projector'.")
  }
  if (!is.matrix(U_true) || !is.numeric(U_true)) {
    stop("`U_true` must be a numeric matrix.")
  }
  if (!is.numeric(anchor_indices_true) || !is.vector(anchor_indices_true)) {
    stop("`anchor_indices_true` must be a numeric vector.")
  }
  if (any(anchor_indices_true > nrow(U_true)) || any(anchor_indices_true <= 0)) {
    stop("Invalid `anchor_indices_true`.")
  }
  if (is.null(hatsa_object$T_anchor_final)) {
    stop("`hatsa_object$T_anchor_final` is NULL. Cannot compute recovery.")
  }

  T_anchor_estimated <- hatsa_object$T_anchor_final
  U_true_anchors <- U_true[anchor_indices_true, , drop = FALSE]

  if (nrow(T_anchor_estimated) != nrow(U_true_anchors) || ncol(T_anchor_estimated) != ncol(U_true_anchors)) {
    stop(paste(
      "Dimensions of `hatsa_object$T_anchor_final` (", 
      nrow(T_anchor_estimated), "x", ncol(T_anchor_estimated),
      ") and `U_true` at anchor_indices_true (",
      nrow(U_true_anchors), "x", ncol(U_true_anchors),
      ") must match.", sep=""
    ))
  }

  # Align estimated T_anchor_final to U_true_anchors
  procr_fit <- vegan::procrustes(X = U_true_anchors, Y = T_anchor_estimated, symmetric = FALSE, ...)
  T_anchor_aligned <- procr_fit$Yrot

  # Metrics
  correlation <- stats::cor(as.vector(T_anchor_aligned), as.vector(U_true_anchors))
  frobenius_norm_diff <- norm(T_anchor_aligned - U_true_anchors, type = "F")

  return(list(
    correlation = correlation,
    frobenius_norm_diff = frobenius_norm_diff,
    procrustes_result = procr_fit,
    T_anchor_aligned = T_anchor_aligned
  ))
}

#' Compute Rotation Recovery (SO(k) Alignment) Metrics
#'
#' Assesses how accurately the subject-specific rotation matrices (`R_i`)
#' estimated by HATSA match the true rotations.
#'
#' @param hatsa_object A fitted `hatsa_projector` object, which contains
#'   `R_final_list` (list of estimated k x k rotation matrices).
#' @param R_true_list A list of true k x k rotation matrices, corresponding to
#'   each subject in `hatsa_object`.
#' @param ... Additional arguments passed to `misalign_deg` function.
#'
#' @return A numeric vector of misalignment angles in degrees, one for each
#'   subject. Returns `NA` for subjects where either the estimated or true
#'   rotation is missing or invalid for comparison.
#'
#' @export
#' @examples
#' # See core-hatsa-toy-example.Rmd for usage.
compute_rotation_recovery <- function(hatsa_object, R_true_list, ...) {
  if (!inherits(hatsa_object, "hatsa_projector")) {
    stop("`hatsa_object` must be of class 'hatsa_projector'.")
  }
  if (is.null(hatsa_object$R_final_list)) {
    stop("`hatsa_object$R_final_list` is NULL. Cannot compute rotation recovery.")
  }
  if (!is.list(R_true_list)) {
    stop("`R_true_list` must be a list of matrices.")
  }

  # Replace any NULL elements in R_final_list with NA placeholder
  R_est_list <- hatsa_object$R_final_list
  for (i in seq_along(R_est_list)) {
    if (is.null(R_est_list[[i]])) {
      R_est_list[[i]] <- NA
    }
  }
  
  # Find common length for comparison
  num_subjects <- min(length(R_est_list), length(R_true_list))
  
  if (num_subjects == 0) {
    # No subjects to compare
    warning("No subjects available for rotation recovery comparison.")
    return(numeric(0))
  }
  
  # Match lengths - truncate longer list
  R_est_list <- R_est_list[1:num_subjects]
  R_true_list <- R_true_list[1:num_subjects]
  
  misalignment_degrees <- numeric(num_subjects)

  for (i in seq_len(num_subjects)) {
    R_est <- R_est_list[[i]]
    R_true <- R_true_list[[i]]

    if (length(R_est) == 1 && is.na(R_est) || is.null(R_est) || is.null(R_true) || 
        !is.matrix(R_est) || !is.matrix(R_true)) {
      misalignment_degrees[i] <- NA_real_
      warning(sprintf("Missing or non-matrix rotation for subject %d. Setting misalignment to NA.", i), call. = FALSE)
      next
    }
    
    if(!all(dim(R_est) == dim(R_true))) {
      misalignment_degrees[i] <- NA_real_
      warning(sprintf("Dimension mismatch for rotation matrices for subject %d (Est: %s, True: %s). Setting misalignment to NA.", 
                      i, paste(dim(R_est), collapse="x"), paste(dim(R_true), collapse="x")), call. = FALSE)
      next
    }
    
    # Use misalign_deg with namespace to avoid export dependency
    misalignment_degrees[i] <- hatsa::misalign_deg(R_est = R_est, R_true = R_true, ...)
  }

  return(misalignment_degrees)
}

#' Compute Eigenvalue Spectrum Fidelity Metrics
#'
#' Compares subject-specific eigenvalues from HATSA (`Lambda_original_list`)
#' with true eigenvalues (e.g., from a toy data generator).
#'
#' @param hatsa_object A fitted `hatsa_projector` object, which contains
#'   `Lambda_original_list` (list of numeric vectors of eigenvalues).
#' @param true_eigenvalues_list A list of numeric vectors, where each vector contains
#'   the true eigenvalues for the corresponding subject. If a single vector is provided,
#'   it is assumed to be the common true eigenvalues for all subjects.
#' @param k_to_compare Integer or NULL. The number of top eigenvalues to compare.
#'   If NULL (default), all available eigenvalues are compared (up to the length of
#'   the shorter of the estimated or true eigenvalue vectors for that subject).
#'
#' @return A list of lists, one for each subject. Each inner list contains:
#'   \item{correlation}{Pearson correlation between estimated and true eigenvalues.}
#'   \item{mse}{Mean Squared Error between estimated and true eigenvalues.}
#'   \item{num_compared}{The number of eigenvalue pairs actually compared.}
#'   Returns NULL for a subject if inputs are invalid for that subject.
#'
#' @export
#' @examples
#' # See core-hatsa-toy-example.Rmd for usage.
compute_eigenvalue_fidelity <- function(hatsa_object, true_eigenvalues_list, k_to_compare = NULL) {
  if (!inherits(hatsa_object, "hatsa_projector")) {
    stop("`hatsa_object` must be of class 'hatsa_projector'.")
  }
  if (is.null(hatsa_object$Lambda_original_list)) {
    stop("`hatsa_object$Lambda_original_list` is NULL. Cannot compute eigenvalue fidelity.")
  }
  if (!is.list(hatsa_object$Lambda_original_list)) {
    stop("`hatsa_object$Lambda_original_list` must be a list.")
  }
  if (!is.list(true_eigenvalues_list) && !is.numeric(true_eigenvalues_list)){
    stop("`true_eigenvalues_list` must be a list of numeric vectors or a single numeric vector.")
  }
  
  # If true_eigenvalues_list is a single vector, replicate it for all subjects
  if (is.numeric(true_eigenvalues_list) && !is.list(true_eigenvalues_list)) {
    true_eigenvalues_list <- rep(list(true_eigenvalues_list), length(hatsa_object$Lambda_original_list))
  }

  if (length(hatsa_object$Lambda_original_list) != length(true_eigenvalues_list)) {
    stop("Length of `hatsa_object$Lambda_original_list` and `true_eigenvalues_list` must be equal.")
  }
  if (!is.null(k_to_compare) && (!is.numeric(k_to_compare) || length(k_to_compare) != 1 || k_to_compare < 1)){
    stop("`k_to_compare` must be a single positive integer or NULL.")
  }

  num_subjects <- length(hatsa_object$Lambda_original_list)
  results_list <- vector("list", num_subjects)
  names(results_list) <- names(hatsa_object$Lambda_original_list) # Preserve names if any
  if (is.null(names(results_list)) && num_subjects > 0) {
    # Use subject index if no names provided
    names(results_list) <- if(!is.null(names(true_eigenvalues_list))) {
      names(true_eigenvalues_list)
    } else {
      paste0("Subject_", seq_len(num_subjects))
    }
  }

  for (i in seq_len(num_subjects)) {
    est_lambdas <- hatsa_object$Lambda_original_list[[i]]
    true_lambdas <- true_eigenvalues_list[[i]]

    if (is.null(est_lambdas) || !is.numeric(est_lambdas) || 
        is.null(true_lambdas) || !is.numeric(true_lambdas)) {
      warning(sprintf("Missing or non-numeric eigenvalues for subject %d. Skipping.", i), call. = FALSE)
      results_list[[i]] <- list(correlation = NA_real_, mse = NA_real_, num_compared = 0)
      next
    }
    
    len_est <- length(est_lambdas)
    len_true <- length(true_lambdas)
    
    if (len_est == 0 || len_true == 0) {
        warning(sprintf("Empty eigenvalue vector(s) for subject %d. Skipping.", i), call. = FALSE)
        results_list[[i]] <- list(correlation = NA_real_, mse = NA_real_, num_compared = 0)
        next
    }

    current_k <- min(len_est, len_true)
    if (!is.null(k_to_compare)) {
      if (k_to_compare > current_k) {
        warning(sprintf("k_to_compare (%d) is greater than available eigenvalues for subject %d (min_len=%d). Using min_len.", 
                        k_to_compare, i, current_k), call. = FALSE)
      } else {
        current_k <- k_to_compare
      }
    }
    
    if (current_k == 0) { # Should be caught by len_est/len_true checks but as safeguard
        results_list[[i]] <- list(correlation = NA_real_, mse = NA_real_, num_compared = 0)
        next
    }

    est_lambdas_comp <- est_lambdas[1:current_k]
    true_lambdas_comp <- true_lambdas[1:current_k]
    
    # Ensure no NAs that would break cor or mean, replace with warning and NA result
    if(anyNA(est_lambdas_comp) || anyNA(true_lambdas_comp)){
        warning(sprintf("NA values found in eigenvalues for comparison for subject %d. Metrics will be NA.", i), call. = FALSE)
        results_list[[i]] <- list(correlation = NA_real_, mse = NA_real_, num_compared = current_k)
        next
    }
    
    # Avoid correlation with zero variance vectors
    cor_val <- NA_real_
    if (stats::var(est_lambdas_comp) > 1e-9 && stats::var(true_lambdas_comp) > 1e-9) {
        cor_val <- tryCatch(stats::cor(est_lambdas_comp, true_lambdas_comp),
                              error = function(e) NA_real_)
    } else if (identical(est_lambdas_comp, true_lambdas_comp)) {
        cor_val <- 1.0 # Perfect correlation if both are constant and identical
    }

    mse_val <- mean((est_lambdas_comp - true_lambdas_comp)^2)

    results_list[[i]] <- list(
      correlation = cor_val,
      mse = mse_val,
      num_compared = current_k
    )
  }

  return(results_list)
} 