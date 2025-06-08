# Helper function to extract a single subject's aligned sketch
# It should handle cases where U_aligned_list is present or fallback to object$s
#' @importFrom Matrix isSymmetric
NULL

.get_subject_aligned_sketch <- function(object, subject_idx) {
  if (!inherits(object, "hatsa_projector")) {
    stop(".get_subject_aligned_sketch expects a hatsa_projector object")
  }
  if (subject_idx < 1 || subject_idx > object$parameters$N_subjects) {
    stop(sprintf("Invalid subject_idx %d for N_subjects = %d", subject_idx, object$parameters$N_subjects))
  }

  V_p <- object$parameters$V_p
  k <- object$parameters$k

  # Prefer U_aligned_list if available and populated for this subject
  if (!is.null(object$U_aligned_list) && length(object$U_aligned_list) >= subject_idx) {
    sketch <- object$U_aligned_list[[subject_idx]]
    if (!is.null(sketch)) {
      if (is.matrix(sketch) && ncol(sketch) == k && nrow(sketch) == V_p) {
        return(sketch)
      } else {
        warning(sprintf("U_aligned_list[[%d]] is not a valid matrix. Falling back to object$s.", subject_idx))
      }
    }
  }
  
  # Fallback to object$s and block_indices
  if (!is.null(object$s) && !is.null(object$block_indices) && 
      length(object$block_indices) >= subject_idx && 
      !is.null(object$block_indices[[subject_idx]])) {
    
    rows_subj_i <- object$block_indices[[subject_idx]]
    if (length(rows_subj_i) == 0 && V_p > 0) { # Subject might have had no data
        warning(sprintf("Subject %d has no rows in object$s (block_indices are empty). Returning NULL sketch.", subject_idx))
        return(NULL)
    } else if (length(rows_subj_i) == 0 && V_p == 0) {
        return(matrix(NA_real_, nrow=0, ncol=k)) # Consistent 0-row matrix if V_p = 0
    }
    
    if (max(rows_subj_i) > nrow(object$s)) {
        warning(sprintf("Subject %d block_indices exceed nrow(object$s). Returning NULL sketch.", subject_idx))
        return(NULL)
    }

    U_aligned_i_from_s <- object$s[rows_subj_i, , drop = FALSE]

    if (nrow(U_aligned_i_from_s) == V_p && ncol(U_aligned_i_from_s) == k) {
      if (any(is.na(U_aligned_i_from_s))) {
        warning(sprintf("Sketch for subject %d from object$s contains NAs.", subject_idx))
      }
      return(U_aligned_i_from_s)
    } else {
      warning(sprintf("Sketch for subject %d from object$s has unexpected dimensions: %d x %d (expected %d x %d). Returning NULL.",
                      subject_idx, nrow(U_aligned_i_from_s), ncol(U_aligned_i_from_s),
                      V_p, k))
      return(NULL)
    }
  } else {
    warning(sprintf("Could not extract sketch for subject %d using U_aligned_list or object$s.", subject_idx))
    return(NULL)
  }
}

#' Compute Covariance of Aligned Spectral Coefficients
#'
#' For a given subject's aligned spectral sketch `U_aligned_i` (V_p x k),
#' this function computes the `k x k` covariance matrix of these coefficients
#' across the V_p parcels.
#'
#' @param U_aligned_subject A numeric matrix (V_p x k) of aligned spectral coefficients.
#' @return A `k x k` covariance matrix. Returns a matrix of NAs if V_p <= 1 or k = 0, or if input is NULL.
#' @importFrom stats cov
#' @keywords internal
.compute_cov_spectral_coeffs <- function(U_aligned_subject) {
  if (is.null(U_aligned_subject)) {
    return(NULL)
  }
  if (!is.matrix(U_aligned_subject) || !is.numeric(U_aligned_subject)) {
    stop("U_aligned_subject must be a numeric matrix.")
  }
  V_p <- nrow(U_aligned_subject)
  k <- ncol(U_aligned_subject)

  if (k == 0) {
    return(matrix(NA_real_, nrow = 0, ncol = 0))
  }
  if (V_p <= 1) {
    return(matrix(NA_real_, nrow = k, ncol = k))
  }
  return(stats::cov(U_aligned_subject))
}

#' Densify, symmetrize and regularize a connectivity matrix
#'
#' Helper used by `get_spd_representations()` to convert sparse connectivity
#' matrices to dense form if needed, enforce symmetry, and apply SPD
#' regularization.
#'
#' @param M A numeric matrix or `Matrix` object representing a connectivity matrix.
#' @param epsilon Numeric scalar used for SPD regularization.
#' @return A dense symmetric positive-definite matrix, or `NULL` if `M` is `NULL`.
#' @keywords internal
.densify_symmetrize_regularize <- function(M, epsilon) {
  if (is.null(M)) return(NULL)

  if (!Matrix::isSymmetric(M)) {
    warning("Matrix was not symmetric. Symmetrizing via (M + t(M)) / 2.")
    M <- (M + Matrix::t(M)) / 2
  }

  dense_M <- if (inherits(M, "sparseMatrix")) {
    if (prod(dim(M)) > 1e6) {
      warning("Coercing large sparse matrix to dense may be memory-intensive.")
    }
    as.matrix(M)
  } else {
    M
  }

  .regularize_spd(dense_M, epsilon)
}


#' Get SPD Matrix Representations from a HATSA Projector Object
#'
#' Extracts or computes various types of Symmetric Positive-Definite (SPD) matrix
#' representations for subjects from a `hatsa_projector` object.
#'
#' @param object A `hatsa_projector` or `task_hatsa_projector` object.
#' @param type Character string indicating the type of SPD representation to use.
#'   Currently supported for `hatsa_projector`:
#'   `"cov_coeffs"` (covariance of aligned spectral coefficients).
#'   Other types like `"fc_conn"` might be added or supported by specific methods.
#' @param subject_idx Optional integer or vector of integers. If provided, returns
#'   SPD matrices only for these subjects. If NULL (default), for all subjects.
#' @param regularize_epsilon Small positive value for SPD regularization. Default from RGEOM-001.
#' @param subject_data_list_for_fc Optional list of subject time-series matrices
#'   (T_i x V_p), needed for `type = "fc_conn"`.
#' @param k_conn_params_for_fc List of parameters for `compute_subject_connectivity_graph_sparse`,
#'   e.g., `list(k_conn_pos = 10, k_conn_neg = 10, zscore_type = "abs", ...)`.
#'   Needed for `type = "fc_conn"`.
#' @param ... Additional arguments, potentially passed to specific computation functions.
#' @return If `subject_idx` has length 1, a single SPD matrix is returned.
#'   Otherwise a list of matrices is returned with `NULL` for subjects whose SPD
#'   representation could not be computed.
#' @export
get_spd_representations <- function(object, ...) {
  UseMethod("get_spd_representations")
}

#' @rdname get_spd_representations
#' @export
get_spd_representations.hatsa_projector <- function(object,
                                                  type = c("cov_coeffs", "fc_conn"),
                                                  subject_idx = NULL,
                                                  regularize_epsilon = 1e-6, 
                                                  subject_data_list_for_fc = NULL,
                                                  k_conn_params_for_fc = list(),
                                                  ...) {
  type <- match.arg(type)
  N_subjects <- object$parameters$N_subjects
  
  if (is.null(subject_idx)) {
    subject_indices_to_process <- 1:N_subjects
  } else {
    if (!is.numeric(subject_idx) || !all(subject_idx >= 1 & subject_idx <= N_subjects)) {
      stop(sprintf("Invalid subject_idx. Must be NULL or integers between 1 and %d.", N_subjects))
    }
    subject_indices_to_process <- as.integer(unique(subject_idx))
  }
  
  if (length(subject_indices_to_process) == 0 && N_subjects > 0) {
    return(list())
  } else if (N_subjects == 0) {
    return(list())
  }

  all_spd_matrices <- vector("list", length(subject_indices_to_process))
  names(all_spd_matrices) <- as.character(subject_indices_to_process) 

  for (i in seq_along(subject_indices_to_process)) {
    current_subj_orig_idx <- subject_indices_to_process[i]
    spd_matrix_i <- NULL 

    if (type == "cov_coeffs") {
      U_aligned_i <- .get_subject_aligned_sketch(object, current_subj_orig_idx)
      if (!is.null(U_aligned_i)) {
        cov_coeffs_i <- .compute_cov_spectral_coeffs(U_aligned_i)
        if (!is.null(cov_coeffs_i) && is.matrix(cov_coeffs_i) && all(dim(cov_coeffs_i) > 0) && !all(is.na(cov_coeffs_i))) {
            spd_matrix_i <- .regularize_spd(cov_coeffs_i, regularize_epsilon)
        } else {
            warning(sprintf("Covariance of coefficients for subject %d resulted in NULL, NA, or zero-dim matrix.", current_subj_orig_idx))
        }
      } else {
        warning(sprintf("Could not get aligned sketch for subject %d for cov_coeffs.", current_subj_orig_idx))
      }
    } else if (type == "fc_conn") {
      if (is.null(subject_data_list_for_fc)) {
        stop("`subject_data_list_for_fc` is required for type = 'fc_conn'.")
      }
      if (length(subject_data_list_for_fc) < current_subj_orig_idx || is.null(subject_data_list_for_fc[[current_subj_orig_idx]])) {
        warning(sprintf("No data in `subject_data_list_for_fc` for subject index %d.", current_subj_orig_idx))
        all_spd_matrices[[i]] <- NULL
        next
      }
      subj_ts_data <- subject_data_list_for_fc[[current_subj_orig_idx]]
      
      conn_args <- list(
        X_subject = subj_ts_data,
        k_conn_pos = ifelse(is.null(k_conn_params_for_fc$k_conn_pos), 10, k_conn_params_for_fc$k_conn_pos),
        k_conn_neg = ifelse(is.null(k_conn_params_for_fc$k_conn_neg), 10, k_conn_params_for_fc$k_conn_neg),
        zscore_type = ifelse(is.null(k_conn_params_for_fc$zscore_type), "abs", k_conn_params_for_fc$zscore_type),
        corr_method = ifelse(is.null(k_conn_params_for_fc$corr_method), "pearson", k_conn_params_for_fc$corr_method),
        cache_dir = NULL, 
        verbose = FALSE   
      )
      extra_conn_args <- k_conn_params_for_fc[!names(k_conn_params_for_fc) %in% names(conn_args)]
      final_conn_args <- c(conn_args, extra_conn_args)
      
      W_conn_i_list <- tryCatch(
        do.call(hatsa::compute_subject_connectivity_graph_sparse, final_conn_args),
        error = function(e) {
          warning(sprintf("Error computing W_conn for subject %d: %s", current_subj_orig_idx, e$message))
          NULL
        }
      )
      
      if (!is.null(W_conn_i_list) && !is.null(W_conn_i_list$W_sparse)) {
        W_conn_i <- W_conn_i_list$W_sparse
        spd_matrix_i <- .densify_symmetrize_regularize(W_conn_i, regularize_epsilon)
      } else {
        warning(sprintf("W_conn_i computation failed or returned NULL for subject %d.", current_subj_orig_idx))
      }
    } else {
      warning(sprintf("Type '%s' not yet fully implemented for get_spd_representations.hatsa_projector.", type))
    }
    all_spd_matrices[[i]] <- spd_matrix_i
  }
  
  if (is.numeric(subject_idx) && length(subject_idx) == 1 && length(all_spd_matrices) == 1) {
      return(all_spd_matrices[[1]])
  }
  return(all_spd_matrices)
}

#' @rdname get_spd_representations
#' @export
get_spd_representations.task_hatsa_projector <- function(object,
                                                       type = c("cov_coeffs", "fc_conn", "fc_task", "fc_hybrid"),
                                                       subject_idx = NULL,
                                                       regularize_epsilon = 1e-6, 
                                                       subject_data_list_for_fc = NULL,
                                                       k_conn_params_for_fc = list(),
                                                       lambda_blend_for_hybrid = NULL, 
                                                       ...) {
  type <- match.arg(type)
  N_subjects <- object$parameters$N_subjects

  if (is.null(subject_idx)) {
    subject_indices_to_process <- 1:N_subjects
  } else {
    if (!is.numeric(subject_idx) || !all(subject_idx >= 1 & subject_idx <= N_subjects)) {
      stop(sprintf("Invalid subject_idx. Must be NULL or integers between 1 and %d.", N_subjects))
    }
    subject_indices_to_process <- as.integer(unique(subject_idx))
  }

  if (length(subject_indices_to_process) == 0 && N_subjects > 0) return(list())
  if (N_subjects == 0) return(list())

  all_spd_matrices <- vector("list", length(subject_indices_to_process))
  names(all_spd_matrices) <- as.character(subject_indices_to_process)

  if (type %in% c("cov_coeffs", "fc_conn")) {
    spd_matrices_from_parent <- get_spd_representations.hatsa_projector(
        object = object, 
        type = type, 
        subject_idx = subject_indices_to_process, 
        regularize_epsilon = regularize_epsilon,
        subject_data_list_for_fc = subject_data_list_for_fc,
        k_conn_params_for_fc = k_conn_params_for_fc,
        ...
    )
    if (is.numeric(subject_idx) && length(subject_idx) == 1 && !is.list(spd_matrices_from_parent)){
        all_spd_matrices[[1]] <- spd_matrices_from_parent
    } else {
        all_spd_matrices <- spd_matrices_from_parent 
    }
    
  } else if (type %in% c("fc_task", "fc_hybrid")) {
    for (i in seq_along(subject_indices_to_process)) {
      current_subj_orig_idx <- subject_indices_to_process[i]
      spd_matrix_i <- NULL
      W_matrix <- NULL

      if (type == "fc_task") {
        if (!is.null(object$W_task_list) && length(object$W_task_list) >= current_subj_orig_idx) {
          W_matrix <- object$W_task_list[[current_subj_orig_idx]]
          if(is.null(W_matrix)) warning(sprintf("Stored W_task_list for subject %d is NULL.", current_subj_orig_idx))
        } else {
          warning(sprintf("W_task_list not available or not populated for subject %d.", current_subj_orig_idx))
        }
      } else if (type == "fc_hybrid") {
        if (!is.null(object$W_hybrid_list) && length(object$W_hybrid_list) >= current_subj_orig_idx) {
          W_matrix <- object$W_hybrid_list[[current_subj_orig_idx]]
           if(is.null(W_matrix)) warning(sprintf("Stored W_hybrid_list for subject %d is NULL.", current_subj_orig_idx))
        } else {
          warning(sprintf("W_hybrid_list not available or not populated for subject %d.", current_subj_orig_idx))
        }
      }

      if (!is.null(W_matrix)) {
        spd_matrix_i <- .densify_symmetrize_regularize(W_matrix, regularize_epsilon)
        if (is.null(spd_matrix_i)) {
          warning(sprintf("%s for subject %d resulted in NULL after processing.", type, current_subj_orig_idx))
        }
      } else {
        warning(sprintf("Could not retrieve or compute %s matrix for subject %d.", type, current_subj_orig_idx))
      }
      all_spd_matrices[[as.character(current_subj_orig_idx)]] <- spd_matrix_i 
    }
  } else {
    stop(sprintf("Unsupported type '%s' in get_spd_representations.task_hatsa_projector.", type))
  }
  
  if (is.numeric(subject_idx) && length(subject_idx) == 1 && length(all_spd_matrices) == 1) {
      return(all_spd_matrices[[1]])
  }
  return(all_spd_matrices)
} 