#' Validate inputs for run_hatsa_core
#'
#' Checks the validity of input parameters for the main HATSA function.
#' Stops with an error if validation fails.
#'
#' @param subject_data_list A list of matrices, where each matrix `X_i` is
#'   `T_i x V_p` (time points by parcels) for subject `i`.
#' @param anchor_indices A numeric vector of 1-based indices for anchor parcels.
#' @param spectral_rank_k An integer, the spectral rank `k`. Must be `>= 0`.
#' @param k_conn_pos An integer, number of positive connections for sparsification.
#' @param k_conn_neg An integer, number of negative connections for sparsification.
#' @param n_refine An integer, number of GPA refinement iterations.
#' @param V_p An integer, the number of parcels (columns in `subject_data_list` elements).
#'
#' @return Invisibly returns a list containing potentially modified `anchor_indices`
#'   (uniquified) if validation passes.
#' @keywords internal
validate_hatsa_inputs <- function(subject_data_list, anchor_indices,
                                  spectral_rank_k, k_conn_pos, k_conn_neg,
                                  n_refine, V_p) {
  # (Code mostly unchanged from v0.2.0, but with k > m as error)
  if (!is.list(subject_data_list) || length(subject_data_list) == 0) {
    stop("`subject_data_list` must be a non-empty list.")
  }
  if (!all(sapply(subject_data_list, function(m) inherits(m, "matrix") && !inherits(m, "Matrix")))) {
    stop("All elements of `subject_data_list` must be standard R matrices (not sparse Matrix objects).")
  }
  if (V_p > 0 && !all(sapply(subject_data_list, ncol) == V_p)) {
    stop("All matrices in `subject_data_list` must have the same number of columns (parcels).")
  }
  if (any(sapply(subject_data_list, function(X) any(!is.finite(X))))) {
    stop("All matrices in `subject_data_list` must contain finite values.")
  }

  if (!is.numeric(anchor_indices) || !is.vector(anchor_indices) || any(is.na(anchor_indices))) {
      stop("`anchor_indices` must be a numeric vector without NAs.")
  }
  if (V_p > 0 && (any(anchor_indices <= 0) || any(anchor_indices > V_p))) {
    stop("`anchor_indices` must be valid 1-based parcel indices not exceeding V_p.")
  }
  
  unique_anchor_indices <- unique(anchor_indices)
  if (length(unique_anchor_indices) != length(anchor_indices)) {
      if (interactive()) message("Duplicate anchor indices provided; using unique set of anchors.")
  }
  if (length(unique_anchor_indices) == 0 && V_p > 0) {
    stop("`anchor_indices` (after taking unique) must not be empty if V_p > 0.")
  }
  m <- length(unique_anchor_indices)

  if (!is.numeric(spectral_rank_k) || length(spectral_rank_k) != 1 || spectral_rank_k < 0 || spectral_rank_k != round(spectral_rank_k)) {
    stop("`spectral_rank_k` must be a non-negative integer.")
  }
  if (V_p > 0 && spectral_rank_k == 0 && interactive()) {
      message("`spectral_rank_k` is 0. Output sketches will have 0 columns. Procrustes alignment will be trivial.")
  }
  if (V_p > 0 && spectral_rank_k > V_p) {
    stop("`spectral_rank_k` cannot exceed the number of parcels `V_p`.")
  }
  if (m > 0 && spectral_rank_k > 0 && spectral_rank_k > m) { # k=0 is fine if m=0 or m>0
    stop(sprintf("`spectral_rank_k` (%d) cannot be greater than the number of unique anchors `m` (%d) for stable Procrustes. Anchor matrix would be rank-deficient.", spectral_rank_k, m))
  }

  if (!is.numeric(k_conn_pos) || length(k_conn_pos) != 1 || k_conn_pos < 0 || k_conn_pos != round(k_conn_pos)) {
    stop("`k_conn_pos` must be a non-negative integer.")
  }
  if (!is.numeric(k_conn_neg) || length(k_conn_neg) != 1 || k_conn_neg < 0 || k_conn_neg != round(k_conn_neg)) {
    stop("`k_conn_neg` must be a non-negative integer.")
  }
  if (V_p > 1 && (k_conn_pos + k_conn_neg == 0) && interactive()) {
      message("`k_conn_pos` and `k_conn_neg` are both zero. Connectivity graphs will be empty.")
  }
   if (V_p > 1 && (k_conn_pos + k_conn_neg >= (V_p -1) ) ) { 
    if (interactive()) message(sprintf("Warning: `k_conn_pos` + `k_conn_neg` (%d) is high relative to V_p-1 (%d), may lead to dense graphs.",
                    k_conn_pos + k_conn_neg, V_p-1))
  }

  if (!is.numeric(n_refine) || length(n_refine) != 1 || n_refine < 0 || n_refine != round(n_refine)) {
    stop("`n_refine` must be a non-negative integer.")
  }

  return(invisible(list(unique_anchor_indices = unique_anchor_indices)))
}

#' Safe Z-scoring of non-zero values in a (potentially sparse) matrix
#'
#' Z-scores the non-zero elements of a matrix. If the standard deviation of these
#' non-zero values is 0 (e.g., all non-zero values are identical, or only one
#' non-zero value exists), these values are set to 0.
#' IMPORTANT: For symmetric matrices, ensure that only one triangle (and diagonal)
#' is stored in `x@x` (e.g., by `Matrix::forceSymmetric(x, uplo="U")` prior to calling)
#' to prevent breaking symmetry if `x@x` contained duplicates from both triangles.
#'
#' @param x A numeric matrix (can be a base R matrix or a sparse `Matrix` object).
#'   If sparse, it's assumed `x@x` holds the structurally non-zero values.
#' @return A matrix of the same class and dimensions as `x`, with its non-zero
#'   elements z-scored.
#' @importFrom Matrix nnzero
#' @importFrom stats sd
#' @keywords internal
zscore_nonzero_sparse <- function(x) {
  # (Code mostly unchanged from v0.2.0)
  is_sparse <- inherits(x, "Matrix")
  
  if (is_sparse) {
    nnz_x <- Matrix::nnzero(x)
    if (is.na(nnz_x) || nnz_x == 0) return(x)
    vals <- x@x 
    
    # Clean up non-finite values
    non_finite_vals <- !is.finite(vals)
    if (any(non_finite_vals)) {
      vals[non_finite_vals] <- 0
      x@x[non_finite_vals] <- 0
    }
  } else { 
    non_zero_indices <- which(x != 0)
    if (length(non_zero_indices) == 0) return(x)
    vals <- x[non_zero_indices]
    
    # Clean up non-finite values
    non_finite_vals <- !is.finite(vals)
    if (any(non_finite_vals)) {
      vals[non_finite_vals] <- 0
      x[non_zero_indices[non_finite_vals]] <- 0
    }
  }
  
  # Return early if all values are now zero after removing non-finite values
  if (all(vals == 0)) return(x)
  
  if (length(vals) < 2) {
      if (is_sparse) {
        x@x <- rep(0, length(x@x))
      } else {
        x[non_zero_indices] <- 0
      }
      return(x)
  }

  mean_val <- mean(vals, na.rm = TRUE)
  sd_val <- stats::sd(vals, na.rm = TRUE)
  
  if (is.na(sd_val) || sd_val == 0) {
    if (is_sparse) {
      x@x <- rep(0, length(x@x))
    } else {
      x[non_zero_indices] <- 0
    }
  } else {
    if (is_sparse) {
      x@x <- (vals - mean_val) / sd_val
    } else {
      x[non_zero_indices] <- (vals - mean_val) / sd_val
    }
  }
  return(x)
}


#' Helper for status messages
#'
#' Prints a timestamped stage message when \code{verbose} is TRUE. If
#' \code{interactive_only} is TRUE the message is printed only during interactive sessions.
#'
#' @param message_text Text message to display.
#' @param verbose Whether to print the message.
#' @param interactive_only Whether to print only in interactive mode.
#' @keywords internal
message_stage <- function(message_text, verbose = TRUE, interactive_only = FALSE) {
  if (!verbose) return(invisible(NULL))
  if (interactive_only && !interactive()) return(invisible(NULL))
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", message_text))
}


#' Null-coalescing operator
#'
#' Returns \code{x} if it is not \code{NULL}, otherwise \code{y}.
#'
#' @param x, y Values to compare.
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x



#' Safely evaluate an expression with warning on error
#'
#' Evaluates `expr` and returns its result. If an error occurs, a formatted
#' warning is issued and `default` is returned.
#'
#' @param expr Expression to evaluate.
#' @param msg_fmt A format string passed to `sprintf`. The caught error message
#'   will be inserted via `%s`.
#' @param default Value to return if evaluation fails. Defaults to `NULL`.
#' @return The result of `expr`, or `default` on error.
#' @keywords internal
safe_wrapper <- function(expr, msg_fmt, default = NULL) {
  expr_sub <- substitute(expr)
  tryCatch(
    eval(expr_sub, parent.frame()),
    error = function(e) {
      warning(sprintf(msg_fmt, e$message))
      default
    }
  )
}

