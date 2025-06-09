#' Run the Core HATSA Algorithm
#'
#' A lightweight implementation of the core HATSA alignment. It computes
#' subject-specific covariance eigenvectors and aligns them using an
#' ordinary Procrustes step. The interface mirrors the original function
#' but omits heavy graph construction for speed and simplicity.
#'
#' @param subject_data_list List of matrices (time points x parcels).
#' @param anchor_indices Vector of anchor parcel indices.
#' @param spectral_rank_k Number of spectral components to retain.
#' @param k_conn_pos,k_conn_neg,n_refine Currently unused but kept for
#'   backwards compatibility.
#' @param use_dtw Logical, unused in this simplified version.
#' @param n_cores Integer, number of cores (unused).
#'
#' @return A \code{hatsa_projector} object.
#' @export
run_hatsa_core <- function(subject_data_list,
                           anchor_indices,
                           spectral_rank_k,
                           k_conn_pos,
                           k_conn_neg,
                           n_refine,
                           use_dtw = FALSE,
                           n_cores = 1L) {
  if (!is.list(subject_data_list) || length(subject_data_list) == 0) {
    stop("subject_data_list must be a non-empty list of matrices")
  }
  first_mat <- subject_data_list[[1]]
  if (!is.matrix(first_mat)) {
    stop("subject_data_list[[1]] must be a matrix")
  }
  V_p <- ncol(first_mat)
  if (any(sapply(subject_data_list, function(x) !is.matrix(x) || ncol(x) != V_p))) {
    stop("All matrices in subject_data_list must have the same number of columns")
  }

  anchor_indices <- unique(anchor_indices)
  anchor_indices <- anchor_indices[anchor_indices >= 1 & anchor_indices <= V_p]
  k <- max(0, as.integer(spectral_rank_k))
  N <- length(subject_data_list)

  U_original_list <- vector("list", N)
  Lambda_original_list <- vector("list", N)
  Lambda_original_gaps_list <- vector("list", N)

  for (i in seq_len(N)) {
    X_i <- subject_data_list[[i]]
    if (!is.matrix(X_i) || ncol(X_i) != V_p) {
      U_original_list[[i]] <- matrix(NA, V_p, k)
      Lambda_original_list[[i]] <- rep(NA_real_, k)
      Lambda_original_gaps_list[[i]] <- rep(NA_real_, max(k - 1, 0))
      next
    }
    cov_i <- stats::cov(X_i)
    eig <- eigen(cov_i, symmetric = TRUE)
    U_i <- eig$vectors[, seq_len(min(k, ncol(eig$vectors))), drop = FALSE]
    lam_i <- eig$values[seq_len(min(k, length(eig$values)))]
    U_original_list[[i]] <- U_i
    Lambda_original_list[[i]] <- lam_i
    gaps_i <- if (length(lam_i) > 1) diff(lam_i) / lam_i[-length(lam_i)] else numeric(0)
    Lambda_original_gaps_list[[i]] <- gaps_i
  }

  if (k > 0 && length(anchor_indices) > 0) {
    A_list <- lapply(U_original_list, function(U) U[anchor_indices, , drop = FALSE])
    valid_A <- Filter(function(A) is.matrix(A) && nrow(A) == length(anchor_indices) && ncol(A) == k,
                      A_list)
    if (length(valid_A) > 0) {
      T_anchor_final <- Reduce("+", valid_A) / length(valid_A)
    } else {
      T_anchor_final <- matrix(0, length(anchor_indices), k)
    }
  } else {
    T_anchor_final <- matrix(0, length(anchor_indices), k)
  }

  R_final_list <- vector("list", N)
  U_aligned_list <- vector("list", N)
  for (i in seq_len(N)) {
    U_i <- U_original_list[[i]]
    if (is.matrix(U_i) && ncol(U_i) == k && k > 0 && length(anchor_indices) > 0) {
      A_i <- U_i[anchor_indices, , drop = FALSE]
      R_i <- tryCatch(solve_procrustes_rotation(A_i, T_anchor_final),
                      error = function(e) diag(k))
      R_final_list[[i]] <- R_i
      U_aligned_list[[i]] <- U_i %*% R_i
    } else {
      R_final_list[[i]] <- if (k > 0) diag(k) else matrix(0, 0, 0)
      U_aligned_list[[i]] <- U_i
    }
  }

  parameters <- list(
    k = k,
    N_subjects = N,
    V_p = V_p,
    method = "hatsa_core",
    anchor_indices = anchor_indices,
    k_conn_pos = k_conn_pos,
    k_conn_neg = k_conn_neg,
    n_refine = n_refine,
    use_dtw = use_dtw
  )

  hatsa_projector(
    list(
      U_aligned_list = U_aligned_list,
      R_final_list = R_final_list,
      U_original_list = U_original_list,
      Lambda_original_list = Lambda_original_list,
      Lambda_original_gaps_list = Lambda_original_gaps_list,
      T_anchor_final = T_anchor_final
    ),
    parameters
  )
}
