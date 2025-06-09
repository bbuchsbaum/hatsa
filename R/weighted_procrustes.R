#' @title Solve Weighted Procrustes Rotation Problem
#' @description Finds the optimal rotation matrix `R` that aligns a source matrix `A_source`
#'   to a target matrix `T_target`, minimizing `|| Omega (A_source R - T_target) ||_F^2`,
#'   where `Omega` is a diagonal weighting matrix.
#'
#' @param A_source The numeric source matrix (e.g., augmented anchors for one subject),
#'   dimensions `N x k` (N = total anchors, k = dimensions).
#' @param T_target The numeric target matrix (e.g., augmented template), dimensions `N x k`.
#' @param m_parcel_rows Integer, number of rows in `A_source` (and `T_target`)
#'   corresponding to parcel anchors. Assumed to be the first `m_parcel_rows`.
#' @param m_task_rows Integer, number of rows in `A_source` (and `T_target`)
#'   corresponding to task condition anchors. Assumed to be after parcel rows.
#' @param omega_mode Character string, mode for determining weights:
#'   `"fixed"` (default) or `"adaptive"`.
#' @param fixed_omega_weights List, used if `omega_mode == "fixed"`.
#'   E.g., `list(parcel = 1.0, condition = 0.5)`. If `NULL`, defaults to this.
#' @param reliability_scores Numeric vector of length `m_task_rows`, containing
#'   reliability scores (e.g., R^2_p) for each task condition. Required and used
#'   only if `omega_mode == "adaptive"`.
#' @param scale_omega_trace Logical, whether to rescale the diagonal `Omega` matrix
#'   so that `sum(diag(Omega)) == N`. Default `TRUE`.
#'
#' @details Internally computes the weighted cross-product using
#'   `crossprod(A_source, T_target * omega_diag_vector)` to avoid creating
#'   separate weighted matrices for `A_source` and `T_target`.
#'
#' @return A `k x k` rotation matrix `R`.
#'
#' @details
#' When the spectral dimension `k` is 1, the SVD-based determinant correction is
#' unnecessary. In this case the rotation is simply the sign of the scalar
#' cross-product `M`, i.e. `R <- matrix(sign(M), 1, 1)`.
#'
#' @importFrom Matrix Diagonal
#' @importFrom stats median
#' @export
#' @examples
#' N_parcels <- 5
#' N_tasks <- 3
#' N_total <- N_parcels + N_tasks
#' k_dims <- 4
#' A <- matrix(rnorm(N_total * k_dims), N_total, k_dims)
#' T_true <- matrix(rnorm(N_total * k_dims), N_total, k_dims)
#' R_true <- svd(matrix(rnorm(k_dims*k_dims),k_dims,k_dims))$u %*% 
#'           diag(sample(c(1,1,1,-1),k_dims,replace=TRUE)) %*% 
#'           svd(matrix(rnorm(k_dims*k_dims),k_dims,k_dims))$v
#' A_rotated_perfect <- T_true %*% t(R_true) # A_rotated_perfect R_true should be T_true
#'
#' # Fixed weights (default)
#' R_fixed <- solve_procrustes_rotation_weighted(A_rotated_perfect, T_true, 
#'                                               N_parcels, N_tasks)
#' # print(all.equal(T_true %*% t(R_fixed), A_rotated_perfect)) # Check alignment
#' print(sum((A_rotated_perfect %*% R_fixed - T_true)^2)) # Should be small
#'
#' # Adaptive weights (example with dummy reliabilities)
#' set.seed(123)
#' rel_scores <- runif(N_tasks, -0.2, 0.8)
#' R_adaptive <- solve_procrustes_rotation_weighted(A_rotated_perfect, T_true, 
#'                                                  N_parcels, N_tasks, 
#'                                                  omega_mode = "adaptive", 
#'                                                  reliability_scores = rel_scores)
#' print(sum((A_rotated_perfect %*% R_adaptive - T_true)^2)) # Should be small
#'
#' # No task rows
#' R_no_task <- solve_procrustes_rotation_weighted(A[1:N_parcels,], T_true[1:N_parcels,],
#'                                                 N_parcels, 0)
#'
solve_procrustes_rotation_weighted <- function(A_source, T_target, 
                                             m_parcel_rows, m_task_rows,
                                             omega_mode = "fixed",
                                             fixed_omega_weights = NULL,
                                             reliability_scores = NULL,
                                             scale_omega_trace = TRUE) {

  # --- Input Validation ---
  if (!is.matrix(A_source) || !is.numeric(A_source) || 
      !is.matrix(T_target) || !is.numeric(T_target)) {
    stop("A_source and T_target must be numeric matrices.")
  }
  if (nrow(A_source) != nrow(T_target) || ncol(A_source) != ncol(T_target)) {
    stop("A_source and T_target must have the same dimensions.")
  }
  N_total_rows <- nrow(A_source)
  k_dims <- ncol(A_source)
  if (N_total_rows == 0) { 
      # If k_dims is also 0, diag(0) is matrix(0,0). If k_dims > 0, return KxK identity.
      return(diag(k_dims)) 
  }
  if (m_parcel_rows < 0 || m_task_rows < 0) {
      stop("m_parcel_rows and m_task_rows must be non-negative.")
  }
  if (m_parcel_rows + m_task_rows != N_total_rows) {
    stop("Sum of m_parcel_rows and m_task_rows must equal total rows in A_source/T_target.")
  }
  omega_mode <- match.arg(omega_mode, c("fixed", "adaptive"))

  if (is.null(fixed_omega_weights)) {
    fixed_omega_weights <- list(parcel = 1.0, condition = 0.5)
  }
  if (!is.list(fixed_omega_weights) || 
      !all(c("parcel", "condition") %in% names(fixed_omega_weights)) || 
      !is.numeric(fixed_omega_weights$parcel) || length(fixed_omega_weights$parcel)!=1 ||
      !is.numeric(fixed_omega_weights$condition) || length(fixed_omega_weights$condition)!=1 ){
      stop("fixed_omega_weights must be a list with numeric elements 'parcel' and 'condition'.")
  }

  if (omega_mode == "adaptive") {
    if (m_task_rows > 0 && (is.null(reliability_scores) || !is.numeric(reliability_scores) || length(reliability_scores) != m_task_rows)) {
      stop("If omega_mode is 'adaptive' and m_task_rows > 0, reliability_scores must be a numeric vector of length m_task_rows.")
    }
    if (m_task_rows == 0 && !is.null(reliability_scores)){
        warning("reliability_scores provided but m_task_rows is 0; scores will be ignored.")
    }
  }
  
  # --- Construct Omega diagonal vector (using audit patch logic) ---
  omega_diag_vector <- numeric(N_total_rows)

  if (m_parcel_rows > 0) { # Use > 0 check for safety
      omega_diag_vector[1:m_parcel_rows] <- fixed_omega_weights$parcel
  }
  if (m_task_rows > 0) { # Use > 0 check
      idx <- (m_parcel_rows + 1):N_total_rows
      if (omega_mode == "fixed") {
          omega_diag_vector[idx] <- fixed_omega_weights$condition
      } else { # adaptive
          # Safeguard reliability scores: ensure finite, treat NA/NaN/Inf as 0
          safe_rel_scores <- reliability_scores
          non_finite_idx <- !is.finite(safe_rel_scores)
          if(any(non_finite_idx)){
              safe_rel_scores[non_finite_idx] <- 0
              warning("Non-finite reliability_scores found; treated as 0 for weighting.")
          }
          
          safe_r <- pmax(0, safe_rel_scores) # Ensure non-negative
          med <- stats::median(safe_r[safe_r > 0], na.rm = TRUE) # Median of strictly positive scores
          if (is.na(med)) med <- 0 # Handle case where no scores > 0
          med <- max(med, 1e-2) # Audit clamp: Avoid division by tiny median
          
          w_raw <- fixed_omega_weights$condition * safe_r / med
          omega_diag_vector[idx] <- pmin(w_raw, fixed_omega_weights$parcel) # Cap at parcel weight
      }
  }

  # Audit check: Early return if all weights are zero
  if (all(abs(omega_diag_vector) < 1e-14)) {
      warning("All Omega weights are effectively zero. Rotation is undefined; returning identity matrix.")
      return(diag(k_dims))
  }

  # --- Trace Rescaling (Optional) ---
  if (scale_omega_trace && N_total_rows > 0) {
      current_trace <- sum(omega_diag_vector)
      if (abs(current_trace) > 1e-9) { # Check absolute value for safety
          rescale_factor <- N_total_rows / current_trace
          omega_diag_vector <- omega_diag_vector * rescale_factor
      } else { 
          # If trace is near zero, but not all weights were zero initially (handled above),
          # it implies weights canceled out or were extremely small. Unweighted seems safest.
          warning("Sum of omega weights is near zero after potential adaptive step; proceeding as unweighted Procrustes.")
          omega_diag_vector <- rep(1.0, N_total_rows) 
      }
  }
  
  # --- Weighted cross-product without duplicating A_source ---
  # Compute t(A_source) %*% diag(omega) %*% T_target by scaling T_target rows
  # directly. This avoids creating intermediate weighted copies of both matrices.

  # --- Solve Standard Procrustes using SVD (Patch logic Fix 2-C) ---
  if (k_dims == 0) {
    return(matrix(0,0,0))
  }

  M <- crossprod(A_source, T_target * omega_diag_vector) # k_dims x k_dims matrix
  if (all(abs(M) < 1e-14)) {
      warning("Cross-product matrix M (weighted) is near zero; rotation is ill-defined. Returning identity.")
      return(diag(k_dims))
  }

  if (k_dims == 1) {
      return(matrix(sign(M), 1, 1))
  }
  
  svd_M <- svd(M)
  
  # Standard Procrustes solution with proper determinant correction
  # R = U %*% t(V) ensures A_source %*% R ≈ T_target
  U <- svd_M$u
  V <- svd_M$v
  
  # Use determinant of V and U for numerically stable orientation check
  det_U <- det(U)
  det_V <- det(V)
  sign_det <- sign(det_U * det_V)
  
  # Handle numerical edge cases
  if (is.na(sign_det) || abs(sign_det) < .Machine$double.eps) {
    sign_det <- 1
  }
  
  # Correct orientation by flipping last column of U if needed
  if (sign_det < 0) {
    U[, k_dims] <- -U[, k_dims]
  }
  
  # Compute final rotation
  R <- U %*% t(V)
  
  # Re-orthogonalize for numerical stability
  # This ensures R is exactly orthogonal despite floating-point errors
  sv_R <- svd(R)
  R <- sv_R$u %*% t(sv_R$v)
  
  # Final sanity check on determinant (optional, should be close to 1 now)
  # final_det <- det(R)
  # if (abs(final_det - 1.0) > 1e-6) {
  #    warning(sprintf("Final rotation matrix determinant (%.4f) is not close to 1.", final_det))
  # }

  return(R)
} 