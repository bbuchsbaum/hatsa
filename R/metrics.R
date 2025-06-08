# Helper for safe matrix logarithm, returning NULL on error
# Not exported, used internally by misalign_deg
safe_logm_internal <- function(M) {
  res <- tryCatch(
    expm::logm(M),
    error = function(e) {
      # warning(paste("Matrix logarithm failed:", e$message, 
      #               ". Will use fallback for misalignment calculation."))
      return(NULL)
    }
  )
  return(res)
}

#' Calculate Misalignment Angle Between Two Rotation Matrices
#'
#' Computes the geodesic distance on SO(k) between two k x k rotation matrices,
#' representing the angle of rotation needed to align one to the other.
#' The result is given in degrees.
#'
#' @param R_est A k x k estimated rotation matrix (should be orthogonal with determinant +1).
#' @param R_true A k x k true rotation matrix (should be orthogonal with determinant +1).
#' @param method Character string, the method to use. Default is "geodesic".
#'   Currently, "geodesic" uses the matrix logarithm. 
#'   A fallback to a simpler trace-based formula (most accurate for SO(3)) is used 
#'   if `expm::logm` fails or the `expm` package is not installed.
#'
#' @return The misalignment angle in degrees. Returns 0 for identical matrices.
#'   Returns NA if dimensions are mismatched or matrices are not square.
#' 
#' @details
#' The primary method ("geodesic") calculates the angle as:
#' `||logm(R_true^T R_est)||_F / sqrt(2)`,
#' where `logm` is the matrix logarithm and `||.||_F` is the Frobenius norm.
#' This is a standard geodesic distance on the special orthogonal group SO(k).
#'
#' If the `expm` package is not available, or if `expm::logm` fails (e.g., due to
#' numerical issues if matrices are not perfectly orthogonal), the function
#' falls back to a simpler formula derived from the trace:
#' `acos((trace(R_true^T R_est) - 1) / 2)`.
#' This fallback formula is exact for SO(3) and provides a dissimilarity measure
#' for other k, being 0 for perfect alignment.
#'
#' @examples
#' if (requireNamespace("expm", quietly = TRUE)) {
#'   R1 <- diag(3)
#'   theta <- pi/4
#'   R2 <- matrix(c(cos(theta), -sin(theta), 0,
#'                  sin(theta),  cos(theta), 0,
#'                  0,           0,          1), nrow=3, byrow=TRUE)
#'   misalign_deg(R2, R1) # Should be 45 degrees
#'
#'   # Example for k=2
#'   R_true_2d <- matrix(c(1,0,0,1),2,2)
#'   angle_rad_2d <- pi/6 # 30 degrees
#'   R_est_2d <- matrix(c(cos(angle_rad_2d), -sin(angle_rad_2d),
#'                        sin(angle_rad_2d), cos(angle_rad_2d)), 2, 2)
#'   misalign_deg(R_est_2d, R_true_2d)
#' }
#' 
#' @export
#' @importFrom expm logm
misalign_deg <- function(R_est, R_true, method = "geodesic") {
  if (!is.matrix(R_est) || !is.matrix(R_true)) {
    warning("Inputs must be matrices.")
    return(NA_real_)
  }
  if (!all(dim(R_est) == dim(R_true))) {
    warning("Rotation matrices must have the same dimensions.")
    return(NA_real_)
  }
  if (nrow(R_est) != ncol(R_est)) {
    warning("Matrices must be square.")
    return(NA_real_)
  }
  k_dim <- nrow(R_est)

  # Calculate relative rotation: M_rel = R_true^T %*% R_est
  # This measures how R_est is rotated relative to R_true.
  M_rel <- crossprod(R_true, R_est) 

  # Check for perfect alignment first to avoid numerical issues with logm or acos
  if (isTRUE(all.equal(M_rel, diag(k_dim), tolerance = 1e-7))) {
    return(0)
  }

  angle_rad <- NA_real_
  use_fallback <- TRUE

  if (method == "geodesic") {
    if (requireNamespace("expm", quietly = TRUE)) {
      log_M_rel <- safe_logm_internal(M_rel)
      if (!is.null(log_M_rel)) {
        # Geodesic distance on SO(k) is ||logm(M_rel)||_F / sqrt(2)
        angle_rad <- norm(log_M_rel, type = "F") / sqrt(2)
        use_fallback <- FALSE
      } else {
        warning("Matrix logarithm failed for geodesic method. Using fallback trace-based calculation.", call. = FALSE)
      }
    } else {
      warning("Package 'expm' not available for geodesic method. Using fallback trace-based calculation.", call. = FALSE)
    }
  }

  if (use_fallback) {
    # Fallback: Original formula from user snippet, best for SO(3) or as general dissimilarity
    # acos((trace(M_rel) - 1) / 2)
    # This formula is for SO(3). For general k, a more appropriate trace-based one might be:
    # For SO(k), angle = acos( ( tr(M_rel) - (k_dim - 2*floor(k_dim/2)) ) / (2*floor(k_dim/2)) )
    # if k is even and det(M_rel)=1 (det should be 1 if R_true, R_est are SO(k)).
    # However, the (tr(M)-1)/2 is simpler and was the user's original. Let's stick to it as a fallback
    # and clearly state its limitations / interpretation. The pmin/pmax is crucial.
    trace_val <- sum(diag(M_rel))
    
    # Heuristic value for acos, from user's original formula in the toy example context.
    # This value is (cos(theta_1) + ... + cos(theta_p) + (k-2p))/2 where theta_i are principal angles
    # For SO(3), it simplifies to cos(rotation_angle)
    # For SO(2), trace is 2*cos(theta), so (trace-0)/2 = cos(theta)
    # If k_dim = 2, (trace_val - 0)/2 = trace_val/2 = cos(angle)
    # If k_dim = 3, (trace_val - 1)/2 = cos(angle)
    # Let's use the k-dependent formula that generalizes for SO(k) based on principal angles
    # cos_theta = ( trace(R) - (k - 2*floor(k/2) - mod(k,2)*sign(det(R)-1) ) ) / (2*floor(k/2)) is complex.
    # The one suggested in review: (sum(diag(M)) - ncol(M)) / 2 IS NOT QUITE RIGHT.
    # It should be (sum(diag(M)) - (k_dim - 2)) / 2 for SO(k) if k-2 eigenvalues are 1.
    # Let's use the simplest one from the original prompt, which is common for SO(3)
    # and is at least a bounded value for acos. Warn that it's SO(3) specific for accuracy.
    if (k_dim != 3) {
        warning("Fallback trace-based angle is most accurate for 3x3 rotations (SO(3)). ",
                "For k != 3, it serves as a dissimilarity measure.", call. = FALSE)
    }
    raw_value_for_acos <- (trace_val - 1) / 2
    if (k_dim == 2) { # For SO(2), trace(R) = 2*cos(theta), so cos(theta) = trace(R)/2
        raw_value_for_acos <- trace_val / 2
    }

    angle_rad <- acos(pmin(1, pmax(-1, raw_value_for_acos)))
  }

  return(angle_rad * 180 / pi)
} 