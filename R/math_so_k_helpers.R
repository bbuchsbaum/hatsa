#' Fast Fr\u00e9chet Mean on SO(k)
#'
#' Computes an approximate Fr\u00e9chet mean of a list of rotation matrices
#' using a chordal (SVD) mean followed by an optional single Karcher
#' refinement step. If refinement does not converge, it falls back to the
#' full Karcher mean implementation.
#'
#' @param Rlist A list of k x k rotation matrices.
#' @param refine Logical; if TRUE a single Karcher refinement step is
#'   attempted. Default TRUE.
#' @param tol Numeric tolerance used for convergence checks. Default 1e-8.
#'
#' @return A k x k matrix representing the Fr\u00e9chet mean.
#' @keywords internal
#' @importFrom expm logm expm
frechet_mean_so_fast <- function(Rlist, refine = TRUE, tol = 1e-8) {
  if (!all(sapply(Rlist, function(R) is.matrix(R) && inherits(R, "matrix")))) {
    stop("All elements in Rlist must be standard R matrices.")
  }
  valid_Rlist <- Filter(Negate(is.null), Rlist)
  if (length(valid_Rlist) == 0) {
    stop("Cannot determine dimension for Fr\u00e9chet mean with empty or all-NULL Rlist.")
  }
  k <- nrow(valid_Rlist[[1]])
  if (k == 0) return(matrix(numeric(0), 0, 0))

  N <- length(valid_Rlist)
  M <- matrix(0, nrow = k, ncol = k)
  for (R_i in valid_Rlist) {
    if (is.matrix(R_i) && all(dim(R_i) == c(k, k))) {
      M <- M + R_i
    }
  }
  M <- M / N

  sv <- svd(M)
  S_diag_vals <- rep(1, k)
  if (det(sv$u %*% t(sv$v)) < 0) S_diag_vals[k] <- -1
  Rbar_chordal <- sv$u %*% diag(S_diag_vals, nrow = k, ncol = k) %*% t(sv$v)

  if (!refine) return(Rbar_chordal)

  log_sum_tangent <- matrix(0, nrow = k, ncol = k)
  num_valid_for_refine <- 0
  for (R_i in valid_Rlist) {
    if (is.matrix(R_i) && all(dim(R_i) == c(k, k))) {
      A <- Rbar_chordal %*% t(R_i)
      log_A <- tryCatch(expm::logm(A), error = function(e) {
        warning("logm failed in frechet_mean_so_fast refinement. Using zero matrix.")
        matrix(0, k, k)
      })
      log_sum_tangent <- log_sum_tangent + log_A
      num_valid_for_refine <- num_valid_for_refine + 1
    }
  }

  if (num_valid_for_refine == 0) {
    warning("No valid matrices for refinement step in frechet_mean_so_fast. Returning chordal mean.")
    return(Rbar_chordal)
  }
  mean_update_log <- log_sum_tangent / num_valid_for_refine
  Rbar_refined <- expm::expm(mean_update_log) %*% Rbar_chordal

  if (max(abs(mean_update_log)) < tol) return(Rbar_refined)

  frechet_mean_so_karcher(valid_Rlist, R_init = Rbar_refined, tol = tol)
}

#' Full Karcher Mean on SO(k)
#'
#' Iteratively computes the Karcher mean of a list of rotation matrices.
#'
#' @param Rlist A list of k x k rotation matrices.
#' @param R_init Optional initial estimate. Defaults to chordal mean.
#' @param tol Convergence tolerance. Default 1e-10.
#' @param maxit Maximum number of iterations. Default 100.
#'
#' @return A k x k matrix giving the Karcher mean.
#' @keywords internal
#' @importFrom expm logm expm
frechet_mean_so_karcher <- function(Rlist, R_init = NULL, tol = 1e-10, maxit = 100L) {
  valid_Rlist <- Filter(Negate(is.null), Rlist)
  if (length(valid_Rlist) == 0) {
    k_dim_fallback <- if (!is.null(R_init) && is.matrix(R_init)) nrow(R_init) else 0
    if (k_dim_fallback > 0) return(diag(k_dim_fallback))
    stop("Cannot determine dimension for Karcher mean.")
  }
  k <- nrow(valid_Rlist[[1]])
  if (k == 0) return(matrix(numeric(0), 0, 0))

  Rbar_current <- if (is.null(R_init)) {
    frechet_mean_so_fast(valid_Rlist, refine = FALSE)
  } else {
    R_init
  }
  N <- length(valid_Rlist)

  for (it in seq_len(maxit)) {
    log_sum_tangent <- matrix(0, k, k)
    num_valid_R <- 0
    for (R_i in valid_Rlist) {
      if (is.matrix(R_i) && all(dim(R_i) == c(k, k))) {
        err_rot <- t(Rbar_current) %*% R_i
        log_A_i <- tryCatch(expm::logm(err_rot), error = function(e) {
          warning(sprintf("logm failed in Karcher iteration %d. Using zero matrix.", it))
          matrix(0, k, k)
        })
        log_sum_tangent <- log_sum_tangent + log_A_i
        num_valid_R <- num_valid_R + 1
      }
    }
    if (num_valid_R == 0) {
      warning("No valid matrices left in Rlist during Karcher iteration. Returning current mean.")
      return(Rbar_current)
    }
    mean_tangent_update <- log_sum_tangent / num_valid_R
    Rbar_new <- Rbar_current %*% expm::expm(mean_tangent_update)
    if (max(abs(mean_tangent_update)) < tol) {
      Rbar_current <- Rbar_new
      break
    }
    Rbar_current <- Rbar_new
    if (it == maxit) warning("Karcher mean did not converge in ", maxit, " iterations.")
  }
  Rbar_current
}
