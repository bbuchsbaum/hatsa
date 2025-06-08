#' Solve Generalized Eigenvalue Problem for Laplacians using PRIMME
#'
#' Solves the generalized eigenvalue problem `A v = λ B v` where A and B are
#' typically sparse graph Laplacians (e.g., `L_task` and `L_conn`). It finds
#' eigenvectors corresponding to the smallest magnitude eigenvalues (`λ`).
#'
#' @param A The sparse, symmetric matrix on the left side (e.g., `L_task`, `dgCMatrix`).
#' @param B The sparse, symmetric, positive semi-definite matrix on the right
#'   side (e.g., `L_conn`, `dgCMatrix`). Will be regularized.
#' @param k_request Integer, the number of eigenvalues/vectors to compute (`NEig`).
#' @param lambda_max_thresh Numeric, the maximum absolute eigenvalue (`|λ|`)
#'   to retain. Eigenpairs with `abs(values) >= lambda_max_thresh` are discarded.
#' @param epsilon_reg_B Numeric, small value to add to the diagonal of B for
#'   regularization (`B_reg = B + epsilon_reg_B * I`). Helps ensure B is
#'   positive definite for the solver. Default 1e-6.
#' @param tol Numeric, tolerance for eigenvalue decomposition convergence.
#'   Default 1e-8 (PRIMME's default is 1e-6, using slightly tighter).
#' @param primme_which Character string passed to `PRIMME::eigs_sym` to
#'   control which eigenvalues are computed. Defaults to
#'   "primme_closest_abs" which targets the eigenvalues closest to zero.
#' @param ... Additional arguments passed to `PRIMME::eigs_sym`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\code{vectors}: A dense matrix (`V_p x k_actual`) of filtered eigenvectors.}
#'     \item{\code{values}: A numeric vector of the corresponding filtered eigenvalues.}
#'     \item{\code{n_converged}: The number of eigenpairs PRIMME reports converged.}
#'     \item{\code{n_filtered}: The number of eigenpairs remaining after filtering.}
#'     \item{\code{primme_stats}: The stats list returned by PRIMME.}
#'   }
#'   Throws an error if computation fails. Note: Stability
#'   filtering based on split-half reliability (`r_split`) needs to be applied
#'   separately. Eigenvectors are B-orthogonal.
#'
#' @importFrom Matrix Diagonal t
#' @importFrom PRIMME eigs_sym
#' @importFrom methods is
#' @keywords internal
solve_gev_laplacian_primme <- function(A, B, k_request,
                                         lambda_max_thresh = 0.8,
                                         epsilon_reg_B = 1e-6,
                                         tol = 1e-8,
                                         primme_which = "primme_closest_abs",
                                         ...) {

  if (!requireNamespace("PRIMME", quietly = TRUE)) {
      stop("The 'PRIMME' package is required for GEV solving. Please install it.")
  }
    
  if (!is(A, "sparseMatrix") || !is(B, "sparseMatrix")) {
    stop("Inputs A and B must be sparse matrices.")
  }
  V_p <- nrow(A)
  if (V_p == 0) return(list(vectors=matrix(0,0,0), values=numeric(0), n_converged=0, n_filtered=0, primme_stats=list()))
  if (nrow(B) != V_p || ncol(A) != V_p || ncol(B) != V_p) {
    stop("Input matrices A and B must be square and of the same dimension.")
  }

  # Ensure k is valid
  k_solve <- min(k_request, V_p)
  if (k_solve <= 0) {
    warning("k_request must be positive. Cannot solve GEV.")
    return(list(vectors=matrix(0,V_p,0), values=numeric(0), n_converged=0, n_filtered=0, primme_stats=list()))
  }

  # Regularize B: B_reg = B + epsilon * I
  if (epsilon_reg_B > 0) {
    # Safer alternative: explicitly create B_reg without modifying B in the calling scope
    B_reg <- B + epsilon_reg_B * Matrix::Diagonal(V_p)
  } else {
    B_reg <- B # No regularization if epsilon is zero or negative
  }

  # Solve the generalized eigenvalue problem A*v = lambda*B_reg*v using PRIMME
  message_stage(sprintf("Solving generalized eigenvalue problem with PRIMME for %d components...", k_solve), interactive_only = TRUE)
  gev_result <- tryCatch({
    PRIMME::eigs_sym(
      A = A,
      NEig = k_solve,
      B = B_reg,
      which = primme_which,
      tol = tol,
      ...)
  }, error = function(e) {
    stop(paste("PRIMME generalized eigenvalue decomposition failed:", e$message))
  })

  n_converged <- length(gev_result$values)
  if (n_converged == 0) {
     warning("PRIMME GEV solver did not converge for any eigenpairs.")
     return(list(vectors = matrix(0, V_p, 0), values = numeric(0), n_converged = 0, n_filtered = 0, primme_stats = gev_result$stats))
  }
  message_stage(sprintf("PRIMME GEV solver converged for %d/%d eigenpairs.", sum(!is.na(gev_result$values)), k_solve), interactive_only = TRUE)

  # Filter based on lambda_max_thresh
  # PRIMME might return NAs for non-converged values
  converged_idx <- !is.na(gev_result$values)
  values_raw <- gev_result$values[converged_idx]
  vectors_raw <- gev_result$vectors[, converged_idx, drop = FALSE] # V_p x n_converged_valid

  # Use absolute value for thresholding magnitude
  filter_idx <- which(abs(values_raw) < lambda_max_thresh)

  if (length(filter_idx) == 0) {
      message_stage(sprintf("No converged eigenvalues passed the lambda_max_thresh < %.3f filter.", lambda_max_thresh), interactive_only = TRUE)
      vectors_filtered <- matrix(0, V_p, 0)
      values_filtered <- numeric(0)
  } else {
      vectors_filtered <- vectors_raw[, filter_idx, drop = FALSE]
      values_filtered <- values_raw[filter_idx]

      # Sort filtered results by eigenvalue magnitude (absolute value) for consistency
      sorted_order <- order(abs(values_filtered))
      values_filtered <- values_filtered[sorted_order]
      vectors_filtered <- vectors_filtered[, sorted_order, drop = FALSE]
      
      message_stage(sprintf("%d converged eigenvalues passed the lambda_max_thresh < %.3f filter.", length(filter_idx), lambda_max_thresh), interactive_only = TRUE)
  }
  
  # Note: Further filtering based on stability (e.g., split-half r > 0.6)
  # and splitting into labGEV/gevShared patches based on value ranges
  # should be done by the calling function.

  return(list(
    vectors = vectors_filtered,
    values = values_filtered,
    n_converged = sum(converged_idx), # Number of non-NA eigenvalues returned
    n_filtered = length(filter_idx),
    primme_stats = gev_result$stats  # Pass along PRIMME stats
  ))
}

#' @title Compute GEV Spectrum Diagnostics
#' @description Computes various diagnostic statistics for a vector of Generalized Eigenvalues (GEV).
#'
#' @param Lambda_GEV A numeric vector of eigenvalues obtained from GEV.
#' @param lambda_max_thresh A numeric threshold used to categorize eigenvalues.
#'        Eigenvalues with absolute value less than this are considered 'retained' or 'stable'.
#'
#' @return A list containing the following diagnostic statistics:
#'   \\itemize{
#'     \\item{n_eigenvalues: Total number of eigenvalues.}
#'     \\item{min_eigenvalue: Minimum eigenvalue.}
#'     \\item{max_eigenvalue: Maximum eigenvalue.}
#'     \\item{mean_eigenvalue: Mean of eigenvalues.}
#'     \\item{median_eigenvalue: Median of eigenvalues.}
#'     \\item{sd_eigenvalue: Standard deviation of eigenvalues.}
#'     \\item{n_below_thresh: Number of eigenvalues with absolute value < lambda_max_thresh.}
#'     \\item{prop_below_thresh: Proportion of eigenvalues with absolute value < lambda_max_thresh.}
#'     \\item{n_above_thresh: Number of eigenvalues with absolute value >= lambda_max_thresh.}
#'     \\item{prop_above_thresh: Proportion of eigenvalues with absolute value >= lambda_max_thresh.}
#'   }
#' @export
#' @examples
#'   Lambda_GEV_sample <- c(0.1, 0.5, 0.85, 0.95, 1.2)
#'   compute_gev_spectrum_diagnostics(Lambda_GEV_sample, lambda_max_thresh = 0.8)
compute_gev_spectrum_diagnostics <- function(Lambda_GEV, lambda_max_thresh) {
  if (!is.numeric(Lambda_GEV) || !is.vector(Lambda_GEV)) {
    stop("Lambda_GEV must be a numeric vector.")
  }
  if (!is.numeric(lambda_max_thresh) || length(lambda_max_thresh) != 1) {
    stop("lambda_max_thresh must be a single numeric value.")
  }
  if (length(Lambda_GEV) == 0) {
    return(
      list(
        n_eigenvalues = 0,
        min_eigenvalue = NA_real_,
        max_eigenvalue = NA_real_,
        mean_eigenvalue = NA_real_,
        median_eigenvalue = NA_real_,
        sd_eigenvalue = NA_real_,
        n_below_thresh = 0,
        prop_below_thresh = NA_real_,
        n_above_thresh = 0,
        prop_above_thresh = NA_real_
      )
    )
  }

  n_eigenvalues <- length(Lambda_GEV)
  abs_Lambda_GEV <- abs(Lambda_GEV)

  n_below_thresh <- sum(abs_Lambda_GEV < lambda_max_thresh)
  prop_below_thresh <- n_below_thresh / n_eigenvalues
  n_above_thresh <- n_eigenvalues - n_below_thresh
  prop_above_thresh <- n_above_thresh / n_eigenvalues

  stats <- list(
    n_eigenvalues = n_eigenvalues,
    min_eigenvalue = min(Lambda_GEV, na.rm = TRUE),
    max_eigenvalue = max(Lambda_GEV, na.rm = TRUE),
    mean_eigenvalue = mean(Lambda_GEV, na.rm = TRUE),
    median_eigenvalue = stats::median(Lambda_GEV, na.rm = TRUE),
    sd_eigenvalue = stats::sd(Lambda_GEV, na.rm = TRUE),
    n_below_thresh = n_below_thresh,
    prop_below_thresh = prop_below_thresh,
    n_above_thresh = n_above_thresh,
    prop_above_thresh = prop_above_thresh
  )

  return(stats)
} 