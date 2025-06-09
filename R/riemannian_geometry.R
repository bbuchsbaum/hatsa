#' @keywords internal
.regularize_spd <- function(S, epsilon_spd = 1e-6) {
  if (!is.matrix(S) || !is.numeric(S) || !isSymmetric.matrix(S, tol = 1e-9)) { # Check symmetry
    stop("Input S must be a numeric, symmetric matrix.")
  }
  if (nrow(S) == 0) return(S)

  eig <- eigen(S, symmetric = TRUE)
  vals <- eig$values
  vals_reg <- ifelse(vals < epsilon_spd, epsilon_spd, vals)
  S_reg <- eig$vectors %*% diag(vals_reg) %*% t(eig$vectors)
  S_reg_symm <- (S_reg + t(S_reg)) / 2
  return(S_reg_symm)
}

#' Compute Riemannian Distance Between SPD Matrices (Log-Euclidean Metric)
#'
#' Calculates the Riemannian distance between two symmetric positive-definite (SPD)
#' matrices S1 and S2, using the Log-Euclidean metric as described by
#' Mitteroecker & Bookstein (2008) and Arsigny et al. (2006).
#' `d(S1, S2) = ||logm(S1) - logm(S2)||_F` which, for the M&B formulation, is
#' effectively `sqrt(sum(log(lambda_j(S2^-1 S1))^2))`.
#' This function is a core component of RGEOM-001.
#'
#' @param S1 A numeric, symmetric matrix (p x p).
#' @param S2 A numeric, symmetric matrix (p x p).
#' @param regularize_epsilon A small positive value used for regularization to
#'   ensure matrices are SPD before inversion and eigenvalue computation. Default 1e-6.
#' @param eigenvalue_floor A small positive value to threshold relative eigenvalues
#'   before taking the logarithm, preventing issues with `log(0)`. Default 1e-9.
#'
#' @return The Riemannian distance (a scalar). Returns `NA` if computation fails.
#' @export
#' @examples
#' S1 <- matrix(c(2.3, -0.3, -0.3, 3.6), 2, 2)
#' S2 <- matrix(c(3.7, 1.9, 1.9, 2.8), 2, 2)
#' # riemannian_distance_spd(S1, S2) # Expected: ~1.24156 (Need to import MASS for ginv)
#'
#' # Identical matrices
#' # riemannian_distance_spd(S1, S1) # Expected: 0
riemannian_distance_spd <- function(S1, S2, regularize_epsilon = 1e-6, eigenvalue_floor = 1e-9) {
  if (!is.matrix(S1) || !is.matrix(S2) ||
      !is.numeric(S1) || !is.numeric(S2) ||
      !isSymmetric.matrix(S1, tol = 1e-9) || !isSymmetric.matrix(S2, tol = 1e-9) ||
      nrow(S1) != ncol(S1) || nrow(S2) != ncol(S2) ||
      nrow(S1) != nrow(S2)) {
    stop("S1 and S2 must be numeric, square, symmetric matrices of the same dimension.")
  }
  p <- nrow(S1)
  if (p == 0) return(0)

  S1_reg <- .regularize_spd(S1, regularize_epsilon)
  S2_reg <- .regularize_spd(S2, regularize_epsilon)

  # Compute the Log-Euclidean distance: ||logm(S1) - logm(S2)||_F
  logS1 <- matrix_logm_spd(S1_reg, regularize_epsilon)
  logS2 <- matrix_logm_spd(S2_reg, regularize_epsilon)
  
  # Frobenius norm of the difference
  diff_log <- logS1 - logS2
  distance <- sqrt(sum(diff_log^2))
  
  return(distance)
}

#' Matrix Logarithm of an SPD Matrix
#'
#' Computes the principal matrix logarithm of a symmetric positive-definite (SPD)
#' matrix. The input matrix is first regularized to ensure positive definiteness.
#'
#' @param S A numeric, symmetric matrix.
#' @param regularize_epsilon Epsilon for regularization. Default 1e-6.
#' @return The matrix logarithm of S (a symmetric matrix).
#' @export
#' @importFrom expm logm
#' @examples
#' S1 <- matrix(c(2.3, -0.3, -0.3, 3.6), 2, 2)
#' logS1 <- matrix_logm_spd(S1)
#' print(logS1)
matrix_logm_spd <- function(S, regularize_epsilon = 1e-6) {
  if (!is.matrix(S) || !is.numeric(S) || !isSymmetric.matrix(S, tol = 1e-9)) {
    stop("Input S must be a numeric, symmetric matrix.")
  }
  if (nrow(S) == 0) return(S)  # Empty matrix case
  S_reg <- .regularize_spd(S, regularize_epsilon)
  eig <- eigen(S_reg, symmetric = TRUE)
  vals <- pmax(eig$values, regularize_epsilon)
  log_S <- eig$vectors %*% diag(log(vals)) %*% t(eig$vectors)
  log_S_symm <- (log_S + t(log_S)) / 2
  return(log_S_symm)
}

#' Matrix Exponential of a Symmetric Matrix
#'
#' Computes the matrix exponential of a symmetric matrix. The result is an SPD matrix.
#'
#' @param S_symm A numeric, symmetric matrix (e.g., a matrix in the tangent space).
#' @return The matrix exponential of S_symm (an SPD matrix).
#' @export
#' @importFrom expm expm
#' @examples
#' S1 <- matrix(c(2.3, -0.3, -0.3, 3.6), 2, 2)
#' logS1 <- matrix_logm_spd(S1)
#' S1_reconstructed <- matrix_expm_spd(logS1)
#' # all.equal(S1, S1_reconstructed) # Should be TRUE, up to numerical precision
matrix_expm_spd <- function(S_symm) {
  if (!is.matrix(S_symm) || !is.numeric(S_symm) || !isSymmetric.matrix(S_symm, tol = 1e-9)) {
    stop("Input S_symm must be a numeric, symmetric matrix.")
  }
  if (nrow(S_symm) == 0) return(S_symm)  # Empty matrix case
  eig <- eigen(S_symm, symmetric = TRUE)
  exp_S <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
  exp_S_symm <- (exp_S + t(exp_S)) / 2
  return(.regularize_spd(exp_S_symm))
}

#' Matrix Square Root of an SPD Matrix
#'
#' Computes the principal matrix square root of a symmetric positive-definite (SPD)
#' matrix. The input matrix is first regularized to ensure positive definiteness.
#' The matrix square root S_sqrt is such that S_sqrt %*% S_sqrt = S.
#'
#' @param S A numeric, symmetric matrix.
#' @param regularize_epsilon Epsilon for regularization. Default 1e-6.
#' @return The matrix square root of S (an SPD matrix).
#' @export
#' @importFrom expm sqrtm
#' @examples
#' S1 <- matrix(c(2.3, -0.3, -0.3, 3.6), 2, 2)
#' S1_sqrt <- matrix_sqrt_spd(S1)
#' # S1_reconstructed_from_sqrt <- S1_sqrt %*% S1_sqrt
#' # all.equal(S1, S1_reconstructed_from_sqrt) # Should be TRUE
matrix_sqrt_spd <- function(S, regularize_epsilon = 1e-6) {
  if (!is.matrix(S) || !is.numeric(S) || !isSymmetric.matrix(S, tol = 1e-9)) {
    stop("Input S must be a numeric, symmetric matrix.")
  }
  if (nrow(S) == 0) return(S)  # Empty matrix case
  S_reg <- .regularize_spd(S, regularize_epsilon)
  eig <- eigen(S_reg, symmetric = TRUE)
  vals <- pmax(eig$values, regularize_epsilon)
  sqrt_S <- eig$vectors %*% diag(sqrt(vals)) %*% t(eig$vectors)
  sqrt_S_symm <- (Re(sqrt_S) + t(Re(sqrt_S))) / 2
  return(.regularize_spd(sqrt_S_symm, regularize_epsilon))
}

#' Affine-Invariant Riemannian Metric (AIRM) Distance
#'
#' Computes the Affine-Invariant Riemannian Metric (AIRM) distance between two
#' symmetric positive-definite (SPD) matrices S1 and S2.
#' The distance is defined as: `||logm(S1^(-1/2) %*% S2 %*% S1^(-1/2))||_F`,
#' where `||.||_F` is the Frobenius norm.
#'
#' @param S1 A numeric, symmetric positive-definite matrix.
#' @param S2 A numeric, symmetric positive-definite matrix.
#' @param regularize_epsilon Epsilon for regularization. Default 1e-6.
#' @return The AIRM distance (a non-negative scalar).
#' @export
#' @examples
#' S1 <- matrix(c(2.3, -0.3, -0.3, 3.6), 2, 2)
#' S2 <- matrix(c(3.7, 1.9, 1.9, 2.8), 2, 2)
#' # dist_airm <- airm_distance(S1, S2)
#' # print(dist_airm)
#' # dist_logeuclidean <- riemannian_distance_spd(S1, S2)
#' # print(dist_logeuclidean) # Note: AIRM and LogEuclidean are different metrics.
airm_distance <- function(S1, S2, regularize_epsilon = 1e-6) {
  if (!is.matrix(S1) || !is.numeric(S1) || !isSymmetric.matrix(S1, tol = 1e-9) ||
      !is.matrix(S2) || !is.numeric(S2) || !isSymmetric.matrix(S2, tol = 1e-9) ||
      nrow(S1) != ncol(S1) || nrow(S2) != ncol(S2) || nrow(S1) != nrow(S2)) {
    stop("S1 and S2 must be numeric, square, symmetric matrices of the same dimension.")
  }

  # Regularize inputs to ensure they are SPD
  S1_reg <- .regularize_spd(S1, regularize_epsilon)
  S2_reg <- .regularize_spd(S2, regularize_epsilon)

  # Compute S1_reg^(-1/2)
  # S1_sqrt_inv <- inv(S1_reg^(1/2))
  S1_sqrt <- matrix_sqrt_spd(S1_reg, regularize_epsilon)
  
  S1_sqrt_inv <- tryCatch(solve(S1_sqrt),
                          error = function(e) {
                            warning("solve(S1_sqrt) failed, trying pseudo-inverse. AIRM distance may be less stable. Error: ", e$message)
                            if (!requireNamespace("MASS", quietly = TRUE)) {
                                stop("MASS package needed for ginv fallback but not available. Original error for solve(): ", e$message)
                            }
                            MASS::ginv(S1_sqrt)
                          })

  # Compute the argument of logm: S1_sqrt_inv %*% S2_reg %*% S1_sqrt_inv
  # This intermediate matrix should also be SPD.
  inner_matrix <- S1_sqrt_inv %*% S2_reg %*% S1_sqrt_inv
  inner_matrix_reg <- .regularize_spd(inner_matrix, regularize_epsilon) # Regularize before logm
  
  # Compute logm of the inner matrix
  log_inner_matrix <- matrix_logm_spd(inner_matrix_reg, regularize_epsilon)
  
  # Frobenius norm of the result: sqrt(sum(diag(t(X) %*% X))) or sqrt(sum(X^2)) for symmetric X
  # Since log_inner_matrix is symmetric:
  distance <- sqrt(sum(log_inner_matrix^2))
  
  return(distance)
}

#' Log Map for SPD Matrices (Log-Euclidean Metric)
#'
#' Projects SPD matrix S2 to the tangent space of SPD matrix S1
#' using the Log-Euclidean metric.
#' The tangent vector is `logm(S2) - logm(S1)`.
#'
#' @param S1 Reference SPD matrix (point on the manifold where the tangent space is anchored).
#' @param S2 SPD matrix to project to the tangent space at S1.
#' @param regularize_epsilon Epsilon for regularization. Default 1e-6.
#' @return A symmetric matrix representing S2 in the tangent space at S1.
#' @export
logmap_spd_logeuclidean <- function(S1, S2, regularize_epsilon = 1e-6) {
  if (!is.matrix(S1) || !isSymmetric.matrix(S1, tol = 1e-9) ||
      !is.matrix(S2) || !isSymmetric.matrix(S2, tol = 1e-9)) {
    stop("S1 and S2 must be symmetric matrices.")
  }
  logS1 <- matrix_logm_spd(S1, regularize_epsilon)
  logS2 <- matrix_logm_spd(S2, regularize_epsilon)
  return(logS2 - logS1)
}

#' Exp Map for SPD Matrices (Log-Euclidean Metric)
#'
#' Maps a tangent vector V (symmetric matrix) from the tangent space at SPD matrix S1
#' back to the SPD manifold using the Log-Euclidean metric.
#' Result is `expm(logm(S1) + V)`.
#'
#' @param S1 SPD matrix (point on manifold where tangent space is anchored).
#' @param V_at_S1 Symmetric matrix (tangent vector at S1).
#' @param regularize_epsilon Epsilon for regularization of S1. Default 1e-6.
#' @return An SPD matrix on the manifold.
#' @export
expmap_spd_logeuclidean <- function(S1, V_at_S1, regularize_epsilon = 1e-6) {
  if (!is.matrix(S1) || !isSymmetric.matrix(S1, tol = 1e-9) ||
      !is.matrix(V_at_S1) || !isSymmetric.matrix(V_at_S1, tol = 1e-9)) {
    stop("S1 and V_at_S1 must be symmetric matrices.")
  }
  logS1 <- matrix_logm_spd(S1, regularize_epsilon)
  # The argument to expm is logm(S1) + V_at_S1, which must be symmetric.
  # logS1 is symmetric, V_at_S1 is assumed symmetric by input validation.
  return(matrix_expm_spd(logS1 + V_at_S1))
}

#' Log Map for SPD Matrices (AIRM Metric)
#'
#' Projects SPD matrix S2 to the tangent space of SPD matrix S1
#' using the Affine-Invariant Riemannian Metric (AIRM).
#' Tangent vector is `S1^(1/2) %*% logm(S1^(-1/2) %*% S2 %*% S1^(-1/2)) %*% S1^(1/2)`.
#'
#' @param S1 Reference SPD matrix.
#' @param S2 SPD matrix to project.
#' @param regularize_epsilon Epsilon for regularization. Default 1e-6.
#' @return A symmetric matrix (tangent vector at S1).
#' @export
logmap_spd_airm <- function(S1, S2, regularize_epsilon = 1e-6) {
   if (!is.matrix(S1) || !isSymmetric.matrix(S1, tol = 1e-9) ||
      !is.matrix(S2) || !isSymmetric.matrix(S2, tol = 1e-9)) {
    stop("S1 and S2 must be symmetric matrices.")
  }
  S1_reg <- .regularize_spd(S1, regularize_epsilon)
  S2_reg <- .regularize_spd(S2, regularize_epsilon)
  
  S1_sqrt <- matrix_sqrt_spd(S1_reg, regularize_epsilon)
  S1_sqrt_inv <- tryCatch(solve(S1_sqrt), error = function(e) {
      warning("solve(S1_sqrt) failed in logmap_spd_airm, trying pseudo-inverse. Error: ", e$message)
      if(!requireNamespace("MASS", quietly=TRUE)) stop("MASS needed for ginv.")
      MASS::ginv(S1_sqrt)
  })
  
  inner_logm_arg <- S1_sqrt_inv %*% S2_reg %*% S1_sqrt_inv
  log_inner <- matrix_logm_spd(inner_logm_arg, regularize_epsilon)
  
  tangent_vector <- S1_sqrt %*% log_inner %*% S1_sqrt
  # Ensure symmetry for the tangent vector
  return((tangent_vector + t(tangent_vector)) / 2)
}

#' Exp Map for SPD Matrices (AIRM Metric)
#'
#' Maps a tangent vector V (symmetric matrix) from the tangent space at SPD matrix S1
#' back to the SPD manifold using the AIRM.
#' Result is `S1^(1/2) %*% expm(S1^(-1/2) %*% V %*% S1^(-1/2)) %*% S1^(1/2)`.
#'
#' @param S1 SPD matrix (point on manifold).
#' @param V_at_S1 Symmetric matrix (tangent vector at S1).
#' @param regularize_epsilon Epsilon for regularization. Default 1e-6.
#' @return An SPD matrix on the manifold.
#' @export
expmap_spd_airm <- function(S1, V_at_S1, regularize_epsilon = 1e-6) {
  if (!is.matrix(S1) || !isSymmetric.matrix(S1, tol = 1e-9) ||
      !is.matrix(V_at_S1) || !isSymmetric.matrix(V_at_S1, tol = 1e-9)) {
    stop("S1 and V_at_S1 must be symmetric matrices.")
  }
  S1_reg <- .regularize_spd(S1, regularize_epsilon)
  V_at_S1_symm <- (V_at_S1 + t(V_at_S1)) / 2 # Ensure V is symmetric
  
  S1_sqrt <- matrix_sqrt_spd(S1_reg, regularize_epsilon)
  S1_sqrt_inv <- tryCatch(solve(S1_sqrt), error = function(e) {
      warning("solve(S1_sqrt) failed in expmap_spd_airm, trying pseudo-inverse. Error: ", e$message)
      if(!requireNamespace("MASS", quietly=TRUE)) stop("MASS needed for ginv.")
      MASS::ginv(S1_sqrt)
  })
  
  inner_expm_arg <- S1_sqrt_inv %*% V_at_S1_symm %*% S1_sqrt_inv
  # The argument to expm must be symmetric for matrix_expm_spd to guarantee SPD output.
  inner_expm_arg_symm <- (inner_expm_arg + t(inner_expm_arg)) / 2
  
  exp_inner <- matrix_expm_spd(inner_expm_arg_symm)
  
  mapped_matrix <- S1_sqrt %*% exp_inner %*% S1_sqrt
  # Ensure the final result is SPD
  return(.regularize_spd((mapped_matrix + t(mapped_matrix)) / 2, regularize_epsilon))
}

#' Fréchet Mean of SPD Matrices
#'
#' Computes the Fréchet mean (Karcher mean) of a list of Symmetric 
#' Positive-Definite (SPD) matrices.
#' This function can use either the Affine-Invariant Riemannian Metric (AIRM)
#' via the `shapes` package (if available and `metric="airm"`), or an iterative
#' algorithm based on the Log-Euclidean metric (`metric="logeuclidean"`).
#'
#' @param S_list A list of SPD matrices.
#' @param metric Character string, either `"logeuclidean"` (default) or `"airm"`.
#' @param regularize_epsilon Epsilon for regularization. Default 1e-6.
#' @param max_iter Maximum number of iterations for Log-Euclidean algorithm. Default 50.
#' @param tol Tolerance for convergence for Log-Euclidean algorithm (Frobenius norm
#'   of the mean tangent vector). Default 1e-5.
#' @param step_size Step size for the gradient descent in Log-Euclidean. Default 0.5.
#' @param init_method For Log-Euclidean, method to initialize the mean: 
#'  `"euclidean"` (Euclidean mean of matrices) or `"first"` (first matrix in list).
#'  Default `"euclidean"`.
#' @param verbose Logical, if TRUE, prints iteration info for Log-Euclidean. Default FALSE.
#' @return The Fréchet mean (an SPD matrix).
#' @export
#' @importFrom stats median
frechet_mean_spd <- function(S_list,
                             metric = c("logeuclidean", "airm"),
                             regularize_epsilon = 1e-6,
                             max_iter = 50, 
                             tol = 1e-5,
                             step_size = 0.5,
                             init_method = c("euclidean", "first"),
                             verbose = FALSE) {
  metric <- match.arg(metric)
  init_method <- match.arg(init_method)

  if (!is.list(S_list) || length(S_list) == 0) {
    stop("S_list must be a non-empty list of matrices.")
  }
  if (!all(sapply(S_list, function(m) is.matrix(m) && isSymmetric.matrix(m, tol = 1e-9)))) {
    stop("All elements in S_list must be symmetric matrices.")
  }
  # Further ensure all matrices are of the same dimension
  p <- nrow(S_list[[1]])
  if (!all(sapply(S_list, function(m) nrow(m) == p && ncol(m) == p))) {
    stop("All matrices in S_list must have the same dimensions.")
  }
  if (p == 0) return(matrix(0,0,0))

  # Regularize all input matrices first
  S_list_reg <- lapply(S_list, function(S) .regularize_spd(S, regularize_epsilon))
  n_matrices <- length(S_list_reg)

  if (metric == "airm") {
    if (!requireNamespace("shapes", quietly = TRUE)) {
      stop("Package 'shapes' needed for AIRM Fréchet mean. Please install it or use metric='logeuclidean'.")
    }
    if (verbose) message("Using AIRM Fréchet mean via shapes::mediancov().")
    # shapes::mediancov expects an array (p, p, N)
    spd_array <- array(unlist(S_list_reg), dim = c(p, p, n_matrices))
    # shapes::mediancov with type="riemann" calculates the Karcher mean using AIRM
    mean_matrix <- shapes::mediancov(spd_array, type = "riemann")$median
    return(.regularize_spd(mean_matrix, regularize_epsilon)) # Ensure final strict SPD
  
  } else if (metric == "logeuclidean") {
    if (verbose) message("Using iterative Log-Euclidean Fréchet mean.")
    
    # Initialization
    current_mean <- if(init_method == "first") {
      S_list_reg[[1]]
    } else { # euclidean mean
      Reduce("+", S_list_reg) / n_matrices
    }
    current_mean <- .regularize_spd(current_mean, regularize_epsilon)

    for (iter in 1:max_iter) {
      tangent_vectors <- lapply(S_list_reg, function(S_i) {
        logmap_spd_logeuclidean(current_mean, S_i, regularize_epsilon)
      })
      
      mean_tangent_vector <- Reduce("+", tangent_vectors) / n_matrices
      norm_mean_tangent <- sqrt(sum(mean_tangent_vector^2)) # Frobenius norm

      if (verbose) {
        message(sprintf("Iter %d: Norm of mean tangent vector = %.2e", iter, norm_mean_tangent))
      }

      if (norm_mean_tangent < tol) {
        if (verbose) message("Converged.")
        break
      }
      
      # Update mean by moving along the mean tangent vector direction
      current_mean <- expmap_spd_logeuclidean(current_mean, mean_tangent_vector * step_size, regularize_epsilon)
      current_mean <- .regularize_spd(current_mean, regularize_epsilon) # Ensure it stays SPD
      
      if (iter == max_iter && verbose) {
        message("Max iterations reached without convergence.")
      }
    }
    return(current_mean)
  }
}

#' Grassmann Distance between Subspaces
#'
#' Computes the Grassmann distance between two k-dimensional subspaces in R^p,
#' represented by orthonormal basis matrices U and V (both p x k).
#' The distance is `sqrt(sum(theta_i^2))`, where `theta_i` are the principal angles
#' between the subspaces.
#' Principal angles `theta_i = acos(s_i)`, where `s_i` are the singular values of `U^T V`.
#'
#' @param U A numeric matrix (p x k) with orthonormal columns.
#' @param V A numeric matrix (p x k) with orthonormal columns.
#' @param tol Tolerance for checking orthonormality. Default 1e-6.
#' @param sv_floor Smallest value for singular values before acos to avoid `acos(>1)`
#'   due to numerical precision. Values are clamped to `[-1+sv_floor, 1-sv_floor]`.
#'   Default 1e-7.
#' @return The Grassmann distance (a non-negative scalar).
#' @export
#' @examples
#' # Example from P. Edelman, T. A. Arias, and A. S. Edelman. "Geometry of algorithms
#' # with orthogonality constraints." SIAM Journal on Matrix Analysis and Applications, 1998.
#' p <- 3; k <- 2
#' U <- qr.Q(qr(matrix(rnorm(p*k), p, k)))
#' V <- qr.Q(qr(matrix(rnorm(p*k), p, k)))
#' # grassmann_distance(U, V)
#'
#' # Identical subspaces
#' # grassmann_distance(U, U) # Should be 0
#'
#' # Orthogonal subspaces (if k <= p/2)
#' if (k <= p/2) {
#'   # V_ortho <- qr.Q(qr(matrix(rnorm(p*k), p, k)), complete=TRUE)[, (k+1):(2*k), drop=FALSE]
#'   # if (ncol(V_ortho) == k) grassmann_distance(U, V_ortho) # Should be sqrt(k * (pi/2)^2)
#' }
grassmann_distance <- function(U, V, tol = 1e-6, sv_floor = 1e-7) {
  if (!is.matrix(U) || !is.numeric(U) || !is.matrix(V) || !is.numeric(V)) {
    stop("U and V must be numeric matrices.")
  }
  if (nrow(U) != nrow(V) || ncol(U) != ncol(V)) {
    stop("U and V must have the same dimensions.")
  }
  p <- nrow(U)
  k <- ncol(U)

  if (k == 0) return(0) # Distance between 0-dim subspaces is 0
  if (k > p) stop("k (number of columns) cannot be greater than p (number of rows).")

  # Check orthonormality (optional, but good for ensuring valid input)
  # U^T U should be I_k
  utu <- crossprod(U)
  vtv <- crossprod(V)
  if (!all.equal(utu, diag(k), tolerance = tol) || 
      !all.equal(vtv, diag(k), tolerance = tol)) {
    warning("Columns of U and/or V may not be perfectly orthonormal. Results might be approximate.")
    # Optional: Orthonormalize them? e.g., U <- qr.Q(qr(U))
  }

  # Singular values of U^T V are cos(theta_i)
  # Ensure U^T V is k x k
  UtV <- crossprod(U, V) 
  s <- svd(UtV, nu = 0, nv = 0)$d

  # Clamp singular values to avoid acos domain errors due to precision
  # s_i should be <= 1. We clamp to [-(1-floor), 1-floor] to be safe with acos.
  s_clamped <- pmin(pmax(s, -(1-sv_floor)), (1-sv_floor))
  
  # Principal angles theta_i = acos(s_i)
  # If any s_i are very close to 1 (e.g. > 1 - sv_floor^2), theta_i is ~0.
  # If any s_i are very close to 0, theta_i is ~pi/2.
  principal_angles <- acos(s_clamped)
  
  # Grassmann distance = sqrt(sum(theta_i^2))
  distance <- sqrt(sum(principal_angles^2))
  
  return(distance)
}

#' Compute Fréchet Mean of Rotation Matrices on SO(k)
#'
#' Calculates the Fréchet mean (geometric mean) of a list of k x k rotation
#' matrices. The Fréchet mean is the rotation matrix R_bar that minimizes the
#' sum of squared geodesic distances to all rotation matrices in the list:
#' \code{R_bar = argmin_R sum_i d(R, R_i)^2}.
#' This implementation uses an iterative algorithm involving logarithmic and
#' exponential maps on SO(k).
#'
#' @param R_list A list of k x k matrices, expected to be rotation matrices (in SO(k) or O(k)).
#'   The function attempts to work with matrices close to SO(k).
#' @param k_dim Integer, the dimension of the rotation matrices (e.g., 3 for SO(3)).
#'   If NULL (default), it's inferred from the first valid matrix in `R_list`.
#' @param max_iter Integer, the maximum number of iterations for the algorithm (default: 50).
#' @param tol Numeric, the tolerance for convergence. The algorithm stops when the
#'   Frobenius norm of the mean tangent vector is below this tolerance (default: 1e-7).
#' @param initial_mean Optional. A k x k matrix to use as the initial estimate for the mean.
#'   If NULL, the first valid matrix in `R_list` is used, or an identity matrix if none are valid.
#' @param project_to_SOk Logical. If TRUE (default), after each update, the new mean estimate
#'  is projected to the closest matrix in SO(k) using SVD. This helps maintain numerical stability
#'  and ensures the result is indeed in SO(k).
#'
#' @return A k x k matrix representing the Fréchet mean of the input rotation matrices.
#'   Returns an identity matrix of appropriate dimension if `R_list` is empty, contains no valid
#'   matrices, or if `k_dim` cannot be determined.
#'
#' @export
#' @importFrom expm logm expm
#' @examples
#' # Example for SO(3)
#' if (requireNamespace("expm", quietly = TRUE)) {
#'   R1 <- matrix(c(1,0,0, 0,cos(0.1),-sin(0.1), 0,sin(0.1),cos(0.1)), 3, 3)
#'   R2 <- matrix(c(cos(0.2),-sin(0.2),0, sin(0.2),cos(0.2),0, 0,0,1), 3, 3)
#'   R_list_so3 <- list(R1, R2)
#'   # frechet_mean_so_k(R_list_so3) # k_dim inferred
#'   # frechet_mean_so_k(R_list_so3, k_dim = 3)
#' }
#'
#' # Example for SO(2)
#' if (requireNamespace("expm", quietly = TRUE)) {
#'   theta1 <- pi/4
#'   R1_so2 <- matrix(c(cos(theta1), -sin(theta1), sin(theta1), cos(theta1)), 2, 2)
#'   theta2 <- pi/3
#'   R2_so2 <- matrix(c(cos(theta2), -sin(theta2), sin(theta2), cos(theta2)), 2, 2)
#'   # frechet_mean_so_k(list(R1_so2, R2_so2))
#' }
.so_logm_closed_form <- function(R) {
  k <- nrow(R)
  if (k == 2) {
    # SO(2): logm(R) = theta * [0 -1; 1 0], theta = atan2(R[2,1], R[1,1])
    theta <- atan2(R[2,1], R[1,1])
    return(matrix(c(0, -theta, theta, 0), 2, 2))
  } else if (k == 3) {
    # SO(3): logm(R) = (theta/(2*sin(theta))) * (R - t(R)), theta = acos((tr(R)-1)/2)
    trR <- sum(diag(R))
    theta <- acos(pmin(pmax((trR - 1) / 2, -1), 1))
    if (abs(theta) < 1e-10) return(matrix(0, 3, 3))
    logR <- (theta / (2 * sin(theta))) * (R - t(R))
    return(logR)
  } else {
    stop("Closed-form logm only implemented for SO(2) and SO(3)")
  }
}

frechet_mean_so_k <- function(R_list, k_dim = NULL, max_iter = 50, tol = 1e-7, initial_mean = NULL, project_to_SOk = TRUE) {
  if (!requireNamespace("expm", quietly = TRUE)) {
    stop("Package 'expm' is required for frechet_mean_so_k.")
  }

  # Filter out NULLs and non-matrices
  valid_Rs <- Filter(function(R) is.matrix(R) && !is.null(R), R_list)

  if (length(valid_Rs) == 0) {
    if (!is.null(k_dim) && k_dim > 0) return(diag(k_dim))
    warning("R_list is empty or contains no valid matrices, and k_dim not specified. Returning NULL.")
    return(NULL)
  }

  # Determine k_dim if not provided
  if (is.null(k_dim)) {
    k_dim <- nrow(valid_Rs[[1]])
    if (ncol(valid_Rs[[1]]) != k_dim) {
      warning("First valid matrix in R_list is not square. Cannot infer k_dim. Returning NULL.")
      return(NULL)
    }
  }
  if (k_dim == 0) {
      warning("k_dim is 0. Returning NULL or an empty matrix depending on context for SO(0).")
      return(matrix(0,0,0)) # SO(0) is technically the group with one element, an empty map.
  }

  # Further filter for correct dimensions
  valid_Rs <- Filter(function(R) all(dim(R) == c(k_dim, k_dim)), valid_Rs)
  num_matrices <- length(valid_Rs)

  if (num_matrices == 0) {
    warning(sprintf("No valid %d x %d matrices in R_list. Returning identity.", k_dim, k_dim))
    return(diag(k_dim))
  }

  # Initialize R_mean_current
  if (!is.null(initial_mean) && is.matrix(initial_mean) && all(dim(initial_mean) == c(k_dim, k_dim))) {
    R_mean_current <- initial_mean
  } else {
    R_mean_current <- valid_Rs[[1]]
  }
  # Ensure initial R_mean_current is a valid rotation or close to it for stability
  if (project_to_SOk) {
      svd_init <- svd(R_mean_current)
      R_mean_current <- svd_init$u %*% t(svd_init$v)
      if (det(R_mean_current) < 0) { # Ensure it's SO(k) not just O(k)
          V_prime <- svd_init$v
          V_prime[,k_dim] <- -V_prime[,k_dim]
          R_mean_current <- svd_init$u %*% t(V_prime)
      }
  }

  for (iter in 1:max_iter) {
    tangent_vectors_sum <- matrix(0, nrow = k_dim, ncol = k_dim)
    num_valid_tangents <- 0
    logm_failures <- 0

    for (i in 1:num_matrices) {
      R_i <- valid_Rs[[i]]
      arg_logm <- crossprod(R_mean_current, R_i) # t(R_mean_current) %*% R_i
      tangent_i <- NULL
      if (k_dim == 2 || k_dim == 3) {
        # Use closed-form logm for SO(2)/SO(3)
        try({
          tangent_i <- .so_logm_closed_form(arg_logm)
        }, silent = TRUE)
      }
      if (is.null(tangent_i)) {
        # Fallback to expm::logm
        tangent_i <- tryCatch(expm::logm(arg_logm, method="Higham08.b"), error = function(e) NULL)
      }
      if (!is.null(tangent_i) && !all(tangent_i == 0)) {
        tangent_vectors_sum <- tangent_vectors_sum + tangent_i
        num_valid_tangents <- num_valid_tangents + 1
      } else if (max(abs(arg_logm - diag(k_dim))) < tol) {
        num_valid_tangents <- num_valid_tangents + 1
      } else {
        logm_failures <- logm_failures + 1
      }
    }
    if (logm_failures > 0) {
      warning(sprintf("frechet_mean_so_k: %d logm failures in iteration %d (skipped in mean).", logm_failures, iter))
    }
    # If all logm operations failed or resulted in zero (e.g. all R_i are R_mean_current)
    if (num_valid_tangents == 0) {
        if (all(tangent_vectors_sum == 0)) break # Converged if sum is zero
        warning("No valid tangent vectors could be computed in an iteration of frechet_mean_so_k. Returning current mean.")
        break 
    }
    mean_tangent <- tangent_vectors_sum / num_valid_tangents
    # Update R_mean: R_mean_new = R_mean_current %*% expm(mean_tangent)
    R_mean_new <- R_mean_current %*% expm::expm(mean_tangent)
    if (project_to_SOk) {
      # Normalize R_mean_new to be strictly in SO(k) using SVD
      svd_R_mean_new <- svd(R_mean_new)
      R_mean_new_proj <- svd_R_mean_new$u %*% t(svd_R_mean_new$v)
      if (det(R_mean_new_proj) < 0) { # Ensure it's SO(k) not just O(k)
          V_prime <- svd_R_mean_new$v
          V_prime[,k_dim] <- -V_prime[,k_dim] # Flip last column of V
          R_mean_new_proj <- svd_R_mean_new$u %*% t(V_prime)
      }
      R_mean_new <- R_mean_new_proj
    }
    # Check for convergence: norm of the mean_tangent vector
    if (norm(mean_tangent, type = "F") < tol) {
      R_mean_current <- R_mean_new
      break
    }
    R_mean_current <- R_mean_new
    if (iter == max_iter) {
      warning(sprintf("Fréchet mean on SO(%d) did not converge after %d iterations.", k_dim, max_iter))
    }
  }
  return(R_mean_current)
}

#' Compute Bures-Wasserstein Barycenter of SPD Matrices
#'
#' Calculates the Bures-Wasserstein (BW) barycenter of a list of symmetric
#' positive-definite (SPD) matrices. The BW barycenter is the geodesic mean
#' under the BW metric. This function implements an iterative algorithm.
#'
#' @param S_list A list of SPD matrices (p x p).
#' @param weights Optional. A numeric vector of non-negative weights for each matrix
#'   in `S_list`. If NULL (default), uniform weights (1/N) are used. Must sum to 1
#'   if provided, or will be normalized.
#' @param initial_mean Optional. A p x p SPD matrix to use as the initial estimate
#'   for the barycenter. If NULL, the (weighted) arithmetic mean of `S_list`
#'   is used after regularization.
#' @param max_iter Integer, maximum number of iterations for the fixed-point algorithm.
#'   Default: 50.
#' @param tol Numeric, tolerance for convergence. The algorithm stops when the
#'   Frobenius norm of the difference between successive estimates of the barycenter
#'   is below this tolerance. Default: 1e-7.
#' @param regularize_epsilon Numeric, small positive value for regularizing input
#'   matrices and intermediate results to ensure positive definiteness. Default: 1e-6.
#' @param verbose Logical, if TRUE, prints iteration information. Default: FALSE.
#' @param damping Numeric, damping factor for the update step. Default 0.5.
#'
#' @return A p x p SPD matrix representing the Bures-Wasserstein barycenter.
#'   Returns NULL if computation fails or inputs are invalid.
#' @export
#' @examples
#' # S1 <- matrix(c(2,1,1,2), 2,2)
#' # S2 <- matrix(c(3,0,0,3), 2,2)
#' # S_list_bw <- list(S1, S2)
#' # bw_mean <- bures_wasserstein_barycenter(S_list_bw, verbose = TRUE)
#' # print(bw_mean)
#' 
#' # Weighted example
#' # S3 <- matrix(c(1.5,0.5,0.5,1.5), 2,2)
#' # bw_mean_weighted <- bures_wasserstein_barycenter(list(S1,S2,S3), weights=c(0.5,0.25,0.25))
#' # print(bw_mean_weighted)
bures_wasserstein_barycenter <- function(S_list,
                                         weights = NULL,
                                         initial_mean = NULL,
                                         max_iter = 50,
                                         tol = 1e-7,
                                         regularize_epsilon = 1e-6,
                                         verbose = FALSE,
                                         damping = 0.5) {

  # --- Input Validation and Initialization ---
  if (!is.list(S_list) || length(S_list) == 0) {
    stop("S_list must be a non-empty list of matrices.")
  }
  if (!all(sapply(S_list, function(m) is.matrix(m) && isSymmetric.matrix(m, tol = 1e-9)))) {
    stop("All elements in S_list must be symmetric matrices.")
  }
  p <- nrow(S_list[[1]])
  if (p == 0) return(matrix(0,0,0))
  if (!all(sapply(S_list, function(m) nrow(m) == p && ncol(m) == p))) {
    stop("All matrices in S_list must have the same dimensions (p x p).")
  }

  N <- length(S_list)

  # Validate and normalize weights
  if (is.null(weights)) {
    weights <- rep(1/N, N)
  } else {
    if (!is.numeric(weights) || length(weights) != N || any(weights < 0)) {
      stop("Weights must be a numeric vector of non-negative values, same length as S_list.")
    }
    sum_w <- sum(weights)
    if (abs(sum_w - 1.0) > 1e-6 && sum_w > 1e-9) { # Normalize if not sum to 1 (and sum is not zero)
      if (verbose) message("Normalizing weights to sum to 1.")
      weights <- weights / sum_w
    } else if (sum_w < 1e-9) {
      stop("Sum of weights is too close to zero.")
    }
  }

  # Regularize all input matrices
  S_list_reg <- lapply(S_list, function(S) .regularize_spd(S, regularize_epsilon))

  # Initialize barycenter M_current
  M_current <- initial_mean
  if (is.null(M_current)) {
    # Weighted arithmetic mean as initial guess
    M_sum <- matrix(0, p, p)
    for (i in 1:N) {
      M_sum <- M_sum + weights[i] * S_list_reg[[i]]
    }
    M_current <- M_sum # No division by N as weights sum to 1
  } else {
    if (!is.matrix(M_current) || !isSymmetric.matrix(M_current, tol = 1e-9) || 
        nrow(M_current) != p || ncol(M_current) != p) {
      stop("Provided initial_mean is not a valid p x p symmetric matrix.")
    }
  }
  M_current <- .regularize_spd(M_current, regularize_epsilon)

  # --- Iterative Algorithm (Fixed-Point Iteration) ---
  # Based on Bures-Wasserstein barycenter fixed point equation, e.g.,
  # M_new = (sum_i w_i * (M_current^(1/2) S_i M_current^(1/2))^(1/2) ) M_current^(-1/2) ... and then symmetrize + square
  # Or simpler (from Bhatia, Positive Definite Matrices, Chapter X, or Moakher 2005):
  # M_k+1 = M_k^(1/2) ( sum w_i (M_k^(-1/2) S_i M_k^(-1/2))^(1/2) )^2 M_k^(1/2)
  # This can be unstable. A common iterative scheme is:
  # T_k = sum_i w_i * (M_k^(1/2) S_i M_k^(1/2))^(1/2)
  # M_k+1 = M_k^(-1/2) T_k^2 M_k^(1/2)

  if (verbose) message_stage(sprintf("Starting BW Barycenter iteration (max_iter=%d, tol=%.2e)...", max_iter, tol), interactive_only = TRUE)

  for (iter in 1:max_iter) {
    M_old <- M_current
    M_current_sqrt <- matrix_sqrt_spd(M_current, regularize_epsilon)
    M_current_inv_sqrt <- tryCatch(solve(M_current_sqrt), error = function(e) {
        warning("solve(M_current_sqrt) failed, trying pseudo-inverse. BW barycenter may be less stable.")
        if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS package needed for ginv fallback.")
        MASS::ginv(M_current_sqrt)
    })

    if (any(is.na(M_current_inv_sqrt))) {
        warning("M_current_inv_sqrt contains NAs. Stopping BW iteration.")
        return(M_old) # Or NULL
    }

    T_k_sum <- matrix(0, p, p)
    for (i in 1:N) {
      transformed_S_i <- M_current_inv_sqrt %*% S_list_reg[[i]] %*% M_current_inv_sqrt
      transformed_S_i_reg <- .regularize_spd(transformed_S_i, regularize_epsilon) 
      sqrt_transformed_S_i <- matrix_sqrt_spd(transformed_S_i_reg, regularize_epsilon)
      T_k_sum <- T_k_sum + weights[i] * sqrt_transformed_S_i
    }
    
    # Symmetrize T_k_sum before squaring
    T_k_sym <- (T_k_sum + t(T_k_sum)) / 2
    # Update M_current: M_k+1 = M_k^(1/2) T_k_sym^2 M_k^(1/2)
    M_new <- M_current_sqrt %*% (T_k_sym %*% T_k_sym) %*% M_current_sqrt
    M_new <- (M_new + t(M_new)) / 2 # Ensure symmetry
    M_new <- .regularize_spd(M_new, regularize_epsilon) # Ensure SPD
    # Damping update
    M_current <- (1 - damping) * M_old + damping * M_new
    M_current <- (M_current + t(M_current)) / 2
    M_current <- .regularize_spd(M_current, regularize_epsilon)

    # Check for convergence
    diff_norm <- norm(M_current - M_old, type = "F")
    if (verbose) {
      message(sprintf("  BW Iter %d: Change in M = %.3e", iter, diff_norm))
    }
    if (diff_norm < tol) {
      if (verbose) message_stage("BW Barycenter converged.", interactive_only = TRUE)
      break
    }
    if (iter == max_iter && verbose) {
      message_stage("BW Barycenter reached max iterations without specified convergence.", interactive_only = TRUE)
    }
  }

  return(M_current)
} 
