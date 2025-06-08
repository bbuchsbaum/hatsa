#' @title Project Features onto a Spectral Basis
#' @description Projects a feature matrix onto a given spectral basis. Handles potential
#'   transposition of the feature matrix and uses the appropriate orthogonal projection
#'   formula based on whether the basis is orthonormal. Sparse inputs are processed
#'   with `Matrix` methods so that large matrices need not be fully densified. The
#'   returned matrix is always dense.
#'
#' @param feature_matrix A numeric matrix representing features. It can be in
#'   `V_p x C` format (parcels x features/conditions) or `C x V_p` format.
#'   The function will attempt to orient it correctly based on `U_basis`.
#'   If `feature_matrix` is a sparse `Matrix`, operations use sparse methods
#'   when possible (coercion occurs only for small intermediate matrices).
#'   Very large sparse inputs (`prod(dim(.)) > 1e7`) trigger an error to avoid
#'   accidental densification.
#' @param U_basis A numeric matrix representing the spectral basis, typically with
#'   dimensions `V_p x k_dims_basis` (parcels x basis dimensions).
#'   Sparse `U_basis` inputs are handled with `Matrix` methods. Extremely large
#'   sparse bases also trigger an error to prevent densification.
#' @param tol_orthonormal A numeric tolerance to check for orthonormality of `U_basis`.
#'   `max(abs(crossprod(U_basis) - I)) / k_dims_basis` is compared against this tolerance.
#'   Default is `sqrt(.Machine$double.eps)`.
#' @param assume_orthonormal Logical. If `TRUE`, `U_basis` is assumed to be orthonormal,
#'   and the check (including rank check and `UtU` calculation) is skipped, 
#'   directly using the faster projection `crossprod(U_basis, features)`.
#'   Default is `FALSE`.
#'
#' @return A numeric (dense) matrix containing the projected features.
#'   If `feature_matrix` was oriented to `V_p x C` (parcels x features/conditions),
#'   the output will be `k_dims_basis x C` (basis dimensions x features/conditions).
#'   If `feature_matrix` was `C x V_p` and was transposed for projection, the output
#'   will be `C x k_dims_basis` to maintain the original feature orientation as rows.
#'
#' @export
#' @examples
#' V_p <- 10
#' k_dims <- 3
#' C <- 5
#' U_basis_ortho <- svd(matrix(rnorm(V_p * k_dims), V_p, k_dims))$u
#' U_basis_non_ortho <- matrix(rnorm(V_p * k_dims), V_p, k_dims)
#'
#' # Feature matrix: V_p x C (parcels x conditions)
#' features1 <- matrix(rnorm(V_p * C), V_p, C)
#' proj1_ortho <- project_features_to_spectral_space(features1, U_basis_ortho)
#' print(dim(proj1_ortho))
#' proj1_non_ortho <- project_features_to_spectral_space(features1, U_basis_non_ortho)
#' print(dim(proj1_non_ortho))
#'
#' # Sparse input example
#' features_sp <- Matrix::rsparsematrix(V_p, C, 0.2)
#' proj_sp <- project_features_to_spectral_space(features_sp, U_basis_ortho)
#' print(class(proj_sp))
#'
#' # Ambiguous square matrix error
#' V_p_sq <- 7
#' U_basis_sq <- matrix(rnorm(V_p_sq*k_dims), V_p_sq, k_dims)
#' features_sq_ambiguous <- matrix(rnorm(V_p_sq*V_p_sq), V_p_sq, V_p_sq)
#' try(project_features_to_spectral_space(features_sq_ambiguous, U_basis_sq))
#' 
#' # Rank deficient U_basis warning
#' U_rank_def <- matrix(rnorm(V_p * k_dims), V_p, k_dims)
#' if (k_dims > 1) U_rank_def[, k_dims] <- U_rank_def[, 1] # Make last col same as first
#' try(project_features_to_spectral_space(features1, U_rank_def))
#'
project_features_to_spectral_space <- function(feature_matrix, U_basis, 
                                             tol_orthonormal = sqrt(.Machine$double.eps),
                                             assume_orthonormal = FALSE) {
  # Sparse matrix handling with size guard
  sparse_size_threshold <- 1e7 # Heuristic for "large"
  feature_is_sparse <- inherits(feature_matrix, "sparseMatrix")
  basis_is_sparse <- inherits(U_basis, "sparseMatrix")
  if (feature_is_sparse && prod(dim(feature_matrix)) > sparse_size_threshold) {
    stop("feature_matrix is a large sparse Matrix; please convert to dense explicitly if memory allows, or use a sparse-aware projection method.")
  }
  if (basis_is_sparse && prod(dim(U_basis)) > sparse_size_threshold) {
    stop("U_basis is a large sparse Matrix; please convert to dense explicitly if memory allows, or use a sparse-aware projection method.")
  }

  if (!(is.matrix(feature_matrix) || inherits(feature_matrix, "Matrix")) ||
      !is.numeric(feature_matrix)) {
    stop("feature_matrix must be a numeric matrix or Matrix object.")
  }
  if (!(is.matrix(U_basis) || inherits(U_basis, "Matrix")) ||
      !is.numeric(U_basis)) {
    stop("U_basis must be a numeric matrix or Matrix object.")
  }

  V_p_basis <- nrow(U_basis)
  k_dims_basis <- ncol(U_basis)

  if (V_p_basis == 0) {
    warning("U_basis has 0 rows (V_p = 0). Projection is ill-defined or trivial.")
    if (nrow(feature_matrix) == 0) {
        return(matrix(0, nrow = k_dims_basis, ncol = ncol(feature_matrix)))
    } else if (ncol(feature_matrix) == 0) {
        return(matrix(0, nrow = nrow(feature_matrix), ncol = k_dims_basis))
    } else {
        return(matrix(0,0,0))
    }
  }
  
  if (k_dims_basis == 0) {
    warning("U_basis has 0 columns (k_dims_basis = 0). Resulting projection will have 0 dimensions.")
    if (nrow(feature_matrix) == V_p_basis) {
        return(matrix(0, nrow = 0, ncol = ncol(feature_matrix)))
    } else if (ncol(feature_matrix) == V_p_basis) {
        return(matrix(0, nrow = nrow(feature_matrix), ncol = 0))
    } else { 
        stop(sprintf("Dimensions of feature_matrix (%d x %d) are incompatible with U_basis (%d x %d) when k_dims_basis is 0.", 
                     nrow(feature_matrix), ncol(feature_matrix), V_p_basis, k_dims_basis))
    }
  }

  fm_oriented <- NULL
  input_was_transposed <- FALSE
  
  is_rows_match_Vp <- nrow(feature_matrix) == V_p_basis
  is_cols_match_Vp <- ncol(feature_matrix) == V_p_basis

  if (is_rows_match_Vp && is_cols_match_Vp) {
    stop(sprintf("feature_matrix dimensions (%d x %d) are ambiguous as both match V_p_basis (%d). Please ensure feature_matrix is oriented as V_p x C (parcels x features/conditions). If unsure, or if it is C x V_p where C == V_p, please transpose it first before calling.",
                 nrow(feature_matrix), ncol(feature_matrix), V_p_basis))
  } else if (is_rows_match_Vp) {
    fm_oriented <- feature_matrix 
    input_was_transposed <- FALSE
  } else if (is_cols_match_Vp) {
    fm_oriented <- t(feature_matrix) 
    input_was_transposed <- TRUE
  } else {
    stop(sprintf("feature_matrix dimensions (%d x %d) are incompatible with U_basis dimensions (%d x %d). One dimension of feature_matrix must match nrow(U_basis).",
                 nrow(feature_matrix), ncol(feature_matrix), V_p_basis, k_dims_basis))
  }
  
  if (ncol(fm_oriented) == 0) { 
      warning("Oriented feature matrix has 0 columns (C=0). Resulting projection will have 0 columns.")
      if (input_was_transposed) {
          return(matrix(0, nrow = 0, ncol = k_dims_basis))
      } else {
          return(matrix(0, nrow = k_dims_basis, ncol = 0))
      }
  }

  is_orthonormal <- FALSE
  UtU <- NULL # Define UtU outside so it's available for qr.solve if needed

  if (assume_orthonormal) {
    is_orthonormal <- TRUE
  } else {
    # Check U_basis rank first
    qr_U <- qr(U_basis)
    if (qr_U$rank < k_dims_basis) {
        warning(sprintf("U_basis appears rank deficient (rank %d / %d dimensions). Projection results may be unreliable or reflect a lower-dimensional space.",
                        qr_U$rank, k_dims_basis))
    }
    UtU <- crossprod(U_basis)
    Id_k <- diag(k_dims_basis)
    max_dev <- max(abs(as.matrix(UtU) - Id_k))
    # Audit suggestion: Scale tolerance by k_dims_basis
    is_orthonormal <- max_dev < (tol_orthonormal * k_dims_basis)
  }

  projected_features_internal <- NULL 
  if (is_orthonormal) {
    projected_features_internal <- crossprod(U_basis, fm_oriented)
  } else {
    if (is.null(UtU)) UtU <- crossprod(U_basis) # Should be calculated if !assume_orthonormal

    projected_features_internal <- tryCatch({
      # Using UtU in qr.solve for (U^T U)^-1 U^T X form
      qr.solve(as.matrix(UtU), as.matrix(crossprod(U_basis, fm_oriented)),
               tol = tol_orthonormal * k_dims_basis)
    }, error = function(e) {
      stop(paste("Error in qr.solve, likely UtU is singular or ill-conditioned (U_basis may be rank-deficient, or k_dims_basis > V_p_basis and U_basis not full rank). Original error:", e$message))
    })
  }

  if (input_was_transposed) {
    return(as.matrix(t(projected_features_internal)))
  } else {
    return(as.matrix(projected_features_internal))
  }
}

#' @title Residualize Matrix on Subspace
#' @description Residualizes the rows of `matrix_to_residualize` with respect to
#'   the row space of `subspace_basis_matrix`. This is equivalent to projecting
#'   each column of `t(matrix_to_residualize)` onto the column space of
#'   `t(subspace_basis_matrix)` and taking the residuals. Sparse inputs are
#'   processed with `Matrix` routines and the result is returned as a dense
#'   matrix.
#'
#' @param matrix_to_residualize A numeric matrix (e.g., `C x k`) whose rows
#'   are to be residualized. Sparse `Matrix` inputs are handled without full
#'   densification. A warning is issued if the matrix is very large (`C*k > 1e7`).
#' @param subspace_basis_matrix A numeric matrix (e.g., `m x k`) whose rows
#'   span the subspace to project out. Sparse matrices are accepted and processed
#'   with `Matrix` linear algebra.
#'
#' @return A numeric matrix with the same dimensions as `matrix_to_residualize`,
#'   containing the residualized rows.
#'
#' @export
#' @examples
#' # Z_i: 3 conditions, 5-dimensional spectral space (3x5)
#' Z_i <- matrix(rnorm(15), nrow = 3, ncol = 5)
#' # A_parc_i: 2 anchor parcels, 5-dimensional spectral space (2x5)
#' A_parc_i <- matrix(rnorm(10), nrow = 2, ncol = 5)
#' Z_i_res <- residualize_matrix_on_subspace(Z_i, A_parc_i)
#' print(dim(Z_i_res))
#' if (nrow(A_parc_i) > 0 && ncol(A_parc_i) > 0 && qr(t(A_parc_i))$rank > 0) {
#'   # Test orthogonality: Z_i_res %*% t(A_parc_i) should be near zero
#'   print(round(Z_i_res %*% t(A_parc_i), 10))
#' }
#' # Sparse example
#' Z_sp <- Matrix::rsparsematrix(3, 5, 0.3)
#' A_sp <- Matrix::rsparsematrix(2, 5, 0.3)
#' residualize_matrix_on_subspace(Z_sp, A_sp)
#'
residualize_matrix_on_subspace <- function(matrix_to_residualize, subspace_basis_matrix) {
  large_matrix_threshold <- 1e7 # Heuristic for C*k

  matrix_is_sparse <- inherits(matrix_to_residualize, "sparseMatrix")
  basis_is_sparse <- inherits(subspace_basis_matrix, "sparseMatrix")

  if (!(is.matrix(matrix_to_residualize) || inherits(matrix_to_residualize, "Matrix")) ||
      !is.numeric(matrix_to_residualize)) {
    stop("matrix_to_residualize must be a numeric matrix or Matrix object.")
  }
  if (!(is.matrix(subspace_basis_matrix) || inherits(subspace_basis_matrix, "Matrix")) ||
      !is.numeric(subspace_basis_matrix)) {
    stop("subspace_basis_matrix must be a numeric matrix or Matrix object.")
  }

  k_dim_Y <- ncol(matrix_to_residualize)
  k_dim_X <- ncol(subspace_basis_matrix)
  C_dim_Y <- nrow(matrix_to_residualize)
  m_dim_X <- nrow(subspace_basis_matrix)
  
  if (prod(dim(matrix_to_residualize)) > large_matrix_threshold) {
      warning("matrix_to_residualize is large; residualization involves transposes and may be memory intensive.")
  }

  if (k_dim_Y != k_dim_X) {
    stop(sprintf("Number of columns must match: matrix_to_residualize has %d, subspace_basis_matrix has %d.",
                 k_dim_Y, k_dim_X))
  }

  if (k_dim_Y == 0) { 
    return(matrix(0, nrow = C_dim_Y, ncol = 0))
  }
  if (m_dim_X == 0) { 
    return(matrix_to_residualize)
  }
  if (C_dim_Y == 0) { 
    return(matrix(0, nrow = 0, ncol = k_dim_Y))
  }

  Y_eff <- t(matrix_to_residualize)
  X_eff <- t(subspace_basis_matrix) 
  
  # Audit suggestion: NA/Inf check
  if (any(!is.finite(X_eff))) stop("subspace_basis_matrix (after transpose) contains non-finite values.")
  if (any(!is.finite(Y_eff))) stop("matrix_to_residualize (after transpose) contains non-finite values.")
  
  qr_X_eff <- Matrix::qr(X_eff)
  
  if (qr_X_eff$rank == 0) {
      return(matrix_to_residualize) 
  }
  
  if (qr_X_eff$rank < ncol(X_eff)) {
      warning(sprintf("subspace_basis_matrix (after transpose, %d x %d) appears rank deficient (rank %d < %d columns). Projection is onto a lower-dimensional subspace.",
              nrow(X_eff), ncol(X_eff), qr_X_eff$rank, ncol(X_eff)))
  }
  
  # TODO: Consider orientation switch for big k_dim_Y > C_dim_Y if performance critical
  residuals_eff <- Matrix::qr.resid(qr_X_eff, Y_eff)

  return(as.matrix(t(residuals_eff)))
}

#' @title Build Augmented Anchor Matrix
#' @description Combines parcel anchors and projected task features into a single
#'   augmented anchor matrix for Procrustes alignment.
#'   If inputs have different numeric types (e.g. integer and double), the result
#'   is coerced per R's default rules (usually to the more general type).
#'
#' @param A_parcel_anchors A numeric matrix representing parcel anchors,
#'   typically with dimensions `m_parcels x k_dims`.
#' @param Z_task_features_projected A numeric matrix representing projected task
#'   features (e.g., condition means in the spectral space), typically with
#'   dimensions `m_task_features x k_dims`. Can be `NULL` or have 0 rows if
#'   no task features are to be added.
#'
#' @return An augmented numeric matrix with dimensions
#'   `(m_parcels + m_task_features) x k_dims`. 
#'   Row names are combined from inputs (made unique with `make.unique`). 
#'   Column names are taken from `A_parcel_anchors` if it has them and rows; 
#'   otherwise from `Z_task_features_projected`. An error is thrown if both have
#'   differing non-NULL column names.
#'   If `Z_task_features_projected` is `NULL` or has 0 rows, `A_parcel_anchors` is returned.
#'   If `A_parcel_anchors` has 0 rows and `Z_task_features_projected` is valid,
#'   `Z_task_features_projected` is returned.
#'
#' @export
#' @examples
#' m_parcels <- 5; m_task <- 3; k_dims <- 4
#' A_p <- matrix(rnorm(m_parcels*k_dims), m_parcels, k_dims, dimnames=list(paste0("p",1:m_parcels), paste0("k",1:k_dims)))
#' Z_t <- matrix(rnorm(m_task*k_dims), m_task, k_dims, dimnames=list(paste0("t",1:m_task), paste0("k",1:k_dims)))
#' A_aug <- build_augmented_anchor_matrix(A_p, Z_t)
#' print(dimnames(A_aug))
#'
#' # Differing colnames should error
#' Z_t_bad_colnames <- Z_t; colnames(Z_t_bad_colnames) <- paste0("dim",1:k_dims)
#' try(build_augmented_anchor_matrix(A_p, Z_t_bad_colnames))
#'
build_augmented_anchor_matrix <- function(A_parcel_anchors, Z_task_features_projected) {

  if (!is.matrix(A_parcel_anchors) || !is.numeric(A_parcel_anchors)) {
    stop("A_parcel_anchors must be a numeric matrix.")
  }

  has_A_parc_rows <- nrow(A_parcel_anchors) > 0
  has_Z_task_rows <- !is.null(Z_task_features_projected) && 
                       is.matrix(Z_task_features_projected) && 
                       nrow(Z_task_features_projected) > 0

  A_parc_rownames <- if (has_A_parc_rows) rownames(A_parcel_anchors) else NULL
  A_parc_colnames <- if (has_A_parc_rows) colnames(A_parcel_anchors) else NULL # Could be NULL even if rows exist
  
  Z_task_rownames <- NULL
  Z_task_colnames <- NULL
  if (has_Z_task_rows) {
    if (!is.numeric(Z_task_features_projected)) {
      stop("Z_task_features_projected, if provided and non-empty, must be a numeric matrix.")
    }
    Z_task_rownames <- rownames(Z_task_features_projected)
    Z_task_colnames <- colnames(Z_task_features_projected) # Could be NULL
    
    if (has_A_parc_rows && (ncol(A_parcel_anchors) != ncol(Z_task_features_projected))) {
        stop(sprintf("Number of columns must match: A_parcel_anchors has %d, Z_task_features_projected has %d.",
                   ncol(A_parcel_anchors), ncol(Z_task_features_projected)))
    }
    # Audit suggestion: Strict colname check
    if (!is.null(A_parc_colnames) && !is.null(Z_task_colnames) && !identical(A_parc_colnames, Z_task_colnames)) {
         stop("Column names/order differ between A_parcel_anchors and Z_task_features_projected; please align before binding.")
    }
  }

  if (!has_Z_task_rows) {
    return(A_parcel_anchors) 
  }
  if (!has_A_parc_rows) {
    return(Z_task_features_projected) 
  }
  
  result <- rbind(A_parcel_anchors, Z_task_features_projected)
  
  # Ensure A_parc_rownames and Z_task_rownames are character vectors
  if (is.null(A_parc_rownames)) {
    A_parc_rownames <- paste0("P", 1:nrow(A_parcel_anchors))
  }
  
  if (is.null(Z_task_rownames)) {
    Z_task_rownames <- paste0("T", 1:nrow(Z_task_features_projected))
  }
  
  # Combine row names and make unique
  final_rownames <- make.unique(c(A_parc_rownames, Z_task_rownames), sep = "_")
  if (length(final_rownames) == nrow(result)) { # make.unique can return shorter if inputs were NULL
    rownames(result) <- final_rownames
  }
  
  final_colnames <- A_parc_colnames 
  if (is.null(final_colnames) && !is.null(Z_task_colnames)) {
      final_colnames <- Z_task_colnames
  }
  if (!is.null(final_colnames) && length(final_colnames) == ncol(result)) {
      colnames(result) <- final_colnames
  }
  
  return(result)
} 