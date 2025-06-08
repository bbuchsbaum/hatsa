#' Basic SVD-based Procrustes rotation solver
#'
#' Computes the cross-product `t(A) %*% T`, performs SVD and applies the
#' standard determinant correction so that the returned rotation lies in
#' `SO(k)`.
#'
#' @param A Numeric matrix.
#' @param T Numeric matrix of the same dimensions as `A`.
#' @return Rotation matrix with `det(R) = 1`.
#' @keywords internal
procrustes_rotation_basic <- function(A, T) {
  M <- crossprod(A, T)
  sv <- svd(M)
  U <- sv$u
  V <- sv$v
  R_raw <- V %*% t(U)
  sign_det <- sign(prod(diag(qr(R_raw)$qr)))
  R <- R_raw
  if (sign_det < 0) {
    j_min <- which.min(sv$d)
    V_corr <- V
    V_corr[, j_min] <- -V_corr[, j_min]
    R <- V_corr %*% t(U)
  }
  R
}

#' Solve Orthogonal Procrustes Problem for `R_i` in SO(k)
#'
#' Finds `R_i in SO(k)` that best aligns `A_orig_subj_anchor` to `T_anchor_group`.
#' This is the standard unweighted version.
#'
#' @param A_orig_subj_anchor Numeric matrix (`m x k`), subject's anchor sketch.
#' @param T_anchor_group Numeric matrix (`m x k`), group anchor template.
#' @return Orthogonal rotation matrix `R_i` (`k x k`) with `det(R_i) = 1`.
#' @keywords internal
solve_procrustes_rotation <- function(A_orig_subj_anchor, T_anchor_group) {
  if (!is.matrix(A_orig_subj_anchor) || !is.numeric(A_orig_subj_anchor)) {
      stop("A_orig_subj_anchor must be a numeric matrix")
  }
  if (!is.matrix(T_anchor_group) || !is.numeric(T_anchor_group)) {
      stop("T_anchor_group must be a numeric matrix")
  }
  if (nrow(A_orig_subj_anchor) != nrow(T_anchor_group) ||
      ncol(A_orig_subj_anchor) != ncol(T_anchor_group)) {
      stop("A_orig_subj_anchor and T_anchor_group must have matching dimensions")
  }
  if (!all(is.finite(A_orig_subj_anchor))) {
      stop("A_orig_subj_anchor contains non-finite values")
  }
  if (!all(is.finite(T_anchor_group))) {
      stop("T_anchor_group contains non-finite values")
  }

  k_dim <- ncol(A_orig_subj_anchor)
  if (k_dim == 0) return(matrix(0,0,0))

  M <- crossprod(A_orig_subj_anchor, T_anchor_group)

  if (all(abs(M) < 1e-14)) {
      warning("Cross-product matrix M in solve_procrustes_rotation is near zero; rotation is ill-defined. Returning identity.")
      return(diag(k_dim))
  }

###<<<<<<< codex/create-procrustes_rotation_basic-helper
  R_i <- procrustes_rotation_basic(A_orig_subj_anchor, T_anchor_group)

##=======
  svd_M <- svd(M) 
  U_svd <- svd_M$u
  V_svd <- svd_M$v 
  
  R_raw <- V_svd %*% base::t(U_svd) 
  
  # Fast reflection fix using determinant sign from SVD
  # sign(det(V) * det(U)) computes the sign of det(R_raw)
  sign_det <- sign(det(svd_M$v) * det(svd_M$u))
  
  R_i <- R_raw
  if (sign_det < 0) {
      # Audit patch: Flip column corresponding to the *smallest* singular value
      j_min_sv <- which.min(svd_M$d)
      V_svd_corrected <- V_svd
      V_svd_corrected[, j_min_sv] <- -V_svd_corrected[, j_min_sv]
      R_i <- V_svd_corrected %*% base::t(U_svd)
      
      # Optional sanity check (can be removed in production)
      # if (abs(det(R_i) - 1.0) > 1e-6) {
      #   warning("Determinant correction failed?")
      # }
  }
  
  # The previous check for det != +/- 1 is removed as the determinant sign check is sufficient
  # and handles the reflection properly. Issues with the matrix M itself
  # (e.g. near singularity) might still lead to unstable rotations, but the 
  # rotation matrix R_i returned will have det = +1.
  
###>>>>>>> main
  return(R_i)
}

#' Perform Generalized Procrustes Analysis (GPA) refinement iterations
#'
#' Iteratively refines subject rotations and a group anchor template.
#' Can use weighted Procrustes if task rows are present and omega_mode is specified.
#' If n_refine = 0, computes initial rotations against the mean template but does not update the template.
#'
#' @param A_originals_list List of `(m_parcels + m_tasks) x k` anchor matrices for `N` subjects.
#' @param n_refine Integer, number of refinement iterations.
#' @param k Integer, spectral rank (number of columns in anchor matrices).
#' @param m_parcel_rows Integer, number of rows corresponding to parcel anchors.
#'   Required if `m_task_rows > 0` for weighted Procrustes.
#' @param m_task_rows Integer, number of rows corresponding to task anchors. Default 0.
#'   If > 0, `solve_procrustes_rotation_weighted` is used.
#' @param omega_mode Character string, mode for `solve_procrustes_rotation_weighted`.
#'   Default `"fixed"`. Ignored if `m_task_rows == 0`.
#' @param fixed_omega_weights List, weights for `"fixed"` mode. Passed to weighted solver.
#'   Default `list(parcel = 1.0, condition = 0.5)` is handled by the weighted solver if NULL.
#' @param reliability_scores_list List (parallel to `A_originals_list`), each element a
#'   numeric vector of reliability scores for task anchors. Used if `omega_mode == "adaptive"`
#'   and `m_task_rows > 0`. If `NULL`, adaptive mode in solver may use default behavior.
#' @param scale_omega_trace Logical, passed to weighted solver. Default `TRUE`.
#'
#' @return List: `R_final_list` (list of `k x k` rotation matrices),
#'         `T_anchor_final` (final `(m_parcels+m_tasks) x k` group anchor template).
#' @keywords internal
perform_gpa_refinement <- function(A_originals_list, n_refine, k,
                                   m_parcel_rows = NULL, m_task_rows = 0,
                                   omega_mode = "fixed",
                                   fixed_omega_weights = NULL,
                                   reliability_scores_list = NULL,
                                   scale_omega_trace = TRUE) {
  num_subjects <- length(A_originals_list)

  omega_mode <- match.arg(omega_mode, c("fixed", "adaptive"))

  if (length(n_refine) != 1L || is.na(n_refine) || n_refine < 0 ||
      n_refine %% 1 != 0) {
    stop("n_refine must be a non-negative integer")
  }
  
  if (is.null(m_parcel_rows)) {
    first_valid_A <- NULL
    for (subj_A in A_originals_list) {
        if (!is.null(subj_A) && is.matrix(subj_A) && ncol(subj_A) == k) {
            first_valid_A <- subj_A
            break
        }
    }
    if (!is.null(first_valid_A)) {
        m_total_rows_inferred <- nrow(first_valid_A)
    } else {
        m_total_rows_inferred <- 0 
    }
    if (m_task_rows > 0) stop("m_parcel_rows must be specified if m_task_rows > 0.")
    m_parcel_rows <- m_total_rows_inferred 
  }
  
  m_total_rows <- m_parcel_rows + m_task_rows

  dimension_check <- vapply(A_originals_list, function(A) {
    is.null(A) || (is.matrix(A) && all(dim(A) == c(m_total_rows, k)))
  }, logical(1))
  if (any(!dimension_check)) {
    stop("Each non-NULL element of A_originals_list must have dimensions (m_parcel_rows + m_task_rows) x k")
  }

  if (k == 0) { 
      T_anchor <- matrix(0, nrow = m_total_rows, ncol = 0)
      R_final_list <- replicate(num_subjects, matrix(0,0,0), simplify=FALSE)
      return(list(R_final_list = R_final_list, T_anchor_final = T_anchor))
  }
  if (m_total_rows == 0 && k > 0){
      T_anchor <- matrix(0, nrow = 0, ncol = k)
      R_final_list <- replicate(num_subjects, diag(k), simplify = FALSE)
      return(list(R_final_list = R_final_list, T_anchor_final = T_anchor))
  }

  valid_A_originals_indices <- sapply(A_originals_list, function(A) {
    !is.null(A) && is.matrix(A) && nrow(A) == m_total_rows && ncol(A) == k
  })
  
  valid_A_originals <- A_originals_list[valid_A_originals_indices]

  if (length(valid_A_originals) > 0) {
      T_anchor <- Reduce(`+`, valid_A_originals) / length(valid_A_originals)
  } else { 
      T_anchor <- matrix(0, nrow = m_total_rows, ncol = k)
      warning("No valid anchor matrices found for initialization. GPA template starts at zero.")
  }

  R_iter_list <- replicate(num_subjects, diag(k), simplify = FALSE)

  # Audit patch: Handle n_refine == 0 by computing initial rotations only
  if (n_refine == 0L) {
      for (i in seq_len(num_subjects)) {
          Ai <- A_originals_list[[i]]
          # Check validity again, although indices were computed above
          if (!valid_A_originals_indices[i]) next 
          
          subj_reliability_scores <- if (!is.null(reliability_scores_list) && length(reliability_scores_list) >= i) reliability_scores_list[[i]] else NULL
          
          if (m_task_rows > 0) {
                 R_iter_list[[i]] <- solve_procrustes_rotation_weighted(
                                       Ai, T_anchor,
                                       m_parcel_rows, m_task_rows,
                                       omega_mode, fixed_omega_weights,
                                       subj_reliability_scores, # Use looked-up value
                                       scale_omega_trace)
          } else {
                 R_iter_list[[i]] <- solve_procrustes_rotation(Ai, T_anchor)
          }
      }
      return(list(R_final_list = R_iter_list, T_anchor_final = T_anchor))
  }
  
  # --- Main Refinement Loop (n_refine > 0) --- 
  for (iter_num in 1:n_refine) {
    rotated_anchors_sum_for_update <- matrix(0, nrow = m_total_rows, ncol = k)
    active_subjects_count_for_update <- 0

    for (i in 1:num_subjects) {
      A_orig_subj <- A_originals_list[[i]]
      
      if (valid_A_originals_indices[i]) { # Use pre-calculated validity
        current_R_i <- NULL
        subj_reliability_scores <- if (!is.null(reliability_scores_list) && length(reliability_scores_list) >= i) reliability_scores_list[[i]] else NULL
          
        if (m_task_rows > 0) {
          current_R_i <- solve_procrustes_rotation_weighted(
                                A_source = A_orig_subj, 
                                T_target = T_anchor, 
                                m_parcel_rows = m_parcel_rows, 
                                m_task_rows = m_task_rows, 
                                omega_mode = omega_mode, 
                                fixed_omega_weights = fixed_omega_weights,
                                reliability_scores = subj_reliability_scores,
                                scale_omega_trace = scale_omega_trace
                              )
        } else {
          current_R_i <- solve_procrustes_rotation(A_orig_subj, T_anchor)
        }
        R_iter_list[[i]] <- current_R_i
        
        rotated_anchors_sum_for_update <- rotated_anchors_sum_for_update + (A_orig_subj %*% current_R_i)
        active_subjects_count_for_update <- active_subjects_count_for_update + 1
      } else {
        # Subject data invalid or missing, R_iter_list[[i]] remains (identity from init or prev iter)
      }
    }
    
    if (active_subjects_count_for_update > 0) {
        T_anchor <- rotated_anchors_sum_for_update / active_subjects_count_for_update
    } else if (m_total_rows > 0 || k > 0) { 
        T_anchor <- matrix(0, nrow = m_total_rows, ncol = k) 
        warning(sprintf("No valid subjects contributed to GPA template update in iteration %d. Template reset to zeros.", iter_num))
    }
  }
  
  return(list(R_final_list = R_iter_list, T_anchor_final = T_anchor))
}

#' Perform Geometric Generalized Procrustes Analysis (Geo-GPA) on SO(k)
#'
#' Refines a set of rotation matrices and a group anchor template using either
#' an SVD-based iterative Procrustes approach (default) or a variant where the
#' template update is oriented by the Frechet mean of rotations.
#' The SVD mode minimizes the Frobenius distance between rotated subject anchors
#' and the template. The Riemannian mode aims to make the template's orientation
#' consistent with the Frechet mean of subject rotations, which implicitly relates
#' to minimizing intrinsic SO(k) geodesic distances for the rotations' central tendency.
#'
#' @param A_originals_list A list of subject-specific anchor matrices (m_rows x k).
#'   Each matrix `A_i` contains the rows from `U_original_list[[i]]` corresponding
#'   to the `unique_anchor_indices`.
#' @param n_refine Integer, number of GPA refinement iterations. Default: 10.
#' @param k Integer, the dimensionality of the sketch space (number of columns in A_i and R_i).
#' @param m_rows Integer, number of anchor features (number of rows in A_i).
#' @param tol Numeric, tolerance for convergence based on the relative Frobenius norm
#'   change in the template `T_template`. Default: 1e-7.
#' @param rotation_mode Character string, one of `"svd"` (default) or `"riemannian"`.
#'   - `"svd"`: Uses SVD to find the closest SO(k) rotation for each subject to align
#'     `A_i` with `T_template`. `T_template` is the Euclidean mean of `A_i R_i`.
#'     This mode minimizes the sum of squared Frobenius distances `sum(||A_i R_i - T||_F^2)`.
#'   - `"riemannian"`: Individual rotations `R_i` are updated as in `"svd"` mode.
#'     However, `T_template` is updated by first computing the Frechet mean (`R_bar`)
#'     of the current `R_list`, then averaging `A_i R_i R_bar^T` (configurations aligned
#'     to `R_bar`'s frame), and finally rotating this average back by `R_bar`.
#'     This mode seeks to make the template orientation consistent with the intrinsic
#'     mean of the rotation ensemble.
#' @param frechet_mean_options A list of options to pass to `hatsa::frechet_mean_so_k`,
#'   used if `rotation_mode = "riemannian"` or for the final `R_bar_final` computation.
#'   Example: `list(max_iter = 20, tol = 1e-5)`. Defaults are used if not provided.
#' @param verbose Logical, if TRUE, prints progress messages. Default TRUE.
#' @param initial_R_list Optional list of initial k x k rotation matrices. If NULL,
#'   identity matrices are used.
#' @return A list containing:
#'   - `R_final_list`: List of final subject-specific rotation matrices (k x k).
#'   - `T_anchor_final`: The final group anchor template matrix (m_rows x k).
#'   - `R_bar_final`: The Fréchet mean of the final `R_final_list`.
#' @importFrom stats svd
#' @keywords internal
perform_geometric_gpa_refinement <- function(A_originals_list,
                                             n_refine = 10,
                                             k,
                                             m_rows,
                                             tol = 1e-7,
                                             rotation_mode = c("svd", "riemannian"),
                                             frechet_mean_options = list(),
                                             verbose = TRUE,
                                             initial_R_list = NULL) {

  rotation_mode <- match.arg(rotation_mode)
  N <- length(A_originals_list)

  if (N == 0) {
    warning("A_originals_list is empty. Cannot perform GPA.")
    return(list(R_final_list = list(), T_anchor_final = matrix(NA, nrow = m_rows, ncol = k), R_bar_final = diag(k)))
  }

  # Initialize R_list
  R_list <- if (!is.null(initial_R_list) && length(initial_R_list) == N) {
    initial_R_list
  } else {
    replicate(N, diag(k), simplify = FALSE)
  }
  
  # Ensure initial R_list contains SO(k) matrices
  for(i in 1:N){
    if (!is.null(R_list[[i]]) && all(dim(R_list[[i]]) == c(k,k))) {
        R_list[[i]] <- procrustes_rotation_basic(R_list[[i]], diag(k))
    } else {
        R_list[[i]] <- diag(k) # Fallback if invalid initial matrix
    }
  }

  # Initial template T_template
  T_sum_initial <- matrix(0, nrow = m_rows, ncol = k)
  valid_configs_count_initial <- 0
  for (i in 1:N) {
    if (!is.null(A_originals_list[[i]]) && !is.null(R_list[[i]])) {
      T_sum_initial <- T_sum_initial + (A_originals_list[[i]] %*% R_list[[i]])
      valid_configs_count_initial <- valid_configs_count_initial + 1
    }
  }
  if (valid_configs_count_initial == 0) stop("No valid initial configurations (A_i R_i) to compute initial template.")
  T_template <- T_sum_initial / valid_configs_count_initial
  
  R_bar_current_iter <- NULL # Used for Riemannian mode's initial_mean for frechet_mean_so_k

  for (iter in seq_len(n_refine)) {
    T_old <- T_template

    # --- 1. Update R_i for each subject --- 
    for (i in 1:N) {
      Ai <- A_originals_list[[i]]
      if (is.null(Ai) || !is.matrix(Ai) || nrow(Ai) != m_rows || ncol(Ai) != k) {
        if(verbose) message(sprintf("Skipping rotation update for subject %d due to invalid A_i.", i))
        R_list[[i]] <- diag(k) # Keep as identity or previous if A_i is problematic
        next
      }
      R_list[[i]] <- procrustes_rotation_basic(Ai, T_template)
    }

    # --- 2. Update T_template --- 
    if (rotation_mode == "svd") {
      T_sum <- matrix(0, nrow = m_rows, ncol = k)
      valid_configs_count <- 0
      for (i in 1:N) {
          if (!is.null(A_originals_list[[i]]) && !is.null(R_list[[i]])) {
            T_sum <- T_sum + (A_originals_list[[i]] %*% R_list[[i]])
            valid_configs_count <- valid_configs_count + 1
          }
      }
      if (valid_configs_count == 0) {
          warning("No valid configurations to update template in iteration ", iter, ". Using previous template.")
          T_template <- T_old
      } else {
          T_template <- T_sum / valid_configs_count
      }
    } else { # rotation_mode == "riemannian"
      # Compute R_bar = FrechetMean(R_list)
      current_frechet_opts <- list(
          k_dim = k,
          max_iter = frechet_mean_options$max_iter %||% 20, 
          tol = frechet_mean_options$tol %||% 1e-5, 
          initial_mean = frechet_mean_options$initial_mean %||% R_bar_current_iter %||% R_list[[1]],
          project_to_SOk = TRUE
      )
      if(is.null(current_frechet_opts$initial_mean)) current_frechet_opts$initial_mean <- diag(k)

      R_bar <- do.call(frechet_mean_so_k, c(list(R_list = R_list), current_frechet_opts))
      R_bar_current_iter <- R_bar # Save for next iteration's initial_mean hint
      
      T_sum_aligned <- matrix(0, nrow = m_rows, ncol = k)
      valid_configs_count_riem <- 0
      if (!is.null(R_bar) && all(dim(R_bar) == c(k,k))) {
            R_bar_t <- t(R_bar)
            for (i in 1:N) {
                if (!is.null(A_originals_list[[i]]) && !is.null(R_list[[i]])) {
                    T_sum_aligned <- T_sum_aligned + (A_originals_list[[i]] %*% R_list[[i]] %*% R_bar_t)
                    valid_configs_count_riem <- valid_configs_count_riem + 1
                }
            }
            if (valid_configs_count_riem == 0) {
                warning("No valid configurations for Riemannian template update in iter ", iter, ". Using previous template.")
                T_template <- T_old
            } else {
                T_template <- (T_sum_aligned / valid_configs_count_riem) %*% R_bar
            }
      } else {
          warning("R_bar computation failed or invalid in Riemannian mode iter ", iter, ". Using previous template.")
          T_template <- T_old
      }
    }

    # --- 3. Convergence Check --- 
    delta <- norm(T_template - T_old, type = "F") / max(1e-9, norm(T_old, type = "F"))
    if (verbose) {
      message(sprintf("Geo-GPA iter %d/%d (%s mode): Template change (Delta) = %.2e", 
                      iter, n_refine, rotation_mode, delta))
    }
    if (delta < tol) {
      if (verbose) message("Geo-GPA converged.")
      break
    }
    if (iter == n_refine && verbose) {
        message("Geo-GPA reached max iterations without specified convergence.")
    }
  }

  # --- Final R_bar (Fréchet mean of final R_list) --- 
  final_frechet_opts <- list(
      k_dim = k,
      max_iter = frechet_mean_options$max_iter %||% 50, # More iter for final
      tol = frechet_mean_options$tol %||% 1e-7,
      initial_mean = frechet_mean_options$initial_mean %||% R_bar_current_iter %||% (if(length(R_list)>0) R_list[[1]] else diag(k)),
      project_to_SOk = TRUE
  )
   if(is.null(final_frechet_opts$initial_mean)) final_frechet_opts$initial_mean <- diag(k)

  R_bar_final <- if (length(R_list) > 0) {
      do.call(frechet_mean_so_k, c(list(R_list = R_list), final_frechet_opts))
  } else {
      diag(k)
  }
  if (is.null(R_bar_final)) R_bar_final <- diag(k) # Fallback if frechet_mean_so_k fails

  return(list(R_final_list = R_list, 
              T_anchor_final = T_template, 
              R_bar_final = R_bar_final))
}
