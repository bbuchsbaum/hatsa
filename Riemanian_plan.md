# Plan for Riemannian Geometry Metrics in HATSA (Post M&B Proposal)

This document outlines the plan to implement Riemannian geometry-based metrics, particularly focusing on the Mitteroecker & Bookstein (M&B) distance, for analyzing HATSA outputs. This builds upon the "M&B Proposal" and aligns with the RGEOM series of tickets in `Riemannian-tickets.md`.

**Key Objectives:**

1.  **Core Riemannian Distance Function for SPD Matrices in R:** Key functions like `riemannian_distance_spd` (M&B, 2008; Log-Euclidean) and helpers (e.g., `.regularize_spd`) are to be implemented in `R/riemannian_geometry.R` (see RGEOM-001).
2.  **Helper for Covariance of Spectral Coefficients:** A function to compute `Cov_coeff_i = cov(U_aligned_i)`.
3.  **New S3 Methods for `hatsa_projector`:**
    *   `riemannian_distance_matrix_spd()`: Computes all pairwise Riemannian distances between subjects' SPD matrix representations.
    *   `riemannian_dispersion_spd()`: Computes dispersion around a central tendency (e.g., Fréchet mean) using the Riemannian metric.
4.  **Integration into Analysis Workflows:** Provide guidance on using these metrics (rather than directly in `summary()` methods for speed).
5.  **Considerations for Use in Anchor/Parameter Optimization:** Discuss how these new metrics can feed into optimization loops (e.g., for `k`, anchors).

This plan primarily details the implementation for RGEOM-004, building on foundational aspects that would be part of RGEOM-001.

---

**1. Core Riemannian Distance Function for SPD Matrices in R**

Key functions for computing Riemannian distances and related operations on SPD matrices are now located in `R/riemannian_geometry.R`. This includes `.regularize_spd` for ensuring matrices are positive definite and `riemannian_distance_spd` for calculating the Log-Euclidean distance between two SPD matrices (S1, S2). The Log-Euclidean distance is often cited as `d(S1, S2) = ||logm(S1) - logm(S2)||_F` and is equivalent to the Mitteroecker & Bookstein (2008) formulation `sqrt(sum(log(lambda_j(S2^-1 S1))^2))`.

Implementations of these, along with other functions for RGEOM-001 (such as `matrix_logm_spd`, `matrix_expm_spd`, `airm_distance`, `frechet_mean_spd`), will reside in `R/riemannian_geometry.R`.

```R
# Core SPD metric functions like .regularize_spd and riemannian_distance_spd
# are now located in R/riemannian_geometry.R
# See that file for implementation details.
# Further functions for RGEOM-001 will also be added there.
```
*Note on broader context (RGEOM-001): The `riemannian_distance_spd` function implements the Log-Euclidean distance. A comprehensive Riemannian toolkit, as envisioned in RGEOM-001, would also include functions like `matrix_logm_spd`, `matrix_expm_spd`, `airm_distance`, and robust Fréchet mean computations (`frechet_mean_spd`), which are foundational for diverse geometric analyses on SPD manifolds.* This note now primarily serves as a reminder of the scope of RGEOM-001, with implementations in the .R file.

---

**2. Helper for Covariance of Spectral Coefficients**

```R
# R/hatsa_projector.R or R/riemannian_metrics.R

#' Compute Covariance of Aligned Spectral Coefficients
#'
#' For a given subject's aligned spectral sketch `U_aligned_i` (V_p x k),
#' this function computes the `k x k` covariance matrix of these coefficients
#' across the V_p parcels.
#'
#' @param U_aligned_subject A numeric matrix (V_p x k) of aligned spectral coefficients.
#' @return A `k x k` covariance matrix. Returns a matrix of NAs if V_p <= 1 or k = 0.
#' @importFrom stats cov
#' @keywords internal
.compute_cov_spectral_coeffs <- function(U_aligned_subject) {
  if (!is.matrix(U_aligned_subject) || !is.numeric(U_aligned_subject)) {
    stop("U_aligned_subject must be a numeric matrix.")
  }
  V_p <- nrow(U_aligned_subject)
  k <- ncol(U_aligned_subject)

  if (k == 0) {
    return(matrix(NA_real_, nrow = 0, ncol = 0))
  }
  if (V_p <= 1) {
    warning("Cannot compute covariance with V_p <= 1. Returning NA matrix.")
    return(matrix(NA_real_, nrow = k, ncol = k))
  }
  return(stats::cov(U_aligned_subject))
}
```

---

**3. New S3 Methods for `hatsa_projector` (Fulfilling RGEOM-004)**

These S3 methods provide convenient wrappers for calculating pairwise distances and dispersion for subjects within a `hatsa_projector` object. They fulfill the requirements of RGEOM-004.

```R
# R/hatsa_projector.R (or a new R/riemannian_methods_hatsa.R)

#' Compute Pairwise Riemannian Distances Between Subject SPD Representations
#'
#' Calculates a matrix of Riemannian distances (Log-Euclidean metric)
#' between subjects, based on specified SPD matrix representations derived from
#' the HATSA object (e.g., covariance of aligned coefficients or various FC matrices).
#'
#' @param object A `hatsa_projector` or `task_hatsa_projector` object.
#' @param type Character string indicating the type of SPD representation to use.
#'   Commonly `"cov_coeffs"` (default) for the `k x k` covariance of aligned
#'   spectral coefficients. Other types (e.g., `"fc_conn"`, `"fc_task"`) would
#'   rely on `get_spd_representations` (from RGEOM-003) to provide the SPD matrices.
#' @param subject_data_list Optional. A list of subject time-series matrices
#'   (T_i x V_p). May be needed by `get_spd_representations` for certain `type` values.
#' @param spd_regularize_epsilon Epsilon for regularizing SPD matrices.
#' @param ... Additional arguments passed to `riemannian_distance_spd` or
#'   to `get_spd_representations` if used.
#' @return A symmetric `N_subjects x N_subjects` matrix of pairwise Riemannian distances.
#'   Diagonal will be zero. `NA` if a distance cannot be computed.
#' @export
riemannian_distance_matrix_spd <- function(object, ...) {
  UseMethod("riemannian_distance_matrix_spd")
}

#' @rdname riemannian_distance_matrix_spd
#' @export
riemannian_distance_matrix_spd.hatsa_projector <- function(object,
                                             type = "cov_coeffs",
                                             subject_data_list = NULL,
                                             spd_regularize_epsilon = 1e-6,
                                             k_conn_params = NULL,
                                             ...) {
  # type argument is now passed directly to get_spd_representations, which will validate it.
  # N_subjects <- object$parameters$N_subjects # N_subjects derived from spd_matrices_list length
  
  # Retrieve all SPD matrices first for all subjects in the object
  # get_spd_representations will handle subject_idx = NULL to get all.
  # Pass relevant params for specific types like fc_conn.
  spd_matrices_list_all_subjects <- get_spd_representations(object,
                                                            type = type,
                                                            subject_idx = NULL, # Get for all subjects
                                                            regularize_epsilon = spd_regularize_epsilon,
                                                            subject_data_list_for_fc = subject_data_list,
                                                            k_conn_params_for_fc = k_conn_params,
                                                            ...)
  
  # Filter out subjects for whom SPD matrix computation failed (NULL entries)
  valid_indices_in_list <- which(!sapply(spd_matrices_list_all_subjects, is.null))
  
  if (length(valid_indices_in_list) < length(spd_matrices_list_all_subjects)) {
    warning(sprintf("Could not obtain valid SPD matrices for %d out of %d subject(s). Distances involving them will be NA.", 
                    length(spd_matrices_list_all_subjects) - length(valid_indices_in_list), length(spd_matrices_list_all_subjects)))
  }
  
  # If no valid matrices or only one, distance matrix is trivial or not well-defined for pairs
  if (length(valid_indices_in_list) < 2) {
    num_total_subjects <- length(spd_matrices_list_all_subjects) # Original N_subjects if available, else from list length
    if (is.null(num_total_subjects) || num_total_subjects == 0) num_total_subjects <- object$parameters$N_subjects

    warning("Need at least 2 subjects with valid SPD matrices to compute a distance matrix.")
    # Return a matrix of NAs or 0s based on total subjects, consistent with original intent
    # The names of spd_matrices_list_all_subjects are original subject indices as characters
    dist_mat_final_size <- matrix(NA_real_, nrow = num_total_subjects, ncol = num_total_subjects)
    if (num_total_subjects > 0) {
        subj_names <- names(spd_matrices_list_all_subjects) 
        if(is.null(subj_names) || length(subj_names) != num_total_subjects) { # Fallback if names are not set as expected
            subj_names <- as.character(1:num_total_subjects)
        }
        rownames(dist_mat_final_size) <- colnames(dist_mat_final_size) <- subj_names
        diag(dist_mat_final_size) <- 0.0
    }
    return(dist_mat_final_size)
  }

  # Proceed with the subset of valid SPD matrices
  spd_matrices_valid_subset <- spd_matrices_list_all_subjects[valid_indices_in_list]
  num_valid_subjects <- length(spd_matrices_valid_subset)
  valid_subj_original_indices <- names(spd_matrices_valid_subset) # These are original indices as characters

  # Pairwise distance matrix for the valid subset
  dist_mat_subset <- matrix(NA_real_, nrow = num_valid_subjects, ncol = num_valid_subjects)
  rownames(dist_mat_subset) <- colnames(dist_mat_subset) <- valid_subj_original_indices
  diag(dist_mat_subset) <- 0.0
  
  for (i in 1:(num_valid_subjects - 1)) {
    idx1_name <- valid_subj_original_indices[i]
    S1 <- spd_matrices_valid_subset[[idx1_name]]
    
    for (j in (i + 1):num_valid_subjects) {
      idx2_name <- valid_subj_original_indices[j]
      S2 <- spd_matrices_valid_subset[[idx2_name]]
      
      current_dist <- tryCatch(
        riemannian_distance_spd(S1, S2, regularize_epsilon = spd_regularize_epsilon, ...),
        error = function(e) {
          warning(sprintf("Riemannian distance computation failed between subject %s and %s: %s", idx1_name, idx2_name, e$message))
          NA_real_
        }
      )
      dist_mat_subset[idx1_name, idx2_name] <- current_dist
      dist_mat_subset[idx2_name, idx1_name] <- current_dist
    }
  }
  
  # Now, create the full-sized distance matrix including NAs for failed subjects
  # This ensures the output matrix dimension matches N_subjects from the hatsa_object
  # and has correct row/col names corresponding to original subject indices.
  N_total_subjects_in_object <- object$parameters$N_subjects
  all_original_subj_indices_char <- as.character(1:N_total_subjects_in_object)
  
  final_dist_mat <- matrix(NA_real_, nrow = N_total_subjects_in_object, ncol = N_total_subjects_in_object)
  rownames(final_dist_mat) <- colnames(final_dist_mat) <- all_original_subj_indices_char
  diag(final_dist_mat) <- 0.0
  
  # Fill in the computed distances for valid subjects
  if (num_valid_subjects > 0) {
      final_dist_mat[valid_subj_original_indices, valid_subj_original_indices] <- dist_mat_subset
  }
  
  return(final_dist_mat)
}


#' Compute Riemannian Dispersion of SPD Matrices on the Manifold
#'
#' Calculates the dispersion of a set of SPD matrices around their Fréchet mean
#' using the Riemannian (Log-Euclidean) distance.
#'
#' @param object A `hatsa_projector` or `task_hatsa_projector` object.
#' @param type Character string, either `"cov_coeffs"` (default) or others like
#'   `"fc_conn"`, `"fc_task"` (relying on RGEOM-003 `get_spd_representations`).
#'   Specifies which matrices to compute dispersion for.
#' @param subject_data_list Optional, may be needed by `get_spd_representations`.
#' @param use_geometric_median Logical, if TRUE (and package `shapes` available), uses
#'   `shapes::mediancov` for Fréchet mean. If FALSE (default), uses a simple
#'   Euclidean mean of SPD matrices as an approximation.
#' @param ... Additional arguments passed to `riemannian_distance_spd` or `get_spd_representations`.
#' @return A list containing: `mean_spd_matrix`, `distances_to_mean`, `mean_dispersion`, `median_dispersion`.
#' @export
riemannian_dispersion_spd <- function(object, ...) {
  UseMethod("riemannian_dispersion_spd")
}

#' @rdname riemannian_dispersion_spd
#' @export
riemannian_dispersion_spd.hatsa_projector <- function(object,
                                          type = "cov_coeffs",
                                          subject_data_list = NULL,
                                          use_geometric_median = FALSE,
                                          k_conn_params = NULL,
                                          spd_regularize_epsilon = 1e-6,
                                          verbose = FALSE,
                                          ...) {
  # N_subjects is derived from the length of the list returned by get_spd_representations or object$parameters
  
  # Retrieve all SPD matrices first for all subjects in the object
  # Pass subject_data_list as subject_data_list_for_fc to get_spd_representations
  # Pass k_conn_params as k_conn_params_for_fc
  spd_matrices_list_all_subjects <- get_spd_representations(object,
                                                            type = type,
                                                            subject_idx = NULL, # Get for all subjects
                                                            regularize_epsilon = spd_regularize_epsilon,
                                                            subject_data_list_for_fc = subject_data_list,
                                                            k_conn_params_for_fc = k_conn_params,
                                                            ...) # Pass other args like specific fc params
  
  # Filter out subjects for whom SPD matrix computation failed (NULL entries)
  # Also filter out any matrices that might be all NA (though get_spd_representations should give NULL for failures)
  valid_spd_matrices_with_names <- Filter(function(m) !is.null(m) && !all(is.na(m)), spd_matrices_list_all_subjects)
  
  # If spd_matrices_list_all_subjects was named by get_spd_representations, valid_spd_matrices_with_names will retain them.
  # If not, names will be NULL.
  
  num_valid_subjects <- length(valid_spd_matrices_with_names)
  original_num_subjects <- object$parameters$N_subjects

  if (num_valid_subjects == 0) {
    warning("No valid SPD matrices found to compute dispersion.")
    return(list(mean_spd_matrix = NULL, 
                distances_to_mean = setNames(numeric(0), character(0)), # Ensure named empty numeric
                num_valid_subjects = 0,
                original_num_subjects = original_num_subjects,
                mean_dispersion = NA_real_, 
                median_dispersion = NA_real_))
  }
  
  # Handle case with only one valid subject differently
  if (num_valid_subjects == 1) {
      # If there was only one subject to begin with, or only one valid one from many
      single_matrix <- valid_spd_matrices_with_names[[1]]
      single_matrix_name <- names(valid_spd_matrices_with_names)[1]
      if (is.null(single_matrix_name)) single_matrix_name <- "1" # Default name if none

      # For a single matrix, mean is itself, distance to mean is 0.
      warning_msg <- if (original_num_subjects > 1) {
          "Only one valid SPD matrix found out of multiple subjects. Dispersion is 0 for this single matrix."
      } else {
          "Only one subject provided. Dispersion is 0."
      }
      warning(warning_msg)
       return(list(mean_spd_matrix = single_matrix, 
                   distances_to_mean = setNames(0, single_matrix_name),
                   num_valid_subjects = 1,
                   original_num_subjects = original_num_subjects,
                   mean_dispersion = 0, 
                   median_dispersion = 0))
  }
  
  # Proceed with computing the mean and dispersion for the valid subset of matrices
  # The frechet_mean_spd function is from R/riemannian_geometry.R
  
  mean_matrix_metric <- if (use_geometric_median && requireNamespace("shapes", quietly = TRUE)) {
    "airm"
  } else {
    "logeuclidean"
  }
  
  if (use_geometric_median && !requireNamespace("shapes", quietly = TRUE) && verbose) {
    message_stage("`shapes` package not available for AIRM Fréchet mean. Using iterative Log-Euclidean mean.", interactive_only = TRUE)
  } else if (use_geometric_median && verbose && interactive()) {
    message_stage(sprintf("Using %s Fréchet mean (via shapes::mediancov for AIRM, or iterative for LogEuclidean if shapes not used/available).",toupper(mean_matrix_metric)), interactive_only=TRUE)
  } else if (verbose && interactive()){
     message_stage("Using iterative Log-Euclidean Fréchet mean.", interactive_only=TRUE)
  }

  # frechet_mean_spd expects a simple list of matrices, not named.
  # Keep valid_spd_matrices_with_names for later to name the distances_to_mean vector.
  current_mean_matrix <- frechet_mean_spd(S_list = unname(valid_spd_matrices_with_names), 
                                          metric = mean_matrix_metric, 
                                          regularize_epsilon = spd_regularize_epsilon,
                                          verbose = verbose, 
                                          ...) # Pass relevant iter solver params if in ...
  
  if (is.null(current_mean_matrix)) {
      warning("Failed to compute mean SPD matrix. Cannot calculate dispersion.")
      return(list(mean_spd_matrix = NULL, 
                  distances_to_mean = setNames(rep(NA_real_, num_valid_subjects), names(valid_spd_matrices_with_names)),
                  num_valid_subjects = num_valid_subjects,
                  original_num_subjects = original_num_subjects,
                  mean_dispersion = NA_real_, 
                  median_dispersion = NA_real_))
  }

  # Distances from each valid SPD matrix to the computed mean matrix.
  # Use LogEuclidean distance by default, matching riemannian_distance_spd.
  # If AIRM mean was computed, there's a slight metric mismatch here, but often accepted.
  # Could add a 'distance_metric' param if strict consistency is needed.
  distances_to_mean_vals <- sapply(valid_spd_matrices_with_names, function(S_i) {
    tryCatch(
        riemannian_distance_spd(S_i, current_mean_matrix, regularize_epsilon = spd_regularize_epsilon, ...),
        error = function(e) {
            warning(sprintf("Distance calc failed for one subject to mean: %s", e$message)); NA_real_
        })
  })
  # names(distances_to_mean_vals) is already set because valid_spd_matrices_with_names is named.
  
  return(list(
    mean_spd_matrix = current_mean_matrix,
    distances_to_mean = distances_to_mean_vals, # This is now a named vector
    num_valid_subjects = num_valid_subjects,
    original_num_subjects = original_num_subjects,
    mean_dispersion = mean(distances_to_mean_vals, na.rm = TRUE),
    median_dispersion = stats::median(distances_to_mean_vals, na.rm = TRUE)
  ))
}
```

---

**4. S3 Method for Tangent Space Embedding (Fulfilling RGEOM-005)**

This S3 method provides functionality to project SPD matrix representations of subjects into a common tangent space, typically anchored at their Fréchet mean. This allows for applying standard Euclidean multivariate analyses to these Riemannian data.

```R
# R/riemannian_methods_hatsa.R (or R/tangent_space_embedding.R)

#' Project Subject SPD Representations to a Common Tangent Space
#'
#' Computes the Fréchet mean of specified SPD representations for a group of
#' subjects and then projects each subject's SPD matrix to the tangent space
#' anchored at this mean.
#'
#' @param object A `hatsa_projector` or `task_hatsa_projector` object.
#' @param type Character string indicating the type of SPD representation to use
#'   (e.g., `"cov_coeffs"`). Passed to `get_spd_representations`.
#' @param tangent_metric Character string, `"logeuclidean"` (default) or `"airm"`,
#'   specifying which metric's log map to use for projection. This also influences
#'   the metric for Fréchet mean computation if not overridden by `mean_options`.
#' @param subject_data_list Optional. Passed to `get_spd_representations`.
#' @param k_conn_params Optional. Passed to `get_spd_representations`.
#' @param spd_regularize_epsilon Epsilon for regularizing SPD matrices.
#' @param mean_options A list of options for `frechet_mean_spd` (e.g., `max_iter`, `tol`).
#'   The `metric` argument within this list will override `tangent_metric` for the mean computation step if provided.
#' @param verbose Logical, controls verbosity.
#' @param ... Additional arguments passed to `get_spd_representations`.
#' @return A list containing:
#'   - `tangent_vectors`: A list of symmetric matrices (tangent vectors), one per valid subject.
#'                      Names correspond to subject original indices.
#'   - `mean_spd_matrix`: The Fréchet mean SPD matrix used as the tangent point.
#'   - `metric_used`: The metric used for the log map (`"logeuclidean"` or `"airm"`).
#'   - `num_valid_subjects`: Number of subjects for whom tangent vectors were computed.
#'   - `original_num_subjects`: Total number of subjects in the input object.
#' @export
get_tangent_space_coords <- function(object, ...) {
  UseMethod("get_tangent_space_coords")
}

#' @rdname get_tangent_space_coords
#' @export
get_tangent_space_coords.hatsa_projector <- function(object,
                                                     type = "cov_coeffs",
                                                     tangent_metric = c("logeuclidean", "airm"),
                                                     subject_data_list = NULL,
                                                     k_conn_params = NULL,
                                                     spd_regularize_epsilon = 1e-6,
                                                     mean_options = list(),
                                                     verbose = FALSE,
                                                     ...) {
  tangent_metric <- match.arg(tangent_metric)

  # 1. Get SPD matrices for all subjects
  spd_matrices_list <- get_spd_representations(object,
                                               type = type,
                                               subject_idx = NULL, # All subjects
                                               regularize_epsilon = spd_regularize_epsilon,
                                               subject_data_list_for_fc = subject_data_list,
                                               k_conn_params_for_fc = k_conn_params,
                                               ...)
  
  valid_spd_matrices <- Filter(function(m) !is.null(m) && is.matrix(m) && !all(is.na(m)) && all(dim(m) > 0), spd_matrices_list)
  
  num_valid_subjects <- length(valid_spd_matrices)
  original_num_subjects <- object$parameters$N_subjects

  if (num_valid_subjects < 1) {
    if(verbose) message_stage("No valid SPD matrices to compute Fréchet mean or tangent vectors.", interactive_only=TRUE)
    return(list(tangent_vectors = list(), mean_spd_matrix = NULL, metric_used = tangent_metric, 
                num_valid_subjects = 0, original_num_subjects = original_num_subjects))
  }

  # 2. Compute Fréchet mean
  # Determine metric for Fréchet mean: use mean_options$metric if provided, else tangent_metric
  frechet_mean_metric <- mean_options$metric %||% tangent_metric 
  # (Using %||% as a placeholder for: if (!is.null(mean_options$metric)) mean_options$metric else tangent_metric)
  if (is.null(frechet_mean_metric)) frechet_mean_metric <- tangent_metric # Ensure it's set

  # Prepare arguments for frechet_mean_spd, allowing overrides from mean_options
  frechet_args <- list(
    S_list = unname(valid_spd_matrices),
    metric = frechet_mean_metric,
    regularize_epsilon = spd_regularize_epsilon,
    verbose = verbose
  )
  # Merge with user-provided mean_options, giving precedence to mean_options
  # for common parameters like max_iter, tol.
  if (length(mean_options) > 0) {
      for(opt_name in names(mean_options)) {
          frechet_args[[opt_name]] <- mean_options[[opt_name]]
      }
  }
  
  mean_spd <- do.call(frechet_mean_spd, frechet_args)

  if (is.null(mean_spd)) {
    if(verbose) message_stage("Failed to compute Fréchet mean. Cannot compute tangent space coordinates.", interactive_only=TRUE)
    return(list(tangent_vectors = list(), mean_spd_matrix = NULL, metric_used = tangent_metric, 
                num_valid_subjects = num_valid_subjects, original_num_subjects = original_num_subjects))
  }

  # 3. Map each SPD matrix to the tangent space at the mean
  logmap_function <- if (tangent_metric == "logeuclidean") {
    logmap_spd_logeuclidean
  } else { # airm
    logmap_spd_airm
  }

  tangent_vectors_list <- vector("list", num_valid_subjects)
  names(tangent_vectors_list) <- names(valid_spd_matrices) # Preserve subject names/indices

  for (i in seq_along(valid_spd_matrices)) {
    subj_name <- names(valid_spd_matrices)[i]
    S_i <- valid_spd_matrices[[i]]
    tv_i <- tryCatch(
      logmap_function(mean_spd, S_i, regularize_epsilon = spd_regularize_epsilon),
      error = function(e) {
        warning(sprintf("Logmap failed for subject %s: %s", subj_name, e$message))
        NULL
      }
    )
    tangent_vectors_list[[subj_name]] <- tv_i
  }
  
  # Filter out any NULLs from failed logmaps (though ideally logmap should be robust)
  successful_tangent_vectors <- Filter(Negate(is.null), tangent_vectors_list)
  num_successful_projections <- length(successful_tangent_vectors)

  if (num_successful_projections < num_valid_subjects && verbose){
      message_stage(sprintf("Successfully projected %d out of %d valid SPD matrices to tangent space.", 
                          num_successful_projections, num_valid_subjects), interactive_only=TRUE)
  }

  return(list(
    tangent_vectors = successful_tangent_vectors,
    mean_spd_matrix = mean_spd,
    metric_used = tangent_metric,
    num_valid_subjects = num_successful_projections, # Report how many were actually projected
    original_num_subjects = original_num_subjects
  ))
}
```

This outlines the S3 generic and a method for `hatsa_projector`. A similar method would be needed for `task_hatsa_projector` or it could leverage this one via `NextMethod()` if appropriate.

**Potential Output Formats for Tangent Vectors:**

While the above returns a list of symmetric matrices (tangent vectors), for direct use in many R functions (e.g., `prcomp`, `lm`), these tangent vectors often need to be vectorized. A common approach for a symmetric `k x k` matrix is to take its `k*(k+1)/2` unique elements (e.g., upper or lower triangle including diagonal).

An additional helper function or an option within `get_tangent_space_coords` could be added to return a single data matrix where each row is a subject and columns are the vectorized tangent space coordinates.

Example (conceptual helper for vectorization):
```R
.vectorize_symmetric_matrix <- function(symm_matrix, method = "upper_tri_diag") {
  if (!isSymmetric.matrix(symm_matrix)) stop("Matrix must be symmetric.")
  if (method == "upper_tri_diag") {
    return(symm_matrix[upper.tri(symm_matrix, diag = TRUE)])
  }
  # Add other methods if needed (e.g., "lower_tri_diag", "flatten_half_scaled_diag")
  stop("Unsupported vectorization method.")
}
```
This vectorization aspect can be considered an enhancement or a separate utility built upon the primary tangent space projection.

For instance, the `RiemBase` package (Kisung You, CRAN) provides a collection of tools for statistics on manifolds, including implementations for Fréchet means and operations on various manifolds like Grassmann and Stiefel. While HATSA is developing its own tailored geometric functions, `RiemBase` and similar packages can serve as valuable references or potential sources for more advanced/optimized routines if needed in the future.

---

# 6. Future Considerations for Geo-HATSA and Riemannian Toolkit

Building upon the implemented geometric tools and the Geo-HATSA framework, several avenues for future development and refinement exist:

*   **Covariance Normalization in Geometric GPA:**
    *   **Rationale:** The audit feedback (Gap 2e) suggested that geometric GPA often benefits from row-wise `sqrt(degree)` weighting of the anchor matrices (`A_i = U_original_list[[i]][unique_anchor_indices, , drop = FALSE]`) before the Procrustes alignment steps in `perform_geometric_gpa_refinement`.
    *   **Action:** Investigate implementing an optional normalization step, such as `A_i_norm = D_i^(1/2) A_i` (where `D_i` could be derived from row norms or similar weighting scheme), to assess its impact on the stability and quality of the geometric alignment.

*   **Full Riemannian Optimization for SO(k) GPA:**
    *   **Rationale:** The current `perform_geometric_gpa_refinement` (based on the SVD approach) provides a robust and efficient Euclidean approximation for updates on SO(k). For more advanced geometric modeling, especially in cases of high manifold curvature, a full Riemannian optimization approach could be beneficial.
    *   **Action:** Explore replacing the SVD-based rotation update with methods like Riemannian gradient descent or trust-region methods directly on the SO(k) manifold. This would involve minimizing a cost function based on geodesic distances. This could leverage external libraries (e.g., `geomstats`) or require custom implementations of SO(k) manifold operations (log map, exp map, parallel transport).

*   **Comprehensive Benchmarking of Geo-HATSA:**
    *   **Rationale:** To establish the benefits and trade-offs of Geo-HATSA compared to the standard Euclidean HATSA.
    *   **Action:** Conduct thorough benchmarking:
        *   **Synthetic Data:** Use datasets with known ground truth for rotations and templates to assess accuracy and parameter recovery.
        *   **Real-World Data:** Apply to datasets like HCP-15 (as suggested in the audit) to compare alignment consistency (e.g., Inter-Subject Correlation - ICC of aligned time series or features), robustness to noise, and computational performance.
        *   **Comparison Points:** Evaluate against `run_hatsa_core` and potentially other alignment methods. Include comparisons with and without covariance normalization if implemented.

*   **Vectorization of Tangent Space Coordinates:**
    *   **Rationale:** The `get_tangent_space_coords` function returns a list of symmetric matrices (tangent vectors). For direct use in many standard R multivariate analysis functions (e.g., `prcomp`, `lm`, clustering algorithms), these tangent vectors often need to be vectorized.
    *   **Action:** Implement a helper function (e.g., `.vectorize_symmetric_matrix`) or add an option to `get_tangent_space_coords` to convert the list of tangent vector matrices into a single data matrix (subjects x vectorized features). Common vectorization methods include taking the upper or lower triangle elements (including the diagonal). This would enhance the usability of tangent space embeddings for downstream statistical analyses. (This was previously noted in Section 4 of the plan and is reiterated here for completeness in future work).

```R
# R/riemannian_methods_hatsa.R (or R/tangent_space_embedding.R)

#' Project Subject SPD Representations to a Common Tangent Space
#'
#' Computes the Fréchet mean of specified SPD representations for a group of
#' subjects and then projects each subject's SPD matrix to the tangent space
#' anchored at this mean.
#'
#' @param object A `hatsa_projector` or `task_hatsa_projector` object.
#' @param type Character string indicating the type of SPD representation to use
#'   (e.g., `"cov_coeffs"`). Passed to `get_spd_representations`.
#' @param tangent_metric Character string, `"logeuclidean"` (default) or `"airm"`,
#'   specifying which metric's log map to use for projection. This also influences
#'   the metric for Fréchet mean computation if not overridden by `mean_options`.
#' @param subject_data_list Optional. Passed to `get_spd_representations`.
#' @param k_conn_params Optional. Passed to `get_spd_representations`.
#' @param spd_regularize_epsilon Epsilon for regularizing SPD matrices.
#' @param mean_options A list of options for `frechet_mean_spd` (e.g., `max_iter`, `tol`).
#'   The `metric` argument within this list will override `tangent_metric` for the mean computation step if provided.
#' @param verbose Logical, controls verbosity.
#' @param ... Additional arguments passed to `get_spd_representations`.
#' @return A list containing:
#'   - `tangent_vectors`: A list of symmetric matrices (tangent vectors), one per valid subject.
#'                      Names correspond to subject original indices.
#'   - `mean_spd_matrix`: The Fréchet mean SPD matrix used as the tangent point.
#'   - `metric_used`: The metric used for the log map (`"logeuclidean"` or `"airm"`).
#'   - `num_valid_subjects`: Number of subjects for whom tangent vectors were computed.
#'   - `original_num_subjects`: Total number of subjects in the input object.
#' @export
get_tangent_space_coords <- function(object, ...) {
  UseMethod("get_tangent_space_coords")
}

#' @rdname get_tangent_space_coords
#' @export
get_tangent_space_coords.hatsa_projector <- function(object,
                                                     type = "cov_coeffs",
                                                     tangent_metric = c("logeuclidean", "airm"),
                                                     subject_data_list = NULL,
                                                     k_conn_params = NULL,
                                                     spd_regularize_epsilon = 1e-6,
                                                     mean_options = list(),
                                                     verbose = FALSE,
                                                     ...) {
  tangent_metric <- match.arg(tangent_metric)

  # 1. Get SPD matrices for all subjects
  spd_matrices_list <- get_spd_representations(object,
                                               type = type,
                                               subject_idx = NULL, # All subjects
                                               regularize_epsilon = spd_regularize_epsilon,
                                               subject_data_list_for_fc = subject_data_list,
                                               k_conn_params_for_fc = k_conn_params,
                                               ...)
  
  valid_spd_matrices <- Filter(function(m) !is.null(m) && is.matrix(m) && !all(is.na(m)) && all(dim(m) > 0), spd_matrices_list)
  
  num_valid_subjects <- length(valid_spd_matrices)
  original_num_subjects <- object$parameters$N_subjects

  if (num_valid_subjects < 1) {
    if(verbose) message_stage("No valid SPD matrices to compute Fréchet mean or tangent vectors.", interactive_only=TRUE)
    return(list(tangent_vectors = list(), mean_spd_matrix = NULL, metric_used = tangent_metric, 
                num_valid_subjects = 0, original_num_subjects = original_num_subjects))
  }

  # 2. Compute Fréchet mean
  # Determine metric for Fréchet mean: use mean_options$metric if provided, else tangent_metric
  frechet_mean_metric <- mean_options$metric %||% tangent_metric 
  # (Using %||% as a placeholder for: if (!is.null(mean_options$metric)) mean_options$metric else tangent_metric)
  if (is.null(frechet_mean_metric)) frechet_mean_metric <- tangent_metric # Ensure it's set

  # Prepare arguments for frechet_mean_spd, allowing overrides from mean_options
  frechet_args <- list(
    S_list = unname(valid_spd_matrices),
    metric = frechet_mean_metric,
    regularize_epsilon = spd_regularize_epsilon,
    verbose = verbose
  )
  # Merge with user-provided mean_options, giving precedence to mean_options
  # for common parameters like max_iter, tol.
  if (length(mean_options) > 0) {
      for(opt_name in names(mean_options)) {
          frechet_args[[opt_name]] <- mean_options[[opt_name]]
      }
  }
  
  mean_spd <- do.call(frechet_mean_spd, frechet_args)

  if (is.null(mean_spd)) {
    if(verbose) message_stage("Failed to compute Fréchet mean. Cannot compute tangent space coordinates.", interactive_only=TRUE)
    return(list(tangent_vectors = list(), mean_spd_matrix = NULL, metric_used = tangent_metric, 
                num_valid_subjects = num_valid_subjects, original_num_subjects = original_num_subjects))
  }

  # 3. Map each SPD matrix to the tangent space at the mean
  logmap_function <- if (tangent_metric == "logeuclidean") {
    logmap_spd_logeuclidean
  } else { # airm
    logmap_spd_airm
  }

  tangent_vectors_list <- vector("list", num_valid_subjects)
  names(tangent_vectors_list) <- names(valid_spd_matrices) # Preserve subject names/indices

  for (i in seq_along(valid_spd_matrices)) {
    subj_name <- names(valid_spd_matrices)[i]
    S_i <- valid_spd_matrices[[i]]
    tv_i <- tryCatch(
      logmap_function(mean_spd, S_i, regularize_epsilon = spd_regularize_epsilon),
      error = function(e) {
        warning(sprintf("Logmap failed for subject %s: %s", subj_name, e$message))
        NULL
      }
    )
    tangent_vectors_list[[subj_name]] <- tv_i
  }
  
  # Filter out any NULLs from failed logmaps (though ideally logmap should be robust)
  successful_tangent_vectors <- Filter(Negate(is.null), tangent_vectors_list)
  num_successful_projections <- length(successful_tangent_vectors)

  if (num_successful_projections < num_valid_subjects && verbose){
      message_stage(sprintf("Successfully projected %d out of %d valid SPD matrices to tangent space.", 
                          num_successful_projections, num_valid_subjects), interactive_only=TRUE)
  }

  return(list(
    tangent_vectors = successful_tangent_vectors,
    mean_spd_matrix = mean_spd,
    metric_used = tangent_metric,
    num_valid_subjects = num_successful_projections, # Report how many were actually projected
    original_num_subjects = original_num_subjects
  ))
}