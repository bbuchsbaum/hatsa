# S3 Methods for Riemannian Geometry computations on HATSA objects

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
#' @param k_conn_params Parameters for connectivity calculation if type involves FC,
#'   passed to `get_spd_representations`.
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
  # N_subjects derived from spd_matrices_list length or object parameters
  
  spd_matrices_list_all_subjects <- get_spd_representations(object,
                                                            type = type,
                                                            subject_idx = NULL, # Get for all subjects
                                                            regularize_epsilon = spd_regularize_epsilon,
                                                            subject_data_list_for_fc = subject_data_list,
                                                            k_conn_params_for_fc = k_conn_params,
                                                            ...)
  
  valid_indices_in_list <- which(!sapply(spd_matrices_list_all_subjects, is.null))
  
  # Use names from the list returned by get_spd_representations if available
  # These names should correspond to original subject indices as characters
  all_subject_names_in_input_list <- names(spd_matrices_list_all_subjects)
  
  # Determine the total number of subjects based on the input object
  N_total_subjects_in_object <- object$parameters$N_subjects
  # Generate character names for all subjects in the object for the final matrix
  all_original_subj_indices_char <- as.character(1:N_total_subjects_in_object)

  if (length(valid_indices_in_list) < length(spd_matrices_list_all_subjects)) {
    num_failed = length(spd_matrices_list_all_subjects) - length(valid_indices_in_list)
    # Ensure spd_matrices_list_all_subjects is not NULL or empty before trying to get its length
    total_subjects_in_list = if(!is.null(spd_matrices_list_all_subjects)) length(spd_matrices_list_all_subjects) else 0
    
    warning_msg <- sprintf("Could not obtain valid SPD matrices for %d out of %d subject(s) processed by get_spd_representations. Distances involving them will be NA.", 
                           num_failed, total_subjects_in_list)
    if (total_subjects_in_list != N_total_subjects_in_object) {
        warning_msg <- paste(warning_msg, sprintf("Note: get_spd_representations processed %d subjects, but object has %d total subjects.", total_subjects_in_list, N_total_subjects_in_object))
    }
    warning(warning_msg)
  }
  
  # Initialize the final distance matrix with NAs, sized for all subjects in the object
  final_dist_mat <- matrix(NA_real_, nrow = N_total_subjects_in_object, ncol = N_total_subjects_in_object)
  rownames(final_dist_mat) <- colnames(final_dist_mat) <- all_original_subj_indices_char
  diag(final_dist_mat) <- 0.0

  # If no valid matrices or only one, distance matrix computation for pairs is not possible
  if (length(valid_indices_in_list) < 2) {
    warning("Need at least 2 subjects with valid SPD matrices to compute a pairwise distance matrix. Returning matrix with NAs (and 0s on diagonal).")
    return(final_dist_mat) # Returns the initialized NA matrix with 0s on diagonal
  }

  spd_matrices_valid_subset <- spd_matrices_list_all_subjects[valid_indices_in_list]
  num_valid_subjects <- length(spd_matrices_valid_subset)
  valid_subj_original_indices_char <- names(spd_matrices_valid_subset) 

  # Create a temporary distance matrix for the valid subset
  dist_mat_subset <- matrix(NA_real_, nrow = num_valid_subjects, ncol = num_valid_subjects)
  rownames(dist_mat_subset) <- colnames(dist_mat_subset) <- valid_subj_original_indices_char
  diag(dist_mat_subset) <- 0.0
  
  for (i in 1:(num_valid_subjects - 1)) {
    idx1_name <- valid_subj_original_indices_char[i]
    S1 <- spd_matrices_valid_subset[[idx1_name]]
    
    for (j in (i + 1):num_valid_subjects) {
      idx2_name <- valid_subj_original_indices_char[j]
      S2 <- spd_matrices_valid_subset[[idx2_name]]
      
      current_dist <- NA_real_ # Default to NA
      if(!is.null(S1) && !is.null(S2)){ # Ensure matrices themselves are not null (double check after list filter)
         current_dist <- tryCatch(
            riemannian_distance_spd(S1, S2, regularize_epsilon = spd_regularize_epsilon, ...),
            error = function(e) {
              warning(sprintf("Riemannian distance computation failed between subject %s and %s: %s", idx1_name, idx2_name, e$message))
              NA_real_
            }
          )
      } else {
           warning(sprintf("Skipping distance between subject %s and %s due to one or both SPD matrices being NULL.", idx1_name, idx2_name))
      }
      dist_mat_subset[idx1_name, idx2_name] <- current_dist
      dist_mat_subset[idx2_name, idx1_name] <- current_dist # Symmetrize
    }
  }
  
  # Fill in the computed distances for valid subjects into the full-sized matrix
  # This uses character indexing which should align names correctly.
  final_dist_mat[valid_subj_original_indices_char, valid_subj_original_indices_char] <- dist_mat_subset
  
  return(final_dist_mat)
}


#' Compute Riemannian Dispersion of SPD Matrices on the Manifold
#'
#' Calculates the dispersion of a set of SPD matrices around their Fréchet mean
#' using the Riemannian (Log-Euclidean by default) distance.
#'
#' @param object A `hatsa_projector` or `task_hatsa_projector` object.
#' @param type Character string, either `"cov_coeffs"` (default) or others like
#'   `"fc_conn"`, `"fc_task"` (relying on RGEOM-003 `get_spd_representations`).
#'   Specifies which matrices to compute dispersion for.
#' @param subject_data_list Optional, may be needed by `get_spd_representations`.
#' @param use_geometric_median Logical, if TRUE, attempts to use AIRM Fréchet mean 
#'   (via `shapes::mediancov` if available, or iterative Log-Euclidean if not). 
#'   If FALSE (default), uses iterative Log-Euclidean Fréchet mean.
#' @param k_conn_params Parameters for connectivity if type involves FC.
#' @param spd_regularize_epsilon Epsilon for regularizing SPD matrices.
#' @param verbose Logical, controls verbose output during calculations.
#' @param ... Additional arguments passed to `frechet_mean_spd`, 
#'   `riemannian_distance_spd` or `get_spd_representations`.
#' @return A list containing: `mean_spd_matrix`, `distances_to_mean` (named vector), 
#'   `num_valid_subjects`, `original_num_subjects`, `mean_dispersion`, `median_dispersion`.
#' @export
#' @importFrom stats median
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
  
  spd_matrices_list_all_subjects <- get_spd_representations(object,
                                                            type = type,
                                                            subject_idx = NULL, 
                                                            regularize_epsilon = spd_regularize_epsilon,
                                                            subject_data_list_for_fc = subject_data_list,
                                                            k_conn_params_for_fc = k_conn_params,
                                                            ...) 
  
  valid_spd_matrices_with_names <- Filter(function(m) !is.null(m) && is.matrix(m) && !all(is.na(m)) && all(dim(m) > 0), spd_matrices_list_all_subjects)
  
  num_valid_subjects <- length(valid_spd_matrices_with_names)
  original_num_subjects <- object$parameters$N_subjects

  empty_return <- list(mean_spd_matrix = NULL, 
                       distances_to_mean = setNames(numeric(0), character(0)),
                       num_valid_subjects = num_valid_subjects, # Could be 0
                       original_num_subjects = original_num_subjects,
                       mean_dispersion = NA_real_, 
                       median_dispersion = NA_real_)

  if (num_valid_subjects == 0) {
    if(verbose) message_stage("No valid SPD matrices found to compute dispersion.", interactive_only=TRUE)
    # Update num_valid_subjects in the return list for consistency
    empty_return$num_valid_subjects <- 0
    return(empty_return)
  }
  
  if (num_valid_subjects == 1) {
      single_matrix <- valid_spd_matrices_with_names[[1]]
      single_matrix_name <- names(valid_spd_matrices_with_names)[1]
      if (is.null(single_matrix_name) || single_matrix_name == "") single_matrix_name <- "1" 

      warning_msg <- if (original_num_subjects > 1) {
          "Only one valid SPD matrix found out of multiple subjects. Dispersion is 0 for this single matrix."
      } else {
          "Only one subject provided. Dispersion is 0."
      }
      if(verbose) message_stage(warning_msg, interactive_only=TRUE)
      
       return(list(mean_spd_matrix = single_matrix, 
                   distances_to_mean = setNames(0, single_matrix_name),
                   num_valid_subjects = 1,
                   original_num_subjects = original_num_subjects,
                   mean_dispersion = 0, 
                   median_dispersion = 0))
  }
  
  mean_matrix_metric <- "logeuclidean" # Default
  if (use_geometric_median) {
      if (requireNamespace("shapes", quietly = TRUE)) {
          mean_matrix_metric <- "airm"
          if (verbose && interactive()) message_stage("Attempting AIRM Fréchet mean (via shapes::mediancov).", interactive_only=TRUE)
      } else {
          if (verbose && interactive()) message_stage("`shapes` package not available for AIRM Fréchet mean. Using iterative Log-Euclidean mean.", interactive_only=TRUE)
      }
  } else {
      if (verbose && interactive()) message_stage("Using iterative Log-Euclidean Fréchet mean.", interactive_only=TRUE)
  }

  current_mean_matrix <- frechet_mean_spd(S_list = unname(valid_spd_matrices_with_names), 
                                          metric = mean_matrix_metric, 
                                          regularize_epsilon = spd_regularize_epsilon,
                                          verbose = verbose, 
                                          ...) 
  
  if (is.null(current_mean_matrix)) {
      if(verbose) message_stage("Failed to compute mean SPD matrix. Cannot calculate dispersion.", interactive_only=TRUE)
      # Update distances_to_mean to be correctly named NA vector
      empty_return$distances_to_mean <- setNames(rep(NA_real_, num_valid_subjects), names(valid_spd_matrices_with_names))
      return(empty_return)
  }

  # Determine distance metric for individual distances to mean. Default to LogEuclidean.
  # This might differ from mean_matrix_metric if shapes was used for AIRM mean.
  distance_metric_to_mean <- "logeuclidean" 
  
  distances_to_mean_vals <- sapply(valid_spd_matrices_with_names, function(S_i) {
    dist_val <- NA_real_
    if (!is.null(S_i) && !is.null(current_mean_matrix)){
         dist_val <- tryCatch({
            if (distance_metric_to_mean == "airm") {
                airm_distance(S_i, current_mean_matrix, regularize_epsilon = spd_regularize_epsilon, ...)
            } else { # Default to LogEuclidean
                riemannian_distance_spd(S_i, current_mean_matrix, regularize_epsilon = spd_regularize_epsilon, ...)
            }
         }, error = function(e) {
            warning(sprintf("Distance calculation failed for subject '%s' to mean: %s", 
                            # Attempt to get name of S_i if list was named, else use placeholder
                            names(S_i) %||% "unknown", # Using a hypothetical %||% for brevity
                            e$message)); 
            NA_real_
        })
    }
    dist_val
  })
  # sapply should preserve names from valid_spd_matrices_with_names
  
  return(list(
    mean_spd_matrix = current_mean_matrix,
    distances_to_mean = distances_to_mean_vals, 
    num_valid_subjects = num_valid_subjects,
    original_num_subjects = original_num_subjects,
    mean_dispersion = mean(distances_to_mean_vals, na.rm = TRUE),
    median_dispersion = stats::median(distances_to_mean_vals, na.rm = TRUE)
  ))
}

# Placeholder for task_hatsa_projector methods if they are to be different
# If they are identical, we can rely on class inheritance or explicitly define them
# to call the hatsa_projector version.

# For now, assume task_hatsa_projector might have its own nuances or might just call the parent.
# Example:
# riemannian_distance_matrix_spd.task_hatsa_projector <- function(object, ...) {
#   # Could add task-specific logic or simply call the parent method:
#   # NextMethod() 
#   # OR
#   # riemannian_distance_matrix_spd.hatsa_projector(object, ...)
#   # For now, let users define it if specific behavior is needed, or rely on NextMethod if appropriate.
#   message("Using default riemannian_distance_matrix_spd for task_hatsa_projector. Define a specific method if needed.")
#   NextMethod()
# }

# riemannian_dispersion_spd.task_hatsa_projector <- function(object, ...) {
#   message("Using default riemannian_dispersion_spd for task_hatsa_projector. Define a specific method if needed.")
#   NextMethod()
# } 

#' Project Subject SPD Representations to a Common Tangent Space
#'
#' Computes the Fréchet mean of specified SPD representations for a group of
#' subjects and then projects each subject's SPD matrix to the tangent space
#' anchored at this mean. This allows for applying standard Euclidean
#' multivariate analyses to these Riemannian data.
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
#'   The `metric` argument within this list will override `tangent_metric` for the mean
#'   computation step if provided.
#' @param verbose Logical, controls verbosity.
#' @param ... Additional arguments passed to `get_spd_representations` or `frechet_mean_spd`.
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
    if(verbose && interactive()) message_stage("No valid SPD matrices to compute Fréchet mean or tangent vectors.", interactive_only=TRUE)
    return(list(tangent_vectors = list(), mean_spd_matrix = NULL, metric_used = tangent_metric, 
                num_valid_subjects = 0, original_num_subjects = original_num_subjects))
  }

  # 2. Compute Fréchet mean
  frechet_mean_metric <- if (!is.null(mean_options$metric)) mean_options$metric else tangent_metric

  frechet_args <- list(
    S_list = unname(valid_spd_matrices),
    metric = frechet_mean_metric,
    regularize_epsilon = spd_regularize_epsilon,
    verbose = verbose
  )
  # Merge with user-provided mean_options, giving precedence to mean_options
  # for common parameters like max_iter, tol, etc. and also pass further ... args from main call
  custom_frechet_params <- utils::modifyList(mean_options, list(...))
  # Remove params already explicitly set or those not for frechet_mean_spd
  custom_frechet_params$S_list <- NULL 
  custom_frechet_params$metric <- NULL
  custom_frechet_params$regularize_epsilon <- NULL
  custom_frechet_params$verbose <- NULL
  
  final_frechet_args <- utils::modifyList(frechet_args, custom_frechet_params)
  
  mean_spd <- do.call(frechet_mean_spd, final_frechet_args)

  if (is.null(mean_spd)) {
    if(verbose && interactive()) message_stage("Failed to compute Fréchet mean. Cannot compute tangent space coordinates.", interactive_only=TRUE)
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
  names(tangent_vectors_list) <- names(valid_spd_matrices) 

  for (i in seq_along(valid_spd_matrices)) {
    subj_name_for_warning <- names(valid_spd_matrices)[i]
    if(is.null(subj_name_for_warning) || subj_name_for_warning == "") subj_name_for_warning <- as.character(i)

    S_i <- valid_spd_matrices[[i]]
    tv_i <- tryCatch(
      logmap_function(mean_spd, S_i, regularize_epsilon = spd_regularize_epsilon),
      error = function(e) {
        warning(sprintf("Logmap failed for subject %s: %s", subj_name_for_warning, e$message))
        NULL
      }
    )
    # Use the original name from valid_spd_matrices for assignment
    tangent_vectors_list[[names(valid_spd_matrices)[i]]] <- tv_i 
  }
  
  successful_tangent_vectors <- Filter(Negate(is.null), tangent_vectors_list)
  num_successful_projections <- length(successful_tangent_vectors)

  if (num_successful_projections < num_valid_subjects && verbose && interactive()){
      message_stage(sprintf("Successfully projected %d out of %d valid SPD matrices to tangent space.", 
                          num_successful_projections, num_valid_subjects), interactive_only=TRUE)
  }

  return(list(
    tangent_vectors = successful_tangent_vectors,
    mean_spd_matrix = mean_spd,
    metric_used = tangent_metric,
    num_valid_subjects = num_successful_projections, 
    original_num_subjects = original_num_subjects
  ))
}

# Potential helper for vectorizing symmetric matrices (conceptual)
# .vectorize_symmetric_matrix <- function(symm_matrix, method = "upper_tri_diag") {
#   if (!isSymmetric.matrix(symm_matrix)) stop("Matrix must be symmetric.")
#   if (method == "upper_tri_diag") {
#     return(symm_matrix[upper.tri(symm_matrix, diag = TRUE)])
#   }
#   stop("Unsupported vectorization method.")
# } 