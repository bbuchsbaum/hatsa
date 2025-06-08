#' Manifold-Regularized Anchor Selection (MRA-Select)
#'
#' Selects a set of anchor parcels by optimizing a criterion that balances
#' the stability of the mean anchor representation (via condition number) and
#' the dispersion of subject-specific covariance matrices derived from these anchors
#' on the SPD manifold.
#'
#' @param U_original_list_pilot A list of original (unaligned) sketch matrices
#'   (V_p x k) for a set of pilot subjects. These are used to evaluate candidate anchors.
#' @param k_spectral_rank Integer, the spectral rank k used to generate `U_original_list_pilot`.
#' @param m_target Integer, the desired number of anchors to select.
#' @param total_parcels Integer, the total number of parcels (V_p) available for selection.
#'   If NULL, inferred from `U_original_list_pilot`.
#' @param max_kappa Numeric, the maximum allowable condition number for the mean anchor matrix.
#'   Candidates leading to a kappa above this are penalized or excluded. Default: 100.
#' @param weight_inv_kappa Numeric, weight for the inverse condition number term in the score.
#'   Higher values prioritize lower condition numbers (more stable mean). Default: 1.0.
#' @param weight_dispersion Numeric, weight for the negative dispersion term in the score.
#'   Higher values prioritize lower dispersion. Default: 1.0.
#' @param initial_selection Optional. A vector of pre-selected anchor indices (1-based).
#'   The algorithm will try to add `m_target - length(initial_selection)` new anchors.
#' @param candidate_pool Optional. A vector of parcel indices (1-based) from which to select
#'   new anchors. If NULL, all non-initially-selected parcels are candidates.
#' @param parcel_quality_info Optional. Data frame or list providing parcel quality or
#'   network information. Currently unused in the core MRA logic but reserved for future
#'   extensions (e.g., pre-filtering candidates).
#' @param riemannian_dispersion_options A list of options to pass to
#'   `hatsa::riemannian_dispersion_spd` (e.g., `spd_metric`, `use_geometric_median`).
#' @param min_anchors_for_metrics Integer. Minimum number of selected anchors required before
#'   kappa and dispersion metrics are considered stable enough to compute. Default is `k_spectral_rank`.
#' @param verbose Logical. If TRUE, print progress messages. Default: TRUE.
#'
#' @return A vector of selected anchor indices (1-based), sorted.
#' @export
#' @importFrom stats cov sd svd
#'
#' @examples
#' # Conceptual example (requires U_original_list_pilot from actual or synthetic data)
#' # N_pilot <- 5
#' # V_p_total <- 50
#' # k_rank <- 10
#' # U_pilot_list <- replicate(N_pilot,
#' #                           matrix(rnorm(V_p_total * k_rank), V_p_total, k_rank),
#' #                           simplify = FALSE)
#' # selected_anchors_mra <- select_anchors_mra(
#' #   U_original_list_pilot = U_pilot_list,
#' #   k_spectral_rank = k_rank,
#' #   m_target = 15,
#' #   total_parcels = V_p_total,
#' #   verbose = TRUE
#' # )
#' # print(selected_anchors_mra)
select_anchors_mra <- function(U_original_list_pilot,
                               k_spectral_rank,
                               m_target,
                               total_parcels = NULL,
                               max_kappa = 100,
                               weight_inv_kappa = 1.0,
                               weight_dispersion = 1.0,
                               initial_selection = integer(0),
                               candidate_pool = NULL,
                               parcel_quality_info = NULL, # Reserved for future use
                               riemannian_dispersion_options = list(),
                               min_anchors_for_metrics = NULL,
                               verbose = TRUE) {

  # --- Input Validation and Setup ---
  if (!is.list(U_original_list_pilot) || length(U_original_list_pilot) == 0) {
    stop("U_original_list_pilot must be a non-empty list.")
  }
  if (is.null(total_parcels)) {
    if (!is.null(U_original_list_pilot[[1]]) && is.matrix(U_original_list_pilot[[1]])) {
      total_parcels <- nrow(U_original_list_pilot[[1]])
    } else {
      stop("total_parcels not provided and cannot be inferred from U_original_list_pilot.")
    }
  }
  if (m_target <= 0 || m_target > total_parcels) {
    stop("m_target must be positive and not exceed total_parcels.")
  }
  if (is.null(min_anchors_for_metrics)) {
    min_anchors_for_metrics <- k_spectral_rank
  }
  if (min_anchors_for_metrics < 1) min_anchors_for_metrics <- 1
  
  # Validate initial_selection and candidate_pool
  initial_selection <- unique(as.integer(initial_selection))
  if (any(initial_selection < 1 | initial_selection > total_parcels)) {
    stop("initial_selection contains invalid parcel indices.")
  }
  if (length(initial_selection) >= m_target) {
    if (verbose) message("Initial selection already meets or exceeds m_target. Returning initial selection.")
    return(sort(initial_selection))
  }

  if (is.null(candidate_pool)) {
    candidate_pool <- setdiff(1:total_parcels, initial_selection)
  } else {
    candidate_pool <- unique(as.integer(candidate_pool))
    if (any(candidate_pool < 1 | candidate_pool > total_parcels)) {
      stop("candidate_pool contains invalid parcel indices.")
    }
    candidate_pool <- setdiff(candidate_pool, initial_selection)
  }
  if (length(candidate_pool) == 0 && length(initial_selection) < m_target) {
     warning("No candidate anchors available to reach m_target. Returning current selection.")
     return(sort(initial_selection))
  }

  selected_anchors <- initial_selection
  num_to_select_additionally <- m_target - length(selected_anchors)

  if (verbose) {
    message_stage(sprintf("Starting MRA-Select: Target %d anchors (%d pre-selected, %d to find).", 
                          m_target, length(initial_selection), num_to_select_additionally), interactive_only = TRUE)
  }

  # --- Helper: Calculate Kappa ---
  .calculate_kappa <- function(matrix_in) {
    if (is.null(matrix_in) || !is.matrix(matrix_in) || min(dim(matrix_in)) == 0) return(Inf)
    # Kappa is typically for m x k where m >= k.
    # If m < k, it will be rank deficient, min singular value often 0 (or near 0 due to numerics).
    # svd() handles non-square matrices.
    s_vals <- tryCatch(svd(matrix_in, nu=0, nv=0)$d, error = function(e) NULL)
    if (is.null(s_vals) || length(s_vals) == 0) return(Inf)
    min_s_val <- min(s_vals[s_vals > .Machine$double.eps]) # Smallest positive singular value
    if (length(min_s_val) == 0 || min_s_val < .Machine$double.eps) return(Inf)
    return(max(s_vals) / min_s_val)
  }
  
  # --- Helper: Calculate Metrics for a Tentative Set ---
  .evaluate_anchor_set <- function(current_anchor_indices) {
    num_current_anchors <- length(current_anchor_indices)
    if (num_current_anchors < 1) return(list(kappa = Inf, dispersion = Inf, score = -Inf))

    # Cache anchor rows for each subject once
    anchor_rows_list <- lapply(U_original_list_pilot, function(U_subj) {
      if (is.null(U_subj) || nrow(U_subj) < max(current_anchor_indices)) return(NULL)
      U_subj[current_anchor_indices, , drop = FALSE]
    })

    # Metric 1: Kappa of the Euclidean mean of anchor rows from U_original_list_pilot
    kappa_val <- Inf
    if (num_current_anchors >= min_anchors_for_metrics || num_current_anchors >= k_spectral_rank) {
        valid_anchor_matrices <- Filter(Negate(is.null), anchor_rows_list)
        if (length(valid_anchor_matrices) > 0) {
            # Check all matrices have same dim: num_current_anchors x k_spectral_rank
            if (!all(sapply(valid_anchor_matrices, function(m) all(dim(m) == c(num_current_anchors, k_spectral_rank))))) {
                warning("Inconsistent dimensions in pilot anchor matrices for kappa calculation.")
            } else {
                mean_anchor_matrix <- Reduce("+", valid_anchor_matrices) / length(valid_anchor_matrices)
                kappa_val <- .calculate_kappa(mean_anchor_matrix)
            }
        } else {
            warning("No valid pilot anchor matrices for kappa calculation.")
        }
    }

    # Metric 2: Riemannian Dispersion of covariance matrices of anchor rows
    dispersion_val <- Inf
    if (num_current_anchors >= min_anchors_for_metrics || num_current_anchors >= k_spectral_rank) { # Need enough rows for cov
        spd_list_pilot <- lapply(anchor_rows_list, function(anchor_rows_subj) {
            if (is.null(anchor_rows_subj)) return(NULL)
            if (nrow(anchor_rows_subj) < 2 || nrow(anchor_rows_subj) < k_spectral_rank) return(NULL) # cov needs multiple observations
            tryCatch(stats::cov(anchor_rows_subj), error = function(e) NULL)
        })
        valid_spd_list <- Filter(Negate(is.null), spd_list_pilot)
        
        if (length(valid_spd_list) > 1) {
            # Ensure all SPD matrices are k x k
            valid_spd_list <- Filter(function(S) is.matrix(S) && all(dim(S) == c(k_spectral_rank, k_spectral_rank)), valid_spd_list)
            if (length(valid_spd_list) > 1) {
                disp_opts <- riemannian_dispersion_options
                # Construct a temporary projector-like object or pass list directly if API allows
                # For now, assuming riemannian_dispersion_spd can take a raw list of SPDs if object=NULL
                # This requires modification to riemannian_dispersion_spd or a new helper.
                # HACK: Create minimal list structure that riemannian_dispersion_spd might accept if it looks for object$N_subjects
                # This is not ideal. Better: riemannian_dispersion_spd should have a method for lists.
                # For now, let's assume we pass to a conceptual helper: calculate_dispersion_for_spd_list(valid_spd_list, disp_opts)
                # For a quick implementation, let's assume riemannian_dispersion_spd can be adapted or we use a direct calc:
                
                # Placeholder for direct calculation if riemannian_dispersion_spd cannot be easily used:
                # 1. Compute Frechet mean of valid_spd_list (e.g., using hatsa::frechet_mean_spd)
                # 2. Compute sum of squared distances to mean.
                # This is simplified here. Proper use of riemannian_dispersion_spd is preferred.
                
                # Simulate calling riemannian_dispersion_spd (needs a method for raw lists or dummy object)
                # This is a temporary workaround. The actual call would depend on how `riemannian_dispersion_spd` can handle raw SPD lists.
                # We will assume for now it fails gracefully if it cannot. The ticket implies using it.
                dummy_projector_for_disp <- list(
                    parameters = list(N_subjects = length(valid_spd_list)),
                    # get_spd_representations would be mocked to return valid_spd_list
                    # This part needs careful handling based on actual API of riemannian_dispersion_spd
                    .get_spd_representations_output = valid_spd_list # Internal convention for this example
                )
                class(dummy_projector_for_disp) <- "hatsa_projector" # To satisfy S3 dispatch if needed

                current_disp_opts <- riemannian_dispersion_options
                current_disp_opts$object <- NULL # Indicate we are providing spd_matrices_list directly
                current_disp_opts$spd_matrices_list <- valid_spd_list
                current_disp_opts$verbose <- FALSE # Suppress verbose from internal call usually

                # Need to ensure riemannian_dispersion_spd can accept spd_matrices_list
                # if object is NULL. This is a new requirement for that function.
                # For now, let's assume a direct calculation for simplicity of MRA-select draft:
                if (length(valid_spd_list) > 1) {
                    mean_spd <- tryCatch(hatsa::frechet_mean_spd(valid_spd_list, 
                                                                metric = current_disp_opts$spd_metric %||% "logeuclidean",
                                                                tol = current_disp_opts$frechet_tol %||% 1e-7,
                                                                max_iter = current_disp_opts$frechet_max_iter %||% 50,
                                                                verbose = FALSE), error = function(e) NULL)
                    if (!is.null(mean_spd)){
                        distances_sq <- sapply(valid_spd_list, function(S_i) {
                            dist_val <- tryCatch(hatsa::riemannian_distance_spd(S_i, mean_spd, 
                                                                             metric=current_disp_opts$spd_metric %||% "logeuclidean", 
                                                                             epsilon = current_disp_opts$spd_regularize_epsilon %||% 1e-6)^2,
                                                 error = function(e) NA)
                            return(dist_val)
                        })
                        dispersion_val <- mean(na.omit(distances_sq)) # Mean squared distance
                        if(is.nan(dispersion_val)) dispersion_val <- Inf
                    } else { dispersion_val <- Inf }
                } else { dispersion_val <- Inf }
            } else { dispersion_val <- Inf }
        } else { dispersion_val <- Inf }
    }

    # Combine metrics into a score (maximize this score)
    # Score = w1 * (1/kappa) - w2 * dispersion
    # If kappa is Inf or > max_kappa, 1/kappa term becomes 0 or score penalized
    inv_kappa_term <- 0
    if (is.finite(kappa_val) && kappa_val > .Machine$double.eps) {
        if (kappa_val <= max_kappa) {
            inv_kappa_term <- 1 / kappa_val
        } else {
            inv_kappa_term <- 1 / max_kappa # Penalize but don't make it zero if slightly over
                                         # or simply make score -Inf: if (kappa_val > max_kappa) score <- -Inf
        }
    }
    current_score <- weight_inv_kappa * inv_kappa_term - weight_dispersion * ifelse(is.finite(dispersion_val), dispersion_val, Inf)
    if (kappa_val > max_kappa) current_score <- -Inf # Hard constraint

    return(list(kappa = kappa_val, dispersion = dispersion_val, score = current_score))
  }

  # --- Greedy Selection Loop ---
  for (i in 1:num_to_select_additionally) {
    best_current_iteration_score <- -Inf
    best_candidate_this_iteration <- NULL

    if (length(candidate_pool) == 0) {
      if (verbose) message("No more candidates to select from.")
      break
    }

    if (verbose) {
        message(sprintf("  MRA Iteration %d/%d: Evaluating %d candidates to add to current %d anchors...", 
                        i, num_to_select_additionally, length(candidate_pool), length(selected_anchors)))
    }
    
    scores_list <- lapply(candidate_pool, function(candidate_anchor_idx) {
      tentative_anchors <- sort(unique(c(selected_anchors, candidate_anchor_idx)))
      metrics <- .evaluate_anchor_set(tentative_anchors)
      data.frame(candidate_idx = candidate_anchor_idx,
                 score = metrics$score,
                 kappa = metrics$kappa,
                 dispersion = metrics$dispersion)
    })

    iter_scores <- do.call(rbind, scores_list)
    if (nrow(iter_scores) > 0) {
      best_idx <- which.max(iter_scores$score)
      best_current_iteration_score <- iter_scores$score[best_idx]
      best_candidate_this_iteration <- iter_scores$candidate_idx[best_idx]
    }
    
    # Sort iter_scores for review (optional)
    # iter_scores <- iter_scores[order(-iter_scores$score), ]
    # if(verbose && nrow(iter_scores) > 0) {
    #    print(head(iter_scores))
    # }

    if (is.null(best_candidate_this_iteration) || best_current_iteration_score == -Inf) {
      if (verbose) message(sprintf("  MRA Iteration %d: No suitable candidate found to improve score. Stopping.", i))
      break
    }

    selected_anchors <- sort(unique(c(selected_anchors, best_candidate_this_iteration)))
    candidate_pool <- setdiff(candidate_pool, best_candidate_this_iteration)

    if (verbose) {
      sel_metrics <- .evaluate_anchor_set(selected_anchors) # Re-evaluate for the chosen set
      message(sprintf("  MRA Iteration %d: Selected anchor %d. Total selected: %d. Score: %.4f (InvKappa: %.4f, Disp: %.4f)",
                      i, best_candidate_this_iteration, length(selected_anchors), 
                      sel_metrics$score, ifelse(is.finite(sel_metrics$kappa), 1/sel_metrics$kappa, 0), sel_metrics$dispersion))
    }
  }

  if (verbose) {
    final_metrics <- .evaluate_anchor_set(selected_anchors)
    message_stage(sprintf("MRA-Select Finished. Selected %d anchors. Final Score: %.4f (InvKappa: %.4f, Disp: %.4f)", 
                          length(selected_anchors), final_metrics$score, 
                          ifelse(is.finite(final_metrics$kappa), 1/final_metrics$kappa, 0), 
                          final_metrics$dispersion), interactive_only = TRUE)
  }
  
  if (length(selected_anchors) == 0 && m_target > 0 && verbose){
      warning("MRA-Select did not select any anchors.")
  }

  return(sort(selected_anchors))
}

# Helper for null or default in function (already in hatsa_qc_plots.R, define locally or ensure source)
# Consider moving to a utils.R file if used commonly
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
} 