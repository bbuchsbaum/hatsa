#' HATSA Quality Control and Stability Plots
#'
#' This file contains functions for generating QC plots, such as k-stability plots,
#' for HATSA results.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap theme_bw labs scale_x_continuous ggtitle
#' @importFrom stats sd na.omit

# Ensure stats is available for sd, na.omit if not base
# Ensure ggplot2 is available for plotting functions

#' Plot HATSA k-Stability Metrics
#'
#' Generates plots of stability metrics versus the number of components (k)
#' from a list of HATSA results obtained with varying k. This helps in selecting
#' an optimal k.
#'
#' @param projector_list_over_k A list of \code{hatsa_projector} objects.
#'   Each object should be the result of a HATSA run with a different 'k'.
#'   The list can be named with k values, or k will be extracted from each object.
#' @param metrics_to_plot A character vector specifying which metrics to plot.
#'   Possible values are "Hk" (Alignment Homogeneity based on Riemannian dispersion
#'   of covariance of aligned coefficients) and/or "CV_eigen_cov" (Mean Coefficient
#'   of Variation of eigenvalues of covariance of aligned coefficients).
#'   Default is \code{c("Hk", "CV_eigen_cov")}.
#' @param Hk_dispersion_stat Character, either "mean" (default) or "median".
#'   Specifies which statistic from \code{riemannian_dispersion_spd} to use for Hk.
#' @param spd_cov_coeffs_options A list of additional arguments to pass to
#'   \code{hatsa::get_spd_representations(..., type = "cov_coeffs")}.
#'   For example, \code{list(spd_regularize_epsilon = 1e-6)}.
#' @param Hk_options A list of additional arguments to pass to
#'   \code{hatsa::riemannian_dispersion_spd(...)}.
#'   For example, \code{list(tangent_metric = "logeuclidean", use_geometric_median = FALSE)}.
#'   Note that `use_geometric_median` will be overridden by `Hk_dispersion_stat == "median"`.
#' @param verbose Logical, if TRUE, prints progress messages. Default is FALSE.
#'
#' @return A \code{ggplot} object containing the stability plots, or NULL if
#'   no valid metrics can be computed.
#'
#' @export
#' @examples
#' # Conceptual example:
#' # Assume projector_k3, projector_k4, projector_k5 are hatsa_projector objects
#' # k_list <- list("3" = projector_k3, "4" = projector_k4, "5" = projector_k5)
#' # plot_k_stability_hatsa(k_list)
#'
#' # To run if example data were available:
#' # if (requireNamespace("ggplot2", quietly = TRUE) &&
#' #     exists("generate_synthetic_hatsa_output_list")) {
#' #   example_projector_list <- generate_synthetic_hatsa_output_list(
#' #      k_values = c(3, 4, 5, 6), N_subjects = 5, V_p = 20
#' #   )
#' #   plot_k_stability_hatsa(example_projector_list)
#' # }
plot_k_stability_hatsa <- function(projector_list_over_k,
                                   metrics_to_plot = c("Hk", "CV_eigen_cov"),
                                   Hk_dispersion_stat = "mean",
                                   spd_cov_coeffs_options = list(),
                                   Hk_options = list(),
                                   verbose = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!is.list(projector_list_over_k) || length(projector_list_over_k) == 0) {
    stop("`projector_list_over_k` must be a non-empty list of hatsa_projector objects.")
  }
  if (!Hk_dispersion_stat %in% c("mean", "median")) {
    stop("`Hk_dispersion_stat` must be either 'mean' or 'median'.")
  }

  results_df_list <- list()

  for (i in seq_along(projector_list_over_k)) {
    proj_obj <- projector_list_over_k[[i]]

    if (!inherits(proj_obj, "hatsa_projector")) {
      warning(sprintf("Item %d in projector_list_over_k is not a hatsa_projector object. Skipping.", i))
      next
    }
    current_k <- proj_obj$parameters$k
    if (is.null(current_k)) {
      warning(sprintf("Could not determine k for item %d. Skipping.", i))
      next
    }
    if (verbose) message_stage(sprintf("Processing for k = %d...", current_k))

    # --- Calculate Hk (Alignment Homogeneity) ---
    if ("Hk" %in% metrics_to_plot) {
      if (verbose) message_stage("  Calculating Hk (Riemannian Dispersion)...")
      
      current_Hk_options <- Hk_options
      current_Hk_options$object <- proj_obj
      current_Hk_options$type <- "cov_coeffs" # Hk is specifically for cov_coeffs
      current_Hk_options$use_geometric_median <- (Hk_dispersion_stat == "median")
      # Merge get_spd_representations options if needed by riemannian_dispersion_spd indirectly
      # `riemannian_dispersion_spd` itself calls `get_spd_representations`
      # So, spd_cov_coeffs_options should be part of Hk_options if they affect SPD matrices
      # for dispersion (e.g. k_conn_params if type was fc_conn, but here it's cov_coeffs)
      # For cov_coeffs, only spd_regularize_epsilon from spd_cov_coeffs_options might be relevant via `...`
      # to `get_spd_representations` call within `riemannian_dispersion_spd`
      # It might be cleaner if riemannian_dispersion_spd explicitly takes spd_options
      if (!is.null(spd_cov_coeffs_options$spd_regularize_epsilon)) {
           current_Hk_options$spd_regularize_epsilon <- spd_cov_coeffs_options$spd_regularize_epsilon
      }


      disp_results <- tryCatch(
        do.call(hatsa::riemannian_dispersion_spd, current_Hk_options),
        error = function(e) {
          warning(sprintf("Hk calculation failed for k = %d: %s", current_k, e$message))
          NULL
        }
      )
      
      Hk_value <- NA_real_
      if (!is.null(disp_results)) {
        Hk_value <- disp_results[[paste0(Hk_dispersion_stat, "_dispersion")]]
      }
      results_df_list[[length(results_df_list) + 1]] <- data.frame(
        k = current_k,
        metric_name = paste0("Hk (", Hk_dispersion_stat, " disp.)"),
        value = Hk_value
      )
    }

    # --- Calculate CV_eigen_cov (Mean CV of eigenvalues of Cov_coeffs) ---
    if ("CV_eigen_cov" %in% metrics_to_plot) {
      if (verbose) message_stage("  Calculating CV_eigen_cov...")
      
      current_spd_options <- spd_cov_coeffs_options
      current_spd_options$object <- proj_obj
      current_spd_options$type <- "cov_coeffs"
      
      spd_list <- tryCatch(
        do.call(hatsa::get_spd_representations, current_spd_options),
        error = function(e) {
          warning(sprintf("Failed to get SPD representations for k = %d for CV_eigen_cov: %s", current_k, e$message))
          NULL
        }
      )

      mean_cv_value <- NA_real_
      if (!is.null(spd_list)) {
        valid_spd_list <- Filter(Negate(is.null), spd_list)
        if (length(valid_spd_list) > 0) {
          cvs_subject <- vapply(valid_spd_list, function(spd_matrix) {
            if (is.null(spd_matrix) || !is.matrix(spd_matrix) || any(dim(spd_matrix) == 0)) return(NA_real_)
            # Eigenvalues of a kxk covariance matrix. k here is proj_obj$parameters$k.
            if (proj_obj$parameters$k == 0) return(NA_real_)

            eig_vals <- tryCatch(
            eigen(spd_matrix, symmetric = TRUE, only.values = TRUE)$values,
            error = function(e_eig) {warning(sprintf("Eigen decomposition failed for a cov_coeff matrix (k=%d).", current_k)); NULL}
          )
          
          if (is.null(eig_vals) || length(eig_vals) == 0) return(NA_real_)
          
          # Remove numerically zero or negative eigenvalues before CV calculation for stability
          # Though eigenvalues of a (regularized) SPD matrix should be positive.
          eig_vals_clean <- eig_vals[eig_vals > .Machine$double.eps^0.75] # A bit more than double.eps

          if (length(eig_vals_clean) < 1) return(NA_real_) # Needs at least one positive eigenvalue
          if (length(eig_vals_clean) == 1) return(NA_real_) # SD of 1 value is NA, so CV is NA. Correct.
                                                         # Or, define CV as 0 if only one eigenvalue? No, NA is better.

          mean_e <- mean(eig_vals_clean, na.rm = TRUE)
          sd_e <- stats::sd(eig_vals_clean, na.rm = TRUE)

            if (is.na(sd_e) || is.na(mean_e) || abs(mean_e) < .Machine$double.eps) { # check mean_e for near zero
              return(NA_real_)
            }
            return(sd_e / mean_e)
          }, numeric(1))
          mean_cv_value <- mean(stats::na.omit(cvs_subject), na.rm = TRUE)
          if (is.nan(mean_cv_value)) mean_cv_value <- NA_real_ # Handle case where all cvs_subject are NA
        }
      }
      results_df_list[[length(results_df_list) + 1]] <- data.frame(
        k = current_k,
        metric_name = "CV_eigen_cov",
        value = mean_cv_value
      )
    }
  }

  if (length(results_df_list) == 0) {
    if (verbose) message_stage("No metrics could be computed.")
    return(NULL)
  }

  plot_data <- do.call(rbind, results_df_list)
  
  if (nrow(plot_data) == 0 || all(is.na(plot_data$value))) {
    if (verbose) message_stage("No valid data points to plot.")
    return(NULL)
  }
  
  # Ensure k is numeric for plotting scale
  plot_data$k <- as.numeric(plot_data$k)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = value, group = metric_name, color = metric_name)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ metric_name, scales = "free_y", ncol = 1) +
    ggplot2::scale_x_continuous(breaks = unique(plot_data$k)) + # Ensure all k values are ticks
    ggplot2::labs(
      title = "HATSA k-Stability Analysis",
      x = "Number of Components (k)",
      y = "Metric Value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") # Facets make legend redundant

  return(p)
}

# Helper function for generating example data (conceptual)
# To make the example runnable, one might need a synthetic data generator.
# This is a placeholder to illustrate how projector_list_over_k would look.
# generate_synthetic_hatsa_output_list <- function(k_values = c(3,4,5), N_subjects = 5, V_p = 20) {
#   # This function would need to create valid hatsa_projector objects
#   # using synthetic data and calls to run_hatsa_core and hatsa_projector constructor.
#   # For now, it's a conceptual placeholder.
#   # For proper testing, one would mock hatsa_projector objects or use actual test data.
#   
#   # Example structure:
#   # k_list <- lapply(k_values, function(k_val) {
#   #   # ... generate or load a hatsa_projector object for k_val ...
#   #   # For example:
#   #   params <- list(k = k_val, N_subjects = N_subjects, V_p = V_p, method = "hatsa_core", anchor_indices = 1:5)
#   #   core_res <- list(
#   #     U_aligned_list = replicate(N_subjects, matrix(rnorm(V_p*k_val), V_p, k_val), simplify=FALSE),
#   #     R_final_list = replicate(N_subjects, diag(k_val), simplify=FALSE),
#   #     U_original_list = replicate(N_subjects, matrix(rnorm(V_p*k_val), V_p, k_val), simplify=FALSE),
#   #     Lambda_original_list = replicate(N_subjects, sort(runif(k_val, 0.1, 1), decreasing=TRUE), simplify=FALSE),
#   #     Lambda_original_gaps_list = replicate(N_subjects, {
#   #        lams <- sort(runif(k_val, 0.1, 1), decreasing=TRUE); if(k_val>1) (lams[1:(k_val-1)] - lams[2:k_val])/lams[2:k_val] else numeric(0)
#   #        }, simplify=FALSE),
#   #     T_anchor_final = matrix(rnorm(min(5,V_p)*k_val), min(5,V_p), k_val)
#   #   )
#   #   hatsa::hatsa_projector(core_res, params)
#   # })
#   # names(k_list) <- as.character(k_values)
#   # return(k_list)
#   warning("generate_synthetic_hatsa_output_list is a conceptual placeholder and does not produce real data.")
#   return(list())
# }


#' Plot Multidimensional Scaling (MDS) of Subjects based on SPD Matrix Distances
#'
#' Computes Riemannian distances between subjects based on their SPD matrix
#' representations, performs classical MDS, and plots the subjects in a
#' low-dimensional space (typically 2D).
#'
#' @param projector_object A \code{hatsa_projector} or \code{task_hatsa_projector} object.
#' @param k_mds Integer, the number of MDS dimensions to compute (default: 2).
#' @param spd_representation_type Character string, the type of SPD representation
#'   to use for distance calculation (e.g., "cov_coeffs", "fc_conn").
#'   This is passed to \code{hatsa::riemannian_distance_matrix_spd}.
#' @param dist_mat_options A list of additional arguments to pass to
#'   \code{hatsa::riemannian_distance_matrix_spd}. This can include arguments like
#'   \code{spd_metric}, \code{subject_data_list} (if needed for the chosen type),
#'   \code{k_conn_params}, \code{spd_regularize_epsilon}, \code{verbose}, etc.
#' @param subject_info Optional. A data frame with \code{N_subjects} rows.
#'   If provided, it can contain a column named \code{subject_label} for text labels on
#'   the plot, and other columns that can be mapped to ggplot aesthetics (e.g.,
#'   a column named by \code{color_by_column} or \code{shape_by_column}).
#'   Row names of \code{subject_info} should correspond to subject indices (1 to N) or
#'   be actual subject IDs if the distance matrix has them.
#' @param color_by_column Character string. If \code{subject_info} is provided, the name
#'   of the column in \code{subject_info} to use for coloring points.
#' @param shape_by_column Character string. If \code{subject_info} is provided, the name
#'   of the column in \code{subject_info} to use for point shapes.
#' @param plot_labels Logical, whether to plot subject labels near points. Requires
#'   a \code{subject_label} column in \code{subject_info} or uses default labels.
#'   Consider using the \code{ggrepel} package for better label placement if many
#'   points overlap (not directly implemented here to reduce dependencies).
#' @param cmdscale_add Logical, the \code{add} argument for \code{stats::cmdscale}
#'   (default: TRUE). Useful if distances are not perfectly Euclidean.
#' @param verbose Logical, if TRUE, prints progress messages. Default is FALSE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\code{plot}: A \code{ggplot} object for the MDS plot.}
#'     \item{\code{mds_results}: The output from \code{stats::cmdscale}, including
#'           coordinates (\code{points}), eigenvalues (\code{eig}), etc.}
#'     \item{\code{distance_matrix}: The computed Riemannian distance matrix.}
#'     \item{\code{valid_subject_indices}: Indices of subjects included in MDS (after NA removal).}
#'   }
#'   Returns NULL if critical steps fail (e.g., distance matrix computation).
#'
#' @export
#' @importFrom stats cmdscale
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_hline geom_vline labs theme_bw
#' @examples
#' # Conceptual example, assuming 'proj_obj' is a hatsa_projector object
#' # and subject_covariates is a data frame with N_subjects rows
#' # and columns "ID" (for labels), "Group" (for color).
#'
#' # if (requireNamespace("ggplot2", quietly = TRUE) &&
#' #     exists("generate_synthetic_hatsa_output") && # Assuming a single projector gen
#' #     exists("hatsa_projector")) { # Ensure constructor is available
#' #
#' #   # Generate a single projector object
#' #   N_subj_example <- 10
#' #   V_p_example <- 30
#' #   k_example <- 5
#' #   proj_params <- list(k = k_example, N_subjects = N_subj_example, V_p = V_p_example,
#' #                       method="hatsa_core", anchor_indices = 1:5)
#' #   proj_core_res <- list(
#' #      U_aligned_list = replicate(N_subj_example, matrix(rnorm(V_p_example*k_example), V_p_example, k_example), simplify=FALSE),
#' #      R_final_list = replicate(N_subj_example, diag(k_example), simplify=FALSE),
#' #      U_original_list = replicate(N_subj_example, matrix(rnorm(V_p_example*k_example), V_p_example, k_example), simplify=FALSE),
#' #      Lambda_original_list = replicate(N_subj_example, sort(runif(k_example, 0.1, 1), decreasing=TRUE), simplify=FALSE),
#' #      Lambda_original_gaps_list = replicate(N_subj_example, {
#' #         lams <- sort(runif(k_example, 0.1, 1), decreasing=TRUE); if(k_example>1) (lams[1:(k_example-1)] - lams[2:k_example])/lams[2:k_example] else numeric(0)
#' #         }, simplify=FALSE),
#' #      T_anchor_final = matrix(rnorm(min(5,V_p_example)*k_example), min(5,V_p_example), k_example)
#' #    )
#' #   proj_obj <- hatsa::hatsa_projector(proj_core_res, proj_params)
#' #
#' #   # Example subject_info
#' #   subject_covariates <- data.frame(
#' #     subject_label = paste0("S", 1:N_subj_example),
#' #     Group = factor(rep(c("A", "B"), each = N_subj_example / 2)),
#' #     stringsAsFactors = FALSE
#' #   )
#' #
#' #   # Run MDS plot
#' #   mds_plot_result <- plot_mds_spd_subjects(
#' #     projector_object = proj_obj,
#' #     spd_representation_type = "cov_coeffs",
#' #     dist_mat_options = list(spd_metric = "logeuclidean", spd_regularize_epsilon = 1e-6),
#' #     subject_info = subject_covariates,
#' #     color_by_column = "Group",
#' #     plot_labels = TRUE,
#' #     verbose = TRUE
#' #   )
#' #   if (!is.null(mds_plot_result)) print(mds_plot_result$plot)
#' # }
plot_mds_spd_subjects <- function(projector_object,
                                  k_mds = 2,
                                  spd_representation_type = "cov_coeffs",
                                  dist_mat_options = list(),
                                  subject_info = NULL,
                                  color_by_column = NULL,
                                  shape_by_column = NULL,
                                  plot_labels = FALSE,
                                  cmdscale_add = TRUE,
                                  verbose = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }
  if (!inherits(projector_object, "hatsa_projector") && !inherits(projector_object, "task_hatsa_projector")) {
    stop("`projector_object` must be a hatsa_projector or task_hatsa_projector object.")
  }
  if (projector_object$parameters$N_subjects == 0) {
    warning("No subjects in the projector object. Cannot create MDS plot.")
    return(NULL)
  }

  if (verbose) message_stage(sprintf("Computing Riemannian distance matrix (type: %s)...", spd_representation_type))

  dist_args <- dist_mat_options
  dist_args$object <- projector_object
  dist_args$type <- spd_representation_type
  # Ensure verbose from main function is passed if not in dist_mat_options
  if (is.null(dist_args$verbose)) dist_args$verbose <- verbose 

  dist_matrix_full <- tryCatch(
    do.call(hatsa::riemannian_distance_matrix_spd, dist_args),
    error = function(e) {
      warning(sprintf("Riemannian distance matrix computation failed: %s", e$message))
      NULL
    }
  )

  if (is.null(dist_matrix_full)) {
    return(NULL)
  }
  
  N_total_subjects <- projector_object$parameters$N_subjects
  
  # Identify subjects with all NA distances (completely failed) or rows/cols that are all NA
  # These subjects cannot be included in MDS.
  valid_subjects_mask <- apply(dist_matrix_full, 1, function(row) !all(is.na(row))) &
                         apply(dist_matrix_full, 2, function(col) !all(is.na(col)))
  
  num_valid_subjects <- sum(valid_subjects_mask)

  if (num_valid_subjects < 2) {
    warning("Less than 2 subjects have valid distance data. Cannot perform MDS.")
    return(list(plot = NULL, mds_results = NULL, distance_matrix = dist_matrix_full, valid_subject_indices = which(valid_subjects_mask)))
  }
  if (num_valid_subjects < N_total_subjects) {
     warning(sprintf("%d subjects out of %d had NA distances and were excluded from MDS.", 
                     N_total_subjects - num_valid_subjects, N_total_subjects))
  }
  
  dist_matrix_valid <- dist_matrix_full[valid_subjects_mask, valid_subjects_mask, drop = FALSE]
  valid_subject_indices <- which(valid_subjects_mask)


  if (verbose) message_stage(sprintf("Performing classical MDS on %d subjects for %d dimensions...", num_valid_subjects, k_mds))
  
  # cmdscale expects a dist object or a full symmetric matrix with zeros on diagonal
  # Our riemannian_distance_matrix_spd should return this.
  # However, if there are any remaining NAs after subsetting (should not happen if subsetting is correct),
  # cmdscale will fail. We check for this implicitly by tryCatch.
  mds_fit <- tryCatch(
    stats::cmdscale(as.dist(dist_matrix_valid), k = k_mds, eig = TRUE, add = cmdscale_add),
    error = function(e) {
      warning(sprintf("MDS computation (cmdscale) failed: %s", e$message))
      NULL
    }
  )

  if (is.null(mds_fit)) {
    return(list(plot = NULL, mds_results = NULL, distance_matrix = dist_matrix_full, valid_subject_indices = valid_subject_indices))
  }

  mds_coords <- as.data.frame(mds_fit$points)
  colnames(mds_coords) <- paste0("Dim", 1:ncol(mds_coords))
  mds_coords$subject_idx_original <- valid_subject_indices # Store original index

  # Prepare plotting data by merging with subject_info if provided
  plot_df <- mds_coords
  plot_df$subject_label_plot <- paste0("S", plot_df$subject_idx_original) # Default labels

  if (!is.null(subject_info) && inherits(subject_info, "data.frame")) {
      if(nrow(subject_info) == N_total_subjects) {
          subject_info_valid <- subject_info[valid_subjects_mask, , drop = FALSE]
          subject_info_valid$subject_idx_original <- valid_subject_indices # for merging with mds_coords
          
          # Attempt to merge. If subject_info has rownames that are subject IDs from dist_matrix,
          # and dist_matrix had interpretable rownames/colnames, this logic would need adjustment.
          # Current riemannian_distance_matrix_spd aims for 1:N_subjects in rownames/colnames.
          plot_df <- merge(plot_df, subject_info_valid, by = "subject_idx_original", all.x = TRUE)

          if ("subject_label" %in% colnames(plot_df) && plot_labels) {
             plot_df$subject_label_plot <- plot_df$subject_label
          }
      } else {
          warning("`subject_info` provided but row count does not match N_subjects. It will be ignored for merging specifics, but used for column name checks.")
      }
  }
  
  # Set up aesthetics
  plot_aes <- ggplot2::aes(x = Dim1, y = Dim2)
  if (!is.null(color_by_column) && color_by_column %in% colnames(plot_df)) {
    plot_aes$colour <- ggplot2::aes_string(colour = color_by_column)$colour
  }
  if (!is.null(shape_by_column) && shape_by_column %in% colnames(plot_df)) {
    # Ensure the shaping variable is a factor for discrete shapes
    if(!is.factor(plot_df[[shape_by_column]])) plot_df[[shape_by_column]] <- as.factor(plot_df[[shape_by_column]])
    plot_aes$shape <- ggplot2::aes_string(shape = shape_by_column)$shape
  }


  p <- ggplot2::ggplot(plot_df, plot_aes) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = paste("MDS of Subjects based on", spd_representation_type),
      subtitle = paste0("Riemannian Distances (", dist_mat_options$spd_metric %||% "logeuclidean", 
                       "), MDS k=", k_mds, ", GOF=", sprintf("%.2f", mds_fit$GOF[1])),
      x = "MDS Dimension 1",
      y = "MDS Dimension 2",
      caption = if (num_valid_subjects < N_total_subjects) {
                  sprintf("%d/%d subjects plotted due to NA distances.", num_valid_subjects, N_total_subjects)
                } else { "" }
    ) +
    ggplot2::theme_bw()

  if (plot_labels) {
      # Use ggrepel if available, otherwise geom_text
      if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = subject_label_plot), size = 3, 
                                            max.overlaps = Inf, # try to plot all
                                            box.padding = 0.4) 
      } else {
          p <- p + ggplot2::geom_text(ggplot2::aes(label = subject_label_plot), size = 3, nudge_y = 0.05 * diff(range(plot_df$Dim2, na.rm=TRUE)), check_overlap = TRUE)
      }
  }
  
  # Add % variance explained by axes if k_mds=2 and eig available
  if (k_mds >= 2 && !is.null(mds_fit$eig)) {
      total_pos_eig_variance <- sum(mds_fit$eig[mds_fit$eig > 0])
      if (total_pos_eig_variance > 0) {
          var_dim1 <- (mds_fit$eig[1] / total_pos_eig_variance) * 100
          var_dim2 <- (mds_fit$eig[2] / total_pos_eig_variance) * 100
          p <- p + ggplot2::xlab(sprintf("MDS Dimension 1 (%.1f%%)", var_dim1)) +
                   ggplot2::ylab(sprintf("MDS Dimension 2 (%.1f%%)", var_dim2))
      }
  }


  if (verbose) message_stage("MDS plot generated.")

  return(list(
    plot = p,
    mds_results = mds_fit,
    distance_matrix_computed = dist_matrix_full, # The one before NA removal for MDS
    valid_subject_indices_in_mds = valid_subject_indices
  ))
}

# Helper for null or default
'%||%' <- function(a, b) if (is.null(a)) b else a

# Make sure ggplot2 related imports are at the top of the file or in NAMESPACE
# @importFrom ggplot2 ggplot aes geom_point geom_text geom_hline geom_vline labs theme_bw ggtitle
# @importFrom stats cmdscale
# @importFrom ggrepel geom_text_repel # Optional, but nice for labels


#' Plot Multidimensional Scaling (MDS) of Subjects based on SPD Matrix Distances
#'
#' Computes Riemannian distances between subjects based on their SPD matrix
#' representations, performs classical MDS, and plots the subjects in a
#' low-dimensional space (typically 2D).
#'
#' @param projector_object A \code{hatsa_projector} or \code{task_hatsa_projector} object.
#' @param k_mds Integer, the number of MDS dimensions to compute (default: 2).
#' @param spd_representation_type Character string, the type of SPD representation
#'   to use for distance calculation (e.g., "cov_coeffs", "fc_conn").
#'   This is passed to \code{hatsa::riemannian_distance_matrix_spd}.
#' @param dist_mat_options A list of additional arguments to pass to
#'   \code{hatsa::riemannian_distance_matrix_spd}. This can include arguments like
#'   \code{spd_metric}, \code{subject_data_list} (if needed for the chosen type),
#'   \code{k_conn_params}, \code{spd_regularize_epsilon}, \code{verbose}, etc.
#' @param subject_info Optional. A data frame with \code{N_subjects} rows.
#'   If provided, it can contain a column named \code{subject_label} for text labels on
#'   the plot, and other columns that can be mapped to ggplot aesthetics (e.g.,
#'   a column named by \code{color_by_column} or \code{shape_by_column}).
#'   Row names of \code{subject_info} should correspond to subject indices (1 to N) or
#'   be actual subject IDs if the distance matrix has them.
#' @param color_by_column Character string. If \code{subject_info} is provided, the name
#'   of the column in \code{subject_info} to use for coloring points.
#' @param shape_by_column Character string. If \code{subject_info} is provided, the name
#'   of the column in \code{subject_info} to use for point shapes.
#' @param plot_labels Logical, whether to plot subject labels near points. Requires
#'   a \code{subject_label} column in \code{subject_info} or uses default labels.
#'   Consider using the \code{ggrepel} package for better label placement if many
#'   points overlap (not directly implemented here to reduce dependencies).
#' @param cmdscale_add Logical, the \code{add} argument for \code{stats::cmdscale}
#'   (default: TRUE). Useful if distances are not perfectly Euclidean.
#' @param verbose Logical, if TRUE, prints progress messages. Default is FALSE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\code{plot}: A \code{ggplot} object for the MDS plot.}
#'     \item{\code{mds_results}: The output from \code{stats::cmdscale}, including
#'           coordinates (\code{points}), eigenvalues (\code{eig}), etc.}
#'     \item{\code{distance_matrix}: The computed Riemannian distance matrix.}
#'     \item{\code{valid_subject_indices}: Indices of subjects included in MDS (after NA removal).}
#'   }
#'   Returns NULL if critical steps fail (e.g., distance matrix computation).
#'
#' @export
#' @importFrom stats cmdscale
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_hline geom_vline labs theme_bw
#' @examples
#' # Conceptual example, assuming 'proj_obj' is a hatsa_projector object
#' # and subject_covariates is a data frame with N_subjects rows
#' # and columns "ID" (for labels), "Group" (for color).
#'
#' # if (requireNamespace("ggplot2", quietly = TRUE) &&
#' #     exists("generate_synthetic_hatsa_output") && # Assuming a single projector gen
#' #     exists("hatsa_projector")) { # Ensure constructor is available
#' #
#' #   # Generate a single projector object
#' #   N_subj_example <- 10
#' #   V_p_example <- 30
#' #   k_example <- 5
#' #   proj_params <- list(k = k_example, N_subjects = N_subj_example, V_p = V_p_example,
#' #                       method="hatsa_core", anchor_indices = 1:5)
#' #   proj_core_res <- list(
#' #      U_aligned_list = replicate(N_subj_example, matrix(rnorm(V_p_example*k_example), V_p_example, k_example), simplify=FALSE),
#' #      R_final_list = replicate(N_subj_example, diag(k_example), simplify=FALSE),
#' #      U_original_list = replicate(N_subj_example, matrix(rnorm(V_p_example*k_example), V_p_example, k_example), simplify=FALSE),
#' #      Lambda_original_list = replicate(N_subj_example, sort(runif(k_example, 0.1, 1), decreasing=TRUE), simplify=FALSE),
#' #      Lambda_original_gaps_list = replicate(N_subj_example, {
#' #         lams <- sort(runif(k_example, 0.1, 1), decreasing=TRUE); if(k_example>1) (lams[1:(k_example-1)] - lams[2:k_example])/lams[2:k_example] else numeric(0)
#' #         }, simplify=FALSE),
#' #      T_anchor_final = matrix(rnorm(min(5,V_p_example)*k_example), min(5,V_p_example), k_example)
#' #    )
#' #   proj_obj <- hatsa::hatsa_projector(proj_core_res, proj_params)
#' #
#' #   # Example subject_info
#' #   subject_covariates <- data.frame(
#' #     subject_label = paste0("S", 1:N_subj_example),
#' #     Group = factor(rep(c("A", "B"), each = N_subj_example / 2)),
#' #     stringsAsFactors = FALSE
#' #   )
#' #
#' #   # Run MDS plot
#' #   mds_plot_result <- plot_mds_spd_subjects(
#' #     projector_object = proj_obj,
#' #     spd_representation_type = "cov_coeffs",
#' #     dist_mat_options = list(spd_metric = "logeuclidean", spd_regularize_epsilon = 1e-6),
#' #     subject_info = subject_covariates,
#' #     color_by_column = "Group",
#' #     plot_labels = TRUE,
#' #     verbose = TRUE
#' #   )
#' #   if (!is.null(mds_plot_result)) print(mds_plot_result$plot)
#' # } 