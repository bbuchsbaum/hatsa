#' Helper functions for MDS plotting
#'
#' These internal utilities support `plot_mds_spd_subjects` by
#' computing distance matrices, running classical MDS and preparing
#' data frames for plotting.
#'
#' @keywords internal
NULL

#' Compute Riemannian Distance Matrix
#'
#' Wrapper around `riemannian_distance_matrix_spd` used by
#' `plot_mds_spd_subjects`.
#'
#' @param projector_object A `hatsa_projector` or `task_hatsa_projector`.
#' @param spd_representation_type Character string indicating the SPD
#'   representation to use.
#' @param dist_mat_options List of additional arguments for
#'   `riemannian_distance_matrix_spd`.
#' @param verbose Logical controlling message output.
#'
#' @return A numeric matrix or `NULL` on failure.
#' @keywords internal
compute_mds_distance_matrix <- function(projector_object,
                                         spd_representation_type,
                                         dist_mat_options = list(),
                                         verbose = FALSE) {
  dist_args <- dist_mat_options
  dist_args$object <- projector_object
  dist_args$type <- spd_representation_type
  if (is.null(dist_args$verbose)) dist_args$verbose <- verbose

  tryCatch(
    do.call(hatsa::riemannian_distance_matrix_spd, dist_args),
    error = function(e) {
      warning(sprintf("Riemannian distance matrix computation failed: %s", e$message))
      NULL
    }
  )
}

#' Run `cmdscale` Safely
#'
#' Performs classical MDS with error handling.
#'
#' @param dist_matrix Symmetric distance matrix.
#' @param k_mds Number of dimensions to return.
#' @param add Logical `add` argument for `stats::cmdscale`.
#'
#' @return Output of `stats::cmdscale` or `NULL` on failure.
#' @keywords internal
run_cmdscale_safe <- function(dist_matrix, k_mds = 2, add = TRUE) {
  tryCatch(
    stats::cmdscale(as.dist(dist_matrix), k = k_mds, eig = TRUE, add = add),
    error = function(e) {
      warning(sprintf("MDS computation (cmdscale) failed: %s", e$message))
      NULL
    }
  )
}

#' Prepare Data for MDS Plot
#'
#' Creates a data frame of MDS coordinates merged with optional
#' subject information for plotting.
#'
#' @param mds_fit Result from `run_cmdscale_safe`.
#' @param valid_subject_indices Indices of subjects included in the analysis.
#' @param valid_subjects_mask Logical vector identifying valid subjects in the
#'   original distance matrix.
#' @param n_total_subjects Total number of subjects in the projector object.
#' @param subject_info Optional data frame with one row per subject.
#' @param plot_labels Logical indicating if labels should be prepared.
#'
#' @return A data frame suitable for `ggplot2`.
#' @keywords internal
prepare_mds_plot_df <- function(mds_fit,
                                valid_subject_indices,
                                valid_subjects_mask,
                                n_total_subjects,
                                subject_info = NULL,
                                plot_labels = FALSE) {
  mds_coords <- as.data.frame(mds_fit$points)
  colnames(mds_coords) <- paste0("Dim", seq_len(ncol(mds_coords)))
  mds_coords$subject_idx_original <- valid_subject_indices

  plot_df <- mds_coords
  plot_df$subject_label_plot <- paste0("S", plot_df$subject_idx_original)

  if (!is.null(subject_info) && inherits(subject_info, "data.frame")) {
    if (nrow(subject_info) == n_total_subjects) {
      subject_info_valid <- subject_info[valid_subjects_mask, , drop = FALSE]
      subject_info_valid$subject_idx_original <- valid_subject_indices
      plot_df <- merge(plot_df, subject_info_valid,
                       by = "subject_idx_original", all.x = TRUE)
      if ("subject_label" %in% colnames(plot_df) && plot_labels) {
        plot_df$subject_label_plot <- plot_df$subject_label
      }
    } else {
      warning("`subject_info` provided but row count does not match N_subjects. It will be ignored for merging specifics, but used for column name checks.")
    }
  }

  plot_df
}
