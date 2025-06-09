#' Convenience Accessors for HATSA Results
#'
#' @description
#' Functions to extract commonly needed components from HATSA results
#' in a user-friendly format.
#'
#' @name hatsa-accessors
#' @param x A hatsa_projector or task_hatsa_projector object
#' @param subject Optional subject index or "all"
#' @param ... Additional arguments
NULL

#' @describeIn hatsa-accessors Get aligned data for one or all subjects
#' @export
get_aligned_data <- function(x, subject = "all") {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  
  if (identical(subject, "all")) {
    # Return list of aligned data for all subjects
    lapply(seq_along(x$block_indices), function(i) {
      project_block(x, block = i)
    })
  } else {
    # Return aligned data for specific subject
    project_block(x, block = subject)
  }
}

#' @describeIn hatsa-accessors Get the group template (anchor alignment)
#' @export
get_template <- function(x) {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  x$T_anchor_final
}

#' @describeIn hatsa-accessors Get rotation matrices
#' @export
get_rotations <- function(x, subject = "all") {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  
  if (identical(subject, "all")) {
    x$R_final_list
  } else {
    x$R_final_list[[subject]]
  }
}

#' @describeIn hatsa-accessors Get anchor indices used
#' @export
get_anchor_indices <- function(x) {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  x$parameters$anchor_indices
}

#' @describeIn hatsa-accessors Get alignment quality metrics
#' @export
get_quality_metrics <- function(x) {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  
  # Extract key quality metrics
  summ <- summary(x)
  
  list(
    mean_anchor_error = summ$mean_anchor_alignment_error,
    rotation_dispersion = summ$rotation_metric_mitteroecker,
    eigenvalue_gaps = x$Lambda_original_gaps_list,
    n_subjects = x$parameters$N_subjects,
    n_components = x$parameters$k,
    n_voxels = x$parameters$V_p
  )
}

#' @describeIn hatsa-accessors Get original (unaligned) basis functions
#' @export  
get_original_basis <- function(x, subject = "all") {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  
  if (identical(subject, "all")) {
    x$U_original_list
  } else {
    x$U_original_list[[subject]]
  }
}

#' @describeIn hatsa-accessors Get eigenvalues from spectral decomposition
#' @export
get_eigenvalues <- function(x, subject = "all") {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  
  if (identical(subject, "all")) {
    x$Lambda_original_list
  } else {
    x$Lambda_original_list[[subject]]
  }
}

#' Quick Summary of HATSA Results
#'
#' @description
#' Provides a concise overview of HATSA alignment results with
#' key metrics and diagnostic information.
#'
#' @param x A hatsa_projector object
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns a list of summary statistics
#' @export
hatsa_summary <- function(x, ...) {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  
  metrics <- get_quality_metrics(x)
  
  cat("HATSA Alignment Results\n")
  cat(rep("=", 40), "\n", sep = "")
  cat(sprintf("Subjects: %d\n", metrics$n_subjects))
  cat(sprintf("Voxels: %d\n", metrics$n_voxels))
  cat(sprintf("Components: %d\n", metrics$n_components))
  cat(sprintf("Anchors: %d\n", length(get_anchor_indices(x))))
  cat("\n")
  cat("Alignment Quality:\n")
  cat(sprintf("  Mean anchor error: %.4f\n", metrics$mean_anchor_error))
  cat(sprintf("  Rotation dispersion: %.4f\n", metrics$rotation_dispersion))
  cat("\n")
  cat("Component Separation:\n")
  gaps <- unlist(metrics$eigenvalue_gaps)
  if (length(gaps) > 0) {
    cat(sprintf("  Mean eigenvalue gap: %.4f\n", mean(gaps)))
    cat(sprintf("  Min eigenvalue gap: %.4f\n", min(gaps)))
  }
  
  invisible(metrics)
}

#' Plot HATSA Results
#'
#' @description
#' Creates diagnostic plots for HATSA alignment results.
#'
#' @param x A hatsa_projector object
#' @param type Plot type: "eigenvalues", "rotation_quality", or "alignment"
#' @param ... Additional plotting parameters
#'
#' @export
plot_hatsa <- function(x, type = c("eigenvalues", "rotation_quality", "alignment"), ...) {
  if (!inherits(x, "hatsa_projector")) {
    stop("x must be a hatsa_projector object")
  }
  
  type <- match.arg(type)
  
  if (type == "eigenvalues") {
    # Delegate to existing plot_k_stability_hatsa
    plot_k_stability_hatsa(
      Lambda_list = x$Lambda_original_list,
      Lambda_gaps_list = x$Lambda_original_gaps_list,
      ...
    )
  } else if (type == "rotation_quality") {
    # Use existing MDS plot
    plot_mds_spd_subjects(
      proj_list = list(x),
      proj_names = "HATSA",
      ...
    )
  } else {
    # Placeholder for alignment quality plot
    message("Alignment plot not yet implemented")
  }
}