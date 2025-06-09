#' Suggest HATSA Parameters
#'
#' @description
#' Analyzes your data and suggests appropriate HATSA parameters based on
#' data characteristics such as dimensions, temporal structure, and noise levels.
#'
#' @param data List of subject data matrices (time × voxels)
#' @param task_data Optional list of task activation matrices
#' @param verbose Logical, whether to print explanations (default: TRUE)
#'
#' @return A list of suggested parameters with explanations
#'
#' @examples
#' \dontrun{
#' # Get parameter suggestions
#' params <- hatsa_suggest(subject_data)
#' 
#' # Use suggested parameters
#' result <- hatsa(subject_data, 
#'                 components = params$components,
#'                 preset = params$preset)
#' }
#'
#' @export
hatsa_suggest <- function(data, task_data = NULL, verbose = TRUE) {
  
  # Basic data characteristics
  n_subjects <- length(data)
  n_timepoints <- nrow(data[[1]])
  n_voxels <- ncol(data[[1]])
  
  # Analyze temporal autocorrelation (proxy for smoothness)
  temporal_smooth <- mean(sapply(data[1:min(3, n_subjects)], function(X) {
    mean(cor(X[-1,], X[-nrow(X),], use = "pairwise.complete.obs"), na.rm = TRUE)
  }))
  
  # Estimate noise level from high-frequency components
  noise_level <- estimate_noise_level(data[[1]])
  
  suggestions <- list()
  
  # Suggest number of components
  if (n_voxels < 100) {
    suggestions$components <- min(20, round(n_voxels * 0.3))
    comp_reason <- "Small number of voxels - using 30% as components"
  } else if (n_voxels < 1000) {
    suggestions$components <- min(30, round(n_voxels * 0.1))
    comp_reason <- "Medium voxel count - using 10% as components"
  } else {
    suggestions$components <- min(50, round(n_voxels * 0.05))
    comp_reason <- "Large voxel count - using 5% as components (capped at 50)"
  }
  
  # Suggest connectivity parameters
  if (n_voxels < 500) {
    suggestions$k_conn_pos <- 5
    suggestions$k_conn_neg <- 5
    conn_reason <- "Small dataset - using sparse connectivity"
  } else if (temporal_smooth > 0.7) {
    suggestions$k_conn_pos <- 15
    suggestions$k_conn_neg <- 15
    conn_reason <- "Smooth temporal structure - using denser connectivity"
  } else {
    suggestions$k_conn_pos <- 10
    suggestions$k_conn_neg <- 10
    conn_reason <- "Standard connectivity parameters"
  }
  
  # Suggest preset based on data size and quality
  if (n_subjects < 10 || n_timepoints < 100) {
    suggestions$preset <- "accurate"
    preset_reason <- "Small dataset - using accurate preset for stability"
  } else if (n_subjects > 50 && n_voxels > 5000) {
    suggestions$preset <- "fast"
    preset_reason <- "Large dataset - using fast preset for efficiency"
  } else {
    suggestions$preset <- "default"
    preset_reason <- "Medium-sized dataset - using default preset"
  }
  
  # Suggest number of refinement iterations
  if (noise_level > 0.5) {
    suggestions$n_refine <- 5
    refine_reason <- "High noise level - using more refinement iterations"
  } else if (n_subjects < 20) {
    suggestions$n_refine <- 3
    refine_reason <- "Small sample size - using standard refinement"
  } else {
    suggestions$n_refine <- 2
    refine_reason <- "Clean data with many subjects - minimal refinement needed"
  }
  
  # Task-specific suggestions
  if (!is.null(task_data)) {
    n_conditions <- ncol(task_data[[1]])
    
    if (n_conditions > 20) {
      suggestions$task_method <- "gev"
      task_reason <- "Many task conditions - GEV method recommended"
      suggestions$lambda_max_thresh <- 0.9
    } else if (n_conditions < 5) {
      suggestions$task_method <- "augmented"
      task_reason <- "Few task conditions - augmentation method recommended"
      suggestions$augmentation_weight <- 0.3
    } else {
      suggestions$task_method <- "blend"
      task_reason <- "Moderate task conditions - blend method recommended"
      suggestions$lambda_blend <- 0.15
    }
  }
  
  # Anchor selection suggestions
  min_anchors <- suggestions$components * 2
  max_anchors <- min(200, round(n_voxels * 0.2))
  suggestions$n_anchors <- min(max(min_anchors, 50), max_anchors)
  anchor_reason <- sprintf("Using %d anchors (%.1f%% of voxels)", 
                          suggestions$n_anchors,
                          100 * suggestions$n_anchors / n_voxels)
  
  # Print explanations if verbose
  if (verbose) {
    cat("HATSA Parameter Suggestions\n")
    cat(rep("=", 50), "\n", sep = "")
    cat(sprintf("Data: %d subjects, %d timepoints, %d voxels\n", 
                n_subjects, n_timepoints, n_voxels))
    cat(sprintf("Temporal smoothness: %.2f\n", temporal_smooth))
    cat(sprintf("Estimated noise level: %.2f\n\n", noise_level))
    
    cat("Suggested parameters:\n")
    cat(sprintf("  components: %d (%s)\n", suggestions$components, comp_reason))
    cat(sprintf("  preset: '%s' (%s)\n", suggestions$preset, preset_reason))
    cat(sprintf("  connectivity: k_pos=%d, k_neg=%d (%s)\n", 
                suggestions$k_conn_pos, suggestions$k_conn_neg, conn_reason))
    cat(sprintf("  refinement: %d iterations (%s)\n", 
                suggestions$n_refine, refine_reason))
    cat(sprintf("  anchors: %s\n", anchor_reason))
    
    if (!is.null(task_data)) {
      cat(sprintf("\n  Task method: '%s' (%s)\n", 
                  suggestions$task_method, task_reason))
    }
    
    cat("\nExample usage:\n")
    if (is.null(task_data)) {
      cat("  result <- hatsa(data, components = ", suggestions$components,
          ", preset = '", suggestions$preset, "')\n", sep = "")
    } else {
      cat("  result <- hatsa_task(data, task_data, method = '", 
          suggestions$task_method, "')\n", sep = "")
    }
  }
  
  suggestions
}

#' Estimate noise level from data
#' @keywords internal
estimate_noise_level <- function(X) {
  # Simple noise estimation from high-frequency components
  if (nrow(X) < 10) return(0.5)  # Default for very short time series
  
  # Use median absolute deviation of differences
  diffs <- diff(X)
  mad_vals <- apply(diffs, 2, function(x) mad(x, na.rm = TRUE))
  median(mad_vals, na.rm = TRUE) / sd(as.vector(X), na.rm = TRUE)
}

#' Validate HATSA Parameters
#'
#' @description
#' Checks if the provided parameters are valid for your data and suggests
#' corrections if needed.
#'
#' @param data List of subject data matrices
#' @param anchor_indices Anchor indices to validate
#' @param spectral_rank_k Number of components
#' @param k_conn_pos Number of positive edges
#' @param k_conn_neg Number of negative edges
#' @param ... Other parameters to validate
#'
#' @return Invisible TRUE if valid, otherwise stops with informative error
#' @export
hatsa_validate_params <- function(data, 
                                  anchor_indices,
                                  spectral_rank_k,
                                  k_conn_pos = 10,
                                  k_conn_neg = 10,
                                  ...) {
  
  n_voxels <- ncol(data[[1]])
  n_subjects <- length(data)
  
  # Component checks
  if (spectral_rank_k > n_voxels * 0.5) {
    stop(sprintf(
      "spectral_rank_k (%d) is too large for %d voxels. Try %d or fewer components.",
      spectral_rank_k, n_voxels, round(n_voxels * 0.3)
    ))
  }
  
  # Anchor checks  
  if (length(anchor_indices) < spectral_rank_k) {
    stop(sprintf(
      "Not enough anchors (%d) for %d components. Need at least %d anchors.",
      length(anchor_indices), spectral_rank_k, spectral_rank_k
    ))
  }
  
  if (length(anchor_indices) > n_voxels * 0.5) {
    warning(sprintf(
      "Using %d anchors (%.0f%% of voxels) is unusually high. Consider using %d instead.",
      length(anchor_indices), 100 * length(anchor_indices) / n_voxels,
      round(n_voxels * 0.2)
    ))
  }
  
  # Connectivity checks
  max_edges <- min(50, round(n_voxels * 0.1))
  if (k_conn_pos > max_edges || k_conn_neg > max_edges) {
    warning(sprintf(
      "k_conn values (%d, %d) are high for %d voxels. Consider using (%d, %d).",
      k_conn_pos, k_conn_neg, n_voxels,
      min(k_conn_pos, max_edges), min(k_conn_neg, max_edges)
    ))
  }
  
  invisible(TRUE)
}