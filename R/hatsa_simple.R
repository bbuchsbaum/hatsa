#' Simple HATSA Interface
#'
#' @description
#' A streamlined interface to HATSA with sensible defaults and automatic parameter selection.
#'
#' @param data List of subject data matrices (time × voxels)
#' @param anchors Either "auto" for automatic selection, or a vector of anchor indices
#' @param components Number of spectral components (default: 20)
#' @param preset Configuration preset: "default", "fast", or "accurate"
#' @param ... Additional parameters passed to run_hatsa_core
#'
#' @return A hatsa_projector object with convenience methods
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- hatsa(subject_data)
#' 
#' # With more components
#' result <- hatsa(subject_data, components = 30)
#' 
#' # Fast mode for exploration
#' result <- hatsa(subject_data, preset = "fast")
#' }
#'
#' @export
hatsa <- function(data, 
                  anchors = "auto",
                  components = 20,
                  preset = c("default", "fast", "accurate"),
                  ...) {
  
  preset <- match.arg(preset)
  
  # Get preset configuration
  config <- hatsa_preset(preset)
  
  # Handle automatic anchor selection
  if (identical(anchors, "auto")) {
    message("Selecting anchors automatically using MRA method...")
    # Use the existing select_anchors_mra function
    anchors <- select_anchors_auto(data, components, config)
  }
  
  # Merge user parameters with preset
  params <- utils::modifyList(config, list(...))
  
  # Call the core function with cleaned up parameters
  do.call(run_hatsa_core, c(
    list(
      subject_data_list = data,
      anchor_indices = anchors,
      spectral_rank_k = components
    ),
    params
  ))
}

#' Task-Informed HATSA Interface
#'
#' @description
#' Simplified interface for task-informed HATSA with automatic method selection.
#'
#' @param data List of subject data matrices (time × voxels)
#' @param task_data List of task activation matrices
#' @param anchors Either "auto" for automatic selection, or a vector of anchor indices
#' @param components Number of spectral components (default: 20)
#' @param method Task incorporation method: "auto", "blend", "gev", or "augmented"
#' @param preset Configuration preset: "default", "fast", or "accurate"
#' @param ... Additional parameters passed to run_task_hatsa
#'
#' @return A task_hatsa_projector object
#'
#' @export
hatsa_task <- function(data,
                       task_data,
                       anchors = "auto",
                       components = 20,
                       method = c("auto", "blend", "gev", "augmented"),
                       preset = c("default", "fast", "accurate"),
                       ...) {
  
  method <- match.arg(method)
  preset <- match.arg(preset)
  
  # Get preset configuration
  config <- hatsa_task_preset(preset, method)
  
  # Auto-select method based on data characteristics
  if (method == "auto") {
    method <- select_task_method(data, task_data)
    message(sprintf("Selected method: %s", method))
  }
  
  # Handle automatic anchor selection
  if (identical(anchors, "auto")) {
    message("Selecting anchors automatically...")
    anchors <- select_anchors_auto(data, components, config)
  }
  
  # Configure based on method
  params <- configure_task_params(method, config, ...)
  
  # Call the appropriate function
  do.call(run_task_hatsa, c(
    list(
      subject_data_list = data,
      subject_task_activations_list = task_data,
      anchor_indices = anchors,
      spectral_rank_k = components
    ),
    params
  ))
}

#' HATSA Configuration Presets
#' @keywords internal
hatsa_preset <- function(name) {
  presets <- list(
    default = list(
      k_conn_pos = 10,
      k_conn_neg = 10,
      n_refine = 3,
      alpha_lrw = 0.93
    ),
    fast = list(
      k_conn_pos = 5,
      k_conn_neg = 5,
      n_refine = 1,
      alpha_lrw = 0.9
    ),
    accurate = list(
      k_conn_pos = 20,
      k_conn_neg = 20,
      n_refine = 5,
      alpha_lrw = 0.95
    )
  )
  presets[[name]]
}

#' Task-HATSA Configuration Presets
#' @keywords internal
hatsa_task_preset <- function(preset_name, method) {
  base <- hatsa_preset(preset_name)
  
  method_defaults <- list(
    blend = list(
      lambda_blend = 0.15,
      auto_residualize = TRUE
    ),
    gev = list(
      lambda_max_thresh = 0.8,
      min_stable_conditions = 2
    ),
    augmented = list(
      augmentation_weight = 0.5,
      normalize_augmentation = TRUE
    )
  )
  
  c(base, method_defaults[[method]])
}

#' Automatic Anchor Selection Helper
#' @keywords internal
select_anchors_auto <- function(data, k, config) {
  # Use existing select_anchors_mra with sensible defaults
  V_p <- ncol(data[[1]])
  
  # Reasonable defaults for anchor selection
  n_anchors <- max(
    round(0.1 * V_p),  # 10% of voxels
    k * 2,             # At least 2x components
    50                 # Minimum 50 anchors
  )
  n_anchors <- min(n_anchors, V_p, 200)  # Cap at 200
  
  select_anchors_mra(
    U_list = lapply(data, function(X) {
      # Simple SVD for each subject
      svd(scale(X), nu = k, nv = 0)$u
    }),
    k = k,
    n_anchors = n_anchors,
    n_sample = min(1000, V_p)
  )
}

#' Automatic Task Method Selection
#' @keywords internal  
select_task_method <- function(data, task_data) {
  # Simple heuristics for method selection
  n_timepoints <- nrow(data[[1]])
  n_conditions <- ncol(task_data[[1]])
  
  if (n_conditions > 10) {
    "gev"  # Many conditions benefit from GEV
  } else if (n_timepoints < 200) {
    "augmented"  # Short runs benefit from augmentation
  } else {
    "blend"  # Default to blend for typical data
  }
}

#' Configure Task Parameters
#' @keywords internal
configure_task_params <- function(method, base_config, ...) {
  user_params <- list(...)
  
  # Set method-specific parameters
  if (method == "blend") {
    base_config$lambda_blend <- user_params$lambda_blend %||% base_config$lambda_blend
    base_config$auto_residualize <- user_params$auto_residualize %||% TRUE
  } else if (method == "gev") {
    base_config$lambda_blend <- 0  # No blending for pure GEV
    base_config$lambda_max_thresh <- user_params$lambda_max_thresh %||% 0.8
  } else if (method == "augmented") {
    base_config$lambda_blend <- 0
    base_config$augmentation_weight <- user_params$augmentation_weight %||% 0.5
  }
  
  utils::modifyList(base_config, user_params)
}

#' Null-safe or operator
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}