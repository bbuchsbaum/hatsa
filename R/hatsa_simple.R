#' Simple HATSA Interface
#'
#' @description
#' A streamlined interface to HATSA with sensible defaults and automatic parameter selection.
#'
#' @param data List of subject data matrices (time × voxels)
#' @param anchors Either "auto" for automatic selection, or a vector of anchor indices
#' @param components Number of spectral components (default: 20)
#' @param preset Configuration preset: "default", "fast", or "accurate"
#' @param ... Additional parameters passed to the internal HATSA engine
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
  # For now, use task_hatsa with core_hatsa method for basic HATSA
  result <- .task_hatsa_engine(
    subject_data_list = data,
    anchor_indices = anchors,
    spectral_rank_k = components,
    task_method = "core_hatsa",
    task_data_list = NULL,
    k_conn_pos = params$k_conn_pos,
    k_conn_neg = params$k_conn_neg,
    n_refine = params$n_refine,
    alpha_laplacian = params$alpha_lrw %||% 0.93,
    verbose = params$verbose %||% TRUE
  )
  
  return(result)
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
#' @param ... Additional parameters passed to the internal task HATSA engine
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
  
  # Map method names to task_method values
  task_method_map <- list(
    blend = "lambda_blend",
    gev = "gev_patch",
    augmented = "lambda_blend"  # augmented uses lambda_blend with row_augmentation
  )
  
  # Set task_method
  params$task_method <- task_method_map[[method]]
  
  # Set augmentation flag for augmented method
  if (method == "augmented") {
    params$row_augmentation <- TRUE
    params$lambda_blend_value <- 0  # No blending for augmented
  }
  
  # Call the internal engine directly
  do.call(.task_hatsa_engine, c(
    list(
      subject_data_list = data,
      task_data_list = task_data,
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
      n_refine = 3
    ),
    fast = list(
      k_conn_pos = 5,
      k_conn_neg = 5,
      n_refine = 1
    ),
    accurate = list(
      k_conn_pos = 20,
      k_conn_neg = 20,
      n_refine = 5
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
      lambda_blend_value = 0.15,
      check_redundancy = TRUE
    ),
    gev = list(
      gev_lambda_max = 0.8,
      k_gev_dims = 10
    ),
    augmented = list(
      row_augmentation = TRUE,
      lambda_blend_value = 0
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
  n_conditions <- nrow(task_data[[1]])  # Task data is C x V, so rows are conditions
  
  if (n_conditions > 10) {
    "gev"  # Many conditions benefit from GEV
  } else if (n_conditions < 3) {
    "augmented"  # Few conditions benefit from augmentation
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
    base_config$lambda_blend_value <- user_params$lambda_blend_value %||% base_config$lambda_blend_value
    base_config$check_redundancy <- user_params$check_redundancy %||% TRUE
  } else if (method == "gev") {
    base_config$lambda_blend_value <- 0  # No blending for pure GEV
    base_config$gev_lambda_max <- user_params$gev_lambda_max %||% 0.8
  } else if (method == "augmented") {
    base_config$lambda_blend_value <- 0
    # Note: augmentation_weight is not a parameter of .task_hatsa_engine
  }
  
  # Remove any parameters that were incorrectly named
  base_config$auto_residualize <- NULL
  base_config$lambda_blend <- NULL
  base_config$augmentation_weight <- NULL
  
  utils::modifyList(base_config, user_params)
}

