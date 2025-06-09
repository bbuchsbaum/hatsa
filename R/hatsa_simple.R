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


