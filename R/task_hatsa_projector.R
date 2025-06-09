# R/task_hatsa_projector.R

#' Task-HATSA Projector Object
#'
#' An S3 object of class `task_hatsa_projector` that stores the results of a
#' Task-Informed Harmonized Tensors SVD Alignment (task_hatsa) analysis.
#' This object inherits from `hatsa_projector`.
#'
#' @inherit hatsa_projector
#' @field qc_metrics A list per subject containing QC metrics like `rho_redundancy`
#'   (correlation between W_conn and W_task_raw) and `was_residualized` (boolean flag).
#' @field anchor_augmentation_info A list containing information related to anchor augmentation,
#'   such as `m_parcel_rows`, `m_task_rows_effective`, `condition_labels`,
#'   `was_residualized` (for condition anchors), `omega_mode_used`, `omega_weights_params`,
#'   and `trace_scaled`.
#' @field gev_patch_data A list containing GEV-specific outputs if `task_method == "gev_patch"`.
#'   Includes `U_patch_list`, `Lambda_patch_list`, `R_patch_list`, and `diagnostics`.
#'   NULL otherwise.
#'
#' @return A `task_hatsa_projector` object.
#' @seealso `task_hatsa`, `hatsa_projector`
#' @name task_hatsa_projector
NULL

#' Constructor for task_hatsa_projector S3 class
#'
#' Creates a `task_hatsa_projector` object. It takes the full results list from
#' the main `task_hatsa` computation (as prepared by `construct_output` helper)
#' and structures it into an S3 object, inheriting from `hatsa_projector`.
#'
#' @param task_hatsa_results A list containing the outputs from the `task_hatsa` run.
#'   This list is expected to have specific named elements corresponding to various
#'   stages and parameters of the analysis (e.g., `parameters`, `R_final_list`,
#'   `U_original_list`, `Lambda_original_list`, `T_anchor_final`, `U_aligned_list`,
#'   `qc_metrics`, `anchor_augmentation_info`, `gev_patch_data`).
#'
#' @return An object of class `c("task_hatsa_projector", "hatsa_projector", "multiblock_biprojector", "projector", "list")`.
#' @export
task_hatsa_projector <- function(task_hatsa_results) {

  # Extract parameters needed for base hatsa_projector components
  params <- task_hatsa_results$parameters
  if (is.null(params)) stop("task_hatsa_results$parameters is missing.")

  k <- params$spectral_rank_k
  N_subjects <- params$N_subjects
  V_p <- params$V_p

  if (is.null(k) || is.null(N_subjects) || is.null(V_p)) {
    stop("Essential parameters (k, N_subjects, V_p) missing from task_hatsa_results$parameters.")
  }

  # ---- Construct components for the base hatsa_projector ----
  # v: Mean aligned sketch (V_p x k)
  if (!is.list(task_hatsa_results$U_aligned_list) || length(task_hatsa_results$U_aligned_list) == 0) {
    # Allow for cases where U_aligned_list might be NULL if N_subjects is too small for GPA, etc.
    # but if N_subjects > 0, it should ideally be a list (even of NULLs for failed subjects)
    if (N_subjects > 0 && is.null(task_hatsa_results$U_aligned_list)) {
        stop("task_hatsa_results$U_aligned_list must be a list of matrices or NULLs.")
    }
    # Handle cases like k=0 or no valid subjects gracefully for v and s construction
    # If U_aligned_list is completely NULL or empty, v and s will be based on that.
    valid_aligned_sketches <- Filter(Negate(is.null), task_hatsa_results$U_aligned_list)
    if (length(valid_aligned_sketches) > 0) {
        v_sum <- Reduce("+", valid_aligned_sketches)
        v <- v_sum / length(valid_aligned_sketches)
    } else {
        v <- matrix(0, nrow = V_p, ncol = k) # Default if no valid sketches
    }
  } else {
    valid_aligned_sketches <- Filter(Negate(is.null), task_hatsa_results$U_aligned_list)
    if (length(valid_aligned_sketches) > 0) {
        v_sum <- Reduce("+", valid_aligned_sketches)
        v <- v_sum / length(valid_aligned_sketches)
    } else {
        v <- matrix(0, nrow = V_p, ncol = k) 
    }
  }
  
  # s: Stacked aligned sketches ((N_subjects * V_p) x k) or appropriate if some subjects failed
  # Ensure that NULLs in U_aligned_list are handled: create NA matrices or skip
  # For simplicity, we'll rbind, which will fail if dimensions don't match or if types are mixed with NULL directly.
  # A more robust approach might be to pre-allocate and fill, or replace NULLs with NA matrices.
  # Let's assume U_aligned_list contains conformable matrices (or NA matrices for failed subjects if that's the convention)
  # For now, use only valid sketches for 's', block_indices will need to reflect this.
  if (length(valid_aligned_sketches) > 0) {
      s <- do.call(rbind, valid_aligned_sketches)
      # block_indices: needs to map to rows of 's'
      # This assumes valid_aligned_sketches are ordered by original subject index.
      # A more robust way would be to track original indices of valid sketches.
      num_rows_per_valid_sketch <- sapply(valid_aligned_sketches, nrow)
      valid_subject_counts <- rep(1, length(valid_aligned_sketches)) # placeholder for now
      
      # Create block_indices based on the actual subjects that contributed to 's'
      # This is tricky if U_aligned_list contains NULLs for some subjects. 
      # The original hatsa_projector assumes all N_subjects contribute V_p rows each.
      # For task_hatsa, some subjects might fail basis computation.
      # The `scores` matrix `s` should ideally represent all Vp*N_subjects rows, with NAs for failed ones.
      # Let's try to conform to that for multivarious compatibility.
      s_full <- matrix(NA, nrow = N_subjects * V_p, ncol = k)
      current_row <- 1
      for (subj_idx in 1:N_subjects) {
          if (!is.null(task_hatsa_results$U_aligned_list[[subj_idx]])) {
              s_full[current_row:(current_row + V_p - 1), ] <- task_hatsa_results$U_aligned_list[[subj_idx]]
          }
          current_row <- current_row + V_p
      }
      s <- s_full
      block_indices <- split(seq_len(N_subjects * V_p), rep(seq_len(N_subjects), each = V_p))

  } else {
      s <- matrix(NA, nrow = N_subjects * V_p, ncol = k) # Default if no valid sketches
      block_indices <- split(seq_len(N_subjects * V_p), rep(seq_len(N_subjects), each = V_p))
      if (N_subjects == 0) block_indices <- list() # Handle N_subjects = 0 case
  }

  # sdev: Component standard deviations (length k) - default to 1s
  sdev <- rep(1, k)

  # preproc: Standard for multivarious projector
  preproc_obj <- multivarious::prep(multivarious::pass())

  # --- Assemble the final object ----
  obj <- list(
    # multivarious components
    v = v,
    s = s,
    sdev = sdev,
    preproc = preproc_obj,
    block_indices = block_indices,
    
    # Core HATSA components (from task_hatsa_results)
    R_final_list = task_hatsa_results$R_final_list,
    U_original_list = task_hatsa_results$U_original_list,
    Lambda_original_list = task_hatsa_results$Lambda_original_list,
    T_anchor_final = task_hatsa_results$T_anchor_final,
    
    # Parameters and method
    parameters = params, # Already contains task-specific details
    method = params$method, # Should be "task_hatsa"
    
    # Task-HATSA specific slots
    qc_metrics = task_hatsa_results$qc_metrics,
    anchor_augmentation_info = task_hatsa_results$anchor_augmentation_info,
    gev_patch_data = task_hatsa_results$gev_patch_data,
    
    # Add U_aligned_list for compatibility with tests
    U_aligned_list = task_hatsa_results$U_aligned_list,
    
    # Add task_method-specific outputs 
    W_task_list = task_hatsa_results$W_task_list,
    L_task_list = task_hatsa_results$L_task_list,
    W_hybrid_list = task_hatsa_results$W_hybrid_list,
    L_hybrid_list = task_hatsa_results$L_hybrid_list,
    U_task_list = task_hatsa_results$U_task_list,
    Lambda_task_list = task_hatsa_results$Lambda_task_list
  )

  # --- Initialize cache and store Fréchet mean of rotations ---
  obj$._cache <- list()
  valid_Rs_for_mean <- Filter(function(x) is.matrix(x) && !is.null(x), task_hatsa_results$R_final_list)
  if (length(valid_Rs_for_mean) > 0) {
    obj$._cache$R_frechet_mean <- tryCatch(
      frechet_mean_so_fast(valid_Rs_for_mean, refine = TRUE),
      error = function(e) if (k > 0) diag(k) else matrix(0, 0, 0)
    )
  } else {
    obj$._cache$R_frechet_mean <- if (k > 0) diag(k) else matrix(0, 0, 0)
  }

  # Ensure method is correctly set if not in params (it should be)
  if (is.null(obj$method)) obj$method <- "task_hatsa"
  if (is.null(obj$parameters$method)) obj$parameters$method <- "task_hatsa"

  class(obj) <- c("task_hatsa_projector", "hatsa_projector", "multiblock_biprojector", "projector", "list")

  return(obj)
}

#' Print method for task_hatsa_projector objects
#'
#' @param x A `task_hatsa_projector` object.
#' @param ... Additional arguments passed to print.
#' @return Invisibly returns the input object \code{x}.
#' @export
print.task_hatsa_projector <- function(x, ...) {
  # Call the print method for the parent class first
  NextMethod(generic = "print", object = x) # or print.hatsa_projector(x, ...)
  
  # Add Task-HATSA specific information
  cat("\nTask-HATSA Specifics:\n")
  cat("  Task Method: ", x$parameters$task_method, "\n")
  
  if (x$parameters$task_method == "lambda_blend") {
    cat("  Lambda Blend Value: ", x$parameters$lambda_blend_value, "\n")
  }
  
  if (x$parameters$row_augmentation && !is.null(x$anchor_augmentation_info)) {
    cat("  Anchor Augmentation: Enabled (", x$anchor_augmentation_info$m_task_rows_effective, " task rows added)\n")
    cat("    Omega Mode: ", x$anchor_augmentation_info$omega_mode_used, "\n")
    if (x$anchor_augmentation_info$was_residualized) {
        cat("    Condition Anchors Residualized: TRUE\n")
    }
  } else {
    cat("  Anchor Augmentation: Disabled\n")
  }
  
  if (x$parameters$task_method == "gev_patch" && !is.null(x$gev_patch_data)) {
    # Check if any patches were actually computed and stored successfully
    num_subjects_with_patches <- sum(!sapply(x$gev_patch_data$U_patch_list, is.null))
    if (num_subjects_with_patches > 0) {
        # Try to get k_gev_dims from first valid patch, or from parameters
        k_gev <- if (!is.null(x$gev_patch_data$U_patch_list[[which(!sapply(x$gev_patch_data$U_patch_list, is.null))[1]]])) {
                     ncol(x$gev_patch_data$U_patch_list[[which(!sapply(x$gev_patch_data$U_patch_list, is.null))[1]]])
                 } else {
                     x$parameters$k_gev_dims # Fallback to requested
                 }
        cat("  GEV Patch Data: Present (target k_gev_dims: ", x$parameters$k_gev_dims, ", found ~", k_gev, " dims for ", num_subjects_with_patches, " subjects)\n")
    } else {
        cat("  GEV Patch Data: No valid patches computed/stored.\n")
    }
  }
  
  invisible(x)
}

#' Summary method for task_hatsa_projector objects
#'
#' @param object A `task_hatsa_projector` object.
#' @param ... Additional arguments (unused).
#' @return A list object of class `summary.task_hatsa_projector` containing summary statistics.
#' @export
summary.task_hatsa_projector <- function(object, ...) {
  # Call summary method for the parent class to get base summary_info
  # This requires hatsa_projector to have a summary method.
  # If hatsa_projector.R has summary.hatsa_projector, NextMethod() should work.
  # Let's assume summary.hatsa_projector exists and populates basic fields.
  summary_info <- NextMethod(generic = "summary", object = object)
  
  # Add Task-HATSA specific summary information
  summary_info$task_method_used <- object$parameters$task_method
  if (object$parameters$task_method == "lambda_blend") {
    summary_info$lambda_blend_value <- object$parameters$lambda_blend_value
  }

  # QC Metrics Summary
  if (!is.null(object$qc_metrics) && length(object$qc_metrics) > 0) {
    valid_qc_metrics <- Filter(Negate(is.null), object$qc_metrics)
    if (length(valid_qc_metrics) > 0) {
        all_rhos <- sapply(valid_qc_metrics, function(qc) if(!is.null(qc$rho_redundancy) && is.finite(qc$rho_redundancy)) qc$rho_redundancy else NA)
        num_residualized <- sum(sapply(valid_qc_metrics, function(qc) qc$was_residualized), na.rm = TRUE)
        summary_info$qc_avg_rho_redundancy <- mean(all_rhos, na.rm = TRUE)
        summary_info$qc_median_rho_redundancy <- median(all_rhos, na.rm = TRUE)
        summary_info$qc_num_subjects_w_task_residualized <- num_residualized
        summary_info$qc_prop_subjects_w_task_residualized <- num_residualized / length(valid_qc_metrics)
    }
  }
  
  # Anchor Augmentation Summary
  if (!is.null(object$anchor_augmentation_info)) {
    summary_info$anchor_aug_enabled <- object$parameters$row_augmentation
    summary_info$anchor_aug_task_rows_effective <- object$anchor_augmentation_info$m_task_rows_effective
    summary_info$anchor_aug_omega_mode <- object$anchor_augmentation_info$omega_mode_used
    summary_info$anchor_aug_cond_residualized <- object$anchor_augmentation_info$was_residualized
  }
  
  # GEV Patch Summary
  if (object$parameters$task_method == "gev_patch" && !is.null(object$gev_patch_data)) {
    num_subjects_with_patches <- sum(!sapply(object$gev_patch_data$U_patch_list, is.null))
    summary_info$gev_patches_computed_for_subjects <- num_subjects_with_patches
    if (num_subjects_with_patches > 0) {
        actual_gev_dims <- sapply(Filter(Negate(is.null), object$gev_patch_data$U_patch_list), ncol)
        summary_info$gev_dims_found_summary <- summary(actual_gev_dims)
        # Could also summarize from gev_patch_data$diagnostics if that contains per-subject aggregated info
    }
  }
  
  class(summary_info) <- c("summary.task_hatsa_projector", "summary.hatsa_projector", class(summary_info))
  # Ensure correct class order, summary.hatsa_projector might already be there from NextMethod()
  # A safer way if NextMethod() already adds its class: 
  # current_class <- class(summary_info)
  # class(summary_info) <- unique(c("summary.task_hatsa_projector", current_class))
  return(summary_info)
}

#' Print method for summary.task_hatsa_projector objects
#'
#' @param x A `summary.task_hatsa_projector` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.summary.task_hatsa_projector <- function(x, ...) {
  # Call print method for parent summary class
  NextMethod(generic = "print", object = x) 
  
  cat("\n--- Task-HATSA Specific Summary ---\n")
  cat("Task Method Used: ", x$task_method_used, "\n")
  if (x$task_method_used == "lambda_blend" && !is.null(x$lambda_blend_value)) {
    cat("  Lambda Blend Value: ", sprintf("%.2f", x$lambda_blend_value), "\n")
  }
  
  if (!is.null(x$qc_avg_rho_redundancy)) {
      cat("Redundancy (W_conn vs W_task_raw):\n")
      cat("  Avg. Spearman Rho: ", sprintf("%.3f", x$qc_avg_rho_redundancy), "\n")
      cat("  Median Spearman Rho: ", sprintf("%.3f", x$qc_median_rho_redundancy), "\n")
      cat("  W_task Residualized for: ", x$qc_num_subjects_w_task_residualized, " (", 
          sprintf("%.1f%%", x$qc_prop_subjects_w_task_residualized * 100), ") subjects\n")
  }
  
  if (!is.null(x$anchor_aug_enabled)) {
    cat("Anchor Augmentation: ", ifelse(x$anchor_aug_enabled, "Enabled", "Disabled"), "\n")
    if (x$anchor_aug_enabled) {
      cat("  Effective Task Rows Added: ", x$anchor_aug_task_rows_effective, "\n")
      cat("  Omega Mode: ", x$anchor_aug_omega_mode, "\n")
      cat("  Condition Anchors Residualized: ", x$anchor_aug_cond_residualized, "\n")
    }
  }
  
  if (x$task_method_used == "gev_patch" && !is.null(x$gev_patches_computed_for_subjects)) {
    cat("GEV Patch Data:\n")
    cat("  Patches computed for: ", x$gev_patches_computed_for_subjects, " subjects\n")
    if (x$gev_patches_computed_for_subjects > 0 && !is.null(x$gev_dims_found_summary)) {
        cat("  Dimensions found per GEV patch (summary):\n")
        print(x$gev_dims_found_summary)
    }
  }
  invisible(x)
} 