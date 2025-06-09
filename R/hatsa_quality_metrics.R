#' Additional Quality Metrics for HATSA
#'
#' @description
#' Functions to extract and compute additional quality metrics from HATSA results,
#' implementing the remaining tickets from HMET series.
#'
#' @name hatsa-quality-metrics
NULL

#' Get Anchor Sketches
#'
#' @description
#' Extract subject-specific anchor sketch matrices (original or aligned).
#' Implements HMET-003.
#'
#' @param object A hatsa_projector object
#' @param type Either "original" or "aligned"
#' @param subject_idx Optional subject index (NULL for all subjects)
#' @param cache_in_object Whether to cache results in the object
#'
#' @return List of anchor sketch matrices or single matrix if subject_idx specified
#' @export
get_anchor_sketches <- function(object, 
                                type = c("original", "aligned"), 
                                subject_idx = NULL,
                                cache_in_object = FALSE) {
  
  if (!inherits(object, "hatsa_projector")) {
    stop("object must be a hatsa_projector")
  }
  
  type <- match.arg(type)
  anchor_indices <- object$parameters$anchor_indices
  
  if (length(anchor_indices) == 0) {
    warning("No anchor indices found")
    return(if (is.null(subject_idx)) list() else NULL)
  }
  
  # Check cache first
  cache_key <- paste0("A_", type, "_list_anchors")
  if (!is.null(object$._cache[[cache_key]])) {
    cached_result <- object$._cache[[cache_key]]
    if (!is.null(subject_idx)) {
      return(cached_result[[subject_idx]])
    }
    return(cached_result)
  }
  
  # Compute anchor sketches
  if (type == "original") {
    U_list <- object$U_original_list
  } else {
    # Reconstruct aligned list from s matrix
    U_list <- vector("list", object$parameters$N_subjects)
    for (i in seq_len(object$parameters$N_subjects)) {
      if (length(object$block_indices) >= i && !is.null(object$block_indices[[i]])) {
        rows_i <- object$block_indices[[i]]
        U_list[[i]] <- object$s[rows_i, , drop = FALSE]
      }
    }
  }
  
  # Extract anchor rows
  A_list <- lapply(U_list, function(U) {
    if (!is.null(U) && is.matrix(U) && nrow(U) >= max(anchor_indices)) {
      U[anchor_indices, , drop = FALSE]
    } else {
      NULL
    }
  })
  
  # Cache if requested
  if (cache_in_object) {
    object$._cache[[cache_key]] <- A_list
  }
  
  # Return specific subject or all
  if (!is.null(subject_idx)) {
    return(A_list[[subject_idx]])
  }
  A_list
}

#' Compute Parcel-Level Quality Metrics
#'
#' @description
#' Computes quality metrics for each parcel/voxel, including how well
#' each is represented in the aligned space.
#'
#' @param object A hatsa_projector object
#' @param type Type of quality metric: "reconstruction", "variance_captured", or "connectivity_preservation"
#'
#' @return A matrix (parcels × subjects) or vector of quality scores
#' @export
parcel_quality_metrics <- function(object, 
                                   type = c("reconstruction", 
                                            "variance_captured",
                                            "connectivity_preservation")) {
  
  if (!inherits(object, "hatsa_projector")) {
    stop("object must be a hatsa_projector")
  }
  
  type <- match.arg(type)
  
  N_subjects <- object$parameters$N_subjects
  V_p <- object$parameters$V_p
  k <- object$parameters$k
  
  if (type == "reconstruction") {
    # Compute reconstruction error per parcel
    quality_matrix <- matrix(NA, nrow = V_p, ncol = N_subjects)
    
    for (i in seq_len(N_subjects)) {
      U_orig <- object$U_original_list[[i]]
      R_i <- object$R_final_list[[i]]
      
      if (!is.null(U_orig) && !is.null(R_i)) {
        # Reconstruct: U_recon = U_aligned * R^T = (U_orig * R) * R^T
        U_recon <- (U_orig %*% R_i) %*% t(R_i)
        
        # Per-parcel error (row-wise norm)
        for (p in seq_len(V_p)) {
          quality_matrix[p, i] <- norm(U_orig[p, ] - U_recon[p, ], type = "2")
        }
      }
    }
    
    return(quality_matrix)
    
  } else if (type == "variance_captured") {
    # Compute how much variance each parcel contributes to the k components
    quality_matrix <- matrix(NA, nrow = V_p, ncol = N_subjects)
    
    for (i in seq_len(N_subjects)) {
      U_orig <- object$U_original_list[[i]]
      Lambda_i <- object$Lambda_original_list[[i]]
      
      if (!is.null(U_orig) && !is.null(Lambda_i)) {
        # Weight eigenvectors by eigenvalues
        U_weighted <- U_orig %*% diag(sqrt(Lambda_i))
        
        # Per-parcel contribution (row-wise norm squared)
        for (p in seq_len(V_p)) {
          quality_matrix[p, i] <- sum(U_weighted[p, ]^2)
        }
      }
    }
    
    return(quality_matrix)
    
  } else {
    # connectivity_preservation - placeholder
    stop("Connectivity preservation metric not yet implemented")
  }
}

#' Get Condition Numbers
#'
#' @description
#' Compute condition numbers for various matrices in HATSA results.
#' Implements HMET-006 and related metrics.
#'
#' @param object A hatsa_projector object
#' @param what Which condition number(s) to compute
#'
#' @return Named list of condition numbers
#' @export
get_condition_numbers <- function(object, 
                                  what = c("T_anchor_final", 
                                           "anchor_sketches",
                                           "all")) {
  
  if (!inherits(object, "hatsa_projector")) {
    stop("object must be a hatsa_projector")
  }
  
  what <- match.arg(what, several.ok = TRUE)
  if ("all" %in% what) {
    what <- c("T_anchor_final", "anchor_sketches")
  }
  
  results <- list()
  
  # Condition number of final anchor template
  if ("T_anchor_final" %in% what) {
    T_anchor <- object$T_anchor_final
    if (!is.null(T_anchor) && is.matrix(T_anchor) && 
        nrow(T_anchor) > 0 && ncol(T_anchor) > 0) {
      results$kappa_T_anchor_final <- kappa(T_anchor, exact = TRUE)
    } else {
      results$kappa_T_anchor_final <- NA
    }
  }
  
  # Condition numbers of individual anchor sketches
  if ("anchor_sketches" %in% what) {
    A_orig_list <- get_anchor_sketches(object, type = "original")
    
    kappa_values <- sapply(A_orig_list, function(A) {
      if (!is.null(A) && is.matrix(A) && nrow(A) > 0 && ncol(A) > 0) {
        kappa(A, exact = TRUE)
      } else {
        NA
      }
    })
    
    results$kappa_anchor_sketches_per_subject <- kappa_values
    results$kappa_anchor_sketches_mean <- mean(kappa_values, na.rm = TRUE)
    results$kappa_anchor_sketches_max <- max(kappa_values, na.rm = TRUE)
  }
  
  results
}

#' Enhanced Summary for HATSA with All Quality Metrics
#'
#' @description
#' Provides a comprehensive summary including all HMET metrics.
#'
#' @param object A hatsa_projector object
#' @param include_condition_numbers Whether to compute condition numbers
#' @param ... Additional arguments
#'
#' @return Enhanced summary object
#' @export
summary_enhanced <- function(object, include_condition_numbers = TRUE, ...) {
  
  # Get base summary
  base_summary <- summary(object)
  
  # Add condition numbers if requested
  if (include_condition_numbers) {
    cond_nums <- get_condition_numbers(object, what = "all")
    base_summary$condition_numbers <- cond_nums
  }
  
  # Add parcel quality summary
  parcel_recon <- parcel_quality_metrics(object, type = "reconstruction")
  base_summary$mean_parcel_reconstruction_error <- rowMeans(parcel_recon, na.rm = TRUE)
  
  # Find worst reconstructed parcels
  mean_errors <- rowMeans(parcel_recon, na.rm = TRUE)
  worst_parcels <- order(mean_errors, decreasing = TRUE)[1:min(10, length(mean_errors))]
  base_summary$worst_reconstructed_parcels <- worst_parcels
  
  # Add anchor quality specifically
  anchor_indices <- object$parameters$anchor_indices
  if (length(anchor_indices) > 0) {
    base_summary$anchor_reconstruction_quality <- mean_errors[anchor_indices]
  }
  
  class(base_summary) <- c("summary_enhanced_hatsa", class(base_summary))
  base_summary
}

#' Print Enhanced Summary
#'
#' @param x Enhanced summary object
#' @param ... Additional arguments
#' @export
print.summary_enhanced_hatsa <- function(x, ...) {
  
  # Call base print method
  NextMethod()
  
  # Add condition number information
  if (!is.null(x$condition_numbers)) {
    cat("\nCondition Numbers:\n")
    cat(sprintf("  Template (T_anchor_final): %.1f", 
                x$condition_numbers$kappa_T_anchor_final))
    
    if (x$condition_numbers$kappa_T_anchor_final > 1e4) {
      cat(" [WARNING: High condition number!]")
    }
    cat("\n")
    
    if (!is.null(x$condition_numbers$kappa_anchor_sketches_mean)) {
      cat(sprintf("  Anchor sketches (mean): %.1f\n", 
                  x$condition_numbers$kappa_anchor_sketches_mean))
      cat(sprintf("  Anchor sketches (max): %.1f\n", 
                  x$condition_numbers$kappa_anchor_sketches_max))
    }
  }
  
  # Add worst parcels info
  if (!is.null(x$worst_reconstructed_parcels)) {
    cat("\nPoorly Reconstructed Parcels (top 5):\n")
    worst_5 <- head(x$worst_reconstructed_parcels, 5)
    worst_errors <- x$mean_parcel_reconstruction_error[worst_5]
    for (i in seq_along(worst_5)) {
      cat(sprintf("  Parcel %d: error = %.4f\n", worst_5[i], worst_errors[i]))
    }
  }
  
  invisible(x)
}