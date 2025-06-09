#' Anchor Selection Workflow for HATSA
#'
#' @description
#' Functions to support a two-phase workflow: 
#' 1) Anchor discovery/validation on training data
#' 2) Production use with validated anchors
#'
#' @name anchor-selection-workflow
NULL

#' Discover Optimal Anchors from Training Data
#'
#' @description
#' Comprehensive anchor discovery pipeline that evaluates multiple
#' selection strategies and quality metrics on a training dataset.
#'
#' @param training_data List of training subject data matrices
#' @param n_components Number of spectral components to use
#' @param n_anchors Target number of anchors
#' @param n_pilot Number of pilot subjects for initial selection (default: 10)
#' @param methods Character vector of selection methods to try
#' @param evaluate_stability Whether to evaluate anchor stability
#' @param save_results Whether to save results to file
#' @param output_dir Directory for saving results
#' @param verbose Print progress messages
#'
#' @return An anchor_selection_result object containing:
#'   - selected_anchors: Final selected anchor indices
#'   - quality_metrics: Comprehensive quality assessment
#'   - selection_report: Detailed report of selection process
#'
#' @export
discover_optimal_anchors <- function(training_data,
                                     n_components = 20,
                                     n_anchors = 100,
                                     n_pilot = min(10, length(training_data)),
                                     methods = c("mra", "coverage", "stability"),
                                     evaluate_stability = TRUE,
                                     save_results = TRUE,
                                     output_dir = "hatsa_anchors",
                                     verbose = TRUE) {
  
  if (verbose) {
    message("=== HATSA Anchor Discovery Pipeline ===")
    message(sprintf("Training subjects: %d", length(training_data)))
    message(sprintf("Target anchors: %d", n_anchors))
    message(sprintf("Methods to evaluate: %s", paste(methods, collapse = ", ")))
  }
  
  # Step 1: Initial decomposition on pilot subset
  pilot_idx <- sample(length(training_data), n_pilot)
  pilot_data <- training_data[pilot_idx]
  
  if (verbose) message("\n1. Computing pilot decompositions...")
  
  # Run initial HATSA on pilot data to get U matrices
  pilot_result <- hatsa(
    data = pilot_data,
    anchors = seq_len(min(n_anchors * 2, ncol(pilot_data[[1]]))),
    components = n_components,
    preset = "fast"  # Quick initial run
  )
  
  U_pilot_list <- pilot_result$U_original_list
  
  # Step 2: Compute parcel quality metrics
  if (verbose) message("\n2. Computing parcel quality metrics...")
  
  quality_info <- prepare_parcel_quality_info(
    U_pilot_list = U_pilot_list,
    spectral_rank_k = n_components,
    compute_metrics = c("stability", "coverage", "snr")
  )
  
  # Step 3: Try different selection methods
  if (verbose) message("\n3. Evaluating selection methods...")
  
  method_results <- list()
  
  if ("mra" %in% methods) {
    if (verbose) message("  - MRA selection...")
    mra_anchors <- select_anchors_mra(
      U_original_list_pilot = U_pilot_list,
      spectral_rank_k = n_components,
      m_target = n_anchors,
      parcel_quality_info = quality_info,
      verbose = FALSE
    )
    method_results$mra <- mra_anchors
  }
  
  if ("coverage" %in% methods) {
    if (verbose) message("  - Coverage-based selection...")
    coverage_anchors <- quality_info$parcel_idx[
      order(quality_info$coverage, decreasing = TRUE)[1:n_anchors]
    ]
    method_results$coverage <- coverage_anchors
  }
  
  if ("stability" %in% methods) {
    if (verbose) message("  - Stability-based selection...")
    stability_anchors <- quality_info$parcel_idx[
      order(quality_info$stability, decreasing = TRUE)[1:n_anchors]
    ]
    method_results$stability <- stability_anchors
  }
  
  # Step 4: Evaluate each method on full training data
  if (evaluate_stability && verbose) {
    message("\n4. Evaluating anchor sets on full training data...")
  }
  
  evaluation_results <- list()
  
  for (method_name in names(method_results)) {
    if (verbose) message(sprintf("  - Evaluating %s anchors...", method_name))
    
    anchors <- method_results[[method_name]]
    
    # Run HATSA with these anchors
    result <- hatsa(
      data = training_data,
      anchors = anchors,
      components = n_components,
      preset = "default"
    )
    
    # Compute quality metrics
    summary_metrics <- summary(result)
    recon_error <- reconstruction_error(result, type = "non_anchors")
    
    evaluation_results[[method_name]] <- list(
      anchors = anchors,
      mean_anchor_error = summary_metrics$mean_anchor_alignment_error,
      rotation_dispersion = summary_metrics$rotation_dispersion_mean_geo_dist,
      condition_number = summary_metrics$condition_number_T_anchor_final,
      non_anchor_reconstruction = mean(recon_error$per_subject_error),
      quality_scores = quality_info[quality_info$parcel_idx %in% anchors, ]
    )
  }
  
  # Step 5: Select best method
  if (verbose) message("\n5. Selecting best anchor set...")
  
  # Score each method (lower is better)
  method_scores <- sapply(evaluation_results, function(res) {
    score <- 0
    score <- score + res$mean_anchor_error * 1.0
    score <- score + res$rotation_dispersion * 0.5
    score <- score + log10(res$condition_number) * 0.1
    score <- score + res$non_anchor_reconstruction * 0.5
    score
  })
  
  best_method <- names(which.min(method_scores))
  selected_anchors <- method_results[[best_method]]
  
  if (verbose) {
    message(sprintf("\nBest method: %s", best_method))
    message(sprintf("Score: %.4f", method_scores[best_method]))
  }
  
  # Step 6: Create comprehensive report
  report <- list(
    selected_anchors = selected_anchors,
    selection_method = best_method,
    n_components = n_components,
    n_training_subjects = length(training_data),
    n_pilot_subjects = n_pilot,
    timestamp = Sys.time(),
    quality_metrics = evaluation_results[[best_method]],
    all_methods_comparison = data.frame(
      method = names(method_scores),
      score = method_scores,
      mean_anchor_error = sapply(evaluation_results, "[[", "mean_anchor_error"),
      rotation_dispersion = sapply(evaluation_results, "[[", "rotation_dispersion"),
      condition_number = sapply(evaluation_results, "[[", "condition_number")
    ),
    parcel_quality_info = quality_info
  )
  
  class(report) <- c("anchor_selection_result", "list")
  
  # Step 7: Save results
  if (save_results) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Save R object
    saveRDS(report, file.path(output_dir, "anchor_selection_result.rds"))
    
    # Save anchor indices as text
    writeLines(
      as.character(selected_anchors),
      file.path(output_dir, "selected_anchors.txt")
    )
    
    # Save summary report
    sink(file.path(output_dir, "selection_report.txt"))
    print(report)
    sink()
    
    if (verbose) {
      message(sprintf("\nResults saved to: %s", output_dir))
    }
  }
  
  report
}

#' Validate Pre-selected Anchors
#'
#' @description
#' Validates a set of pre-selected anchors on new data to ensure
#' they maintain good performance.
#'
#' @param data New dataset to validate on
#' @param anchor_indices Pre-selected anchor indices
#' @param n_components Number of components (should match selection)
#' @param reference_metrics Optional reference metrics from training
#'
#' @return Validation report
#' @export
validate_anchors <- function(data,
                             anchor_indices,
                             n_components = 20,
                             reference_metrics = NULL) {
  
  # Run HATSA with the anchors
  result <- hatsa(
    data = data,
    anchors = anchor_indices,
    components = n_components
  )
  
  # Get metrics
  summary_metrics <- summary(result)
  recon_error <- reconstruction_error(result, type = "non_anchors")
  
  validation <- list(
    n_subjects = length(data),
    n_anchors = length(anchor_indices),
    mean_anchor_error = summary_metrics$mean_anchor_alignment_error,
    rotation_dispersion = summary_metrics$rotation_dispersion_mean_geo_dist,
    condition_number = summary_metrics$condition_number_T_anchor_final,
    non_anchor_reconstruction = mean(recon_error$per_subject_error),
    timestamp = Sys.time()
  )
  
  # Compare to reference if provided
  if (!is.null(reference_metrics)) {
    validation$comparison <- data.frame(
      metric = c("mean_anchor_error", "rotation_dispersion", "condition_number"),
      training = c(
        reference_metrics$mean_anchor_error,
        reference_metrics$rotation_dispersion,
        reference_metrics$condition_number
      ),
      validation = c(
        validation$mean_anchor_error,
        validation$rotation_dispersion,
        validation$condition_number
      )
    )
    validation$comparison$ratio <- 
      validation$comparison$validation / validation$comparison$training
  }
  
  class(validation) <- c("anchor_validation_result", "list")
  validation
}

#' Load Pre-selected Anchors
#'
#' @description
#' Convenience function to load previously selected anchors.
#'
#' @param path Path to saved anchors (directory or file)
#' @param type Type of file to load ("rds", "txt", or "auto")
#'
#' @return Anchor indices or full selection result
#' @export
load_anchors <- function(path, type = "auto") {
  
  if (dir.exists(path)) {
    # It's a directory, look for files
    rds_file <- file.path(path, "anchor_selection_result.rds")
    txt_file <- file.path(path, "selected_anchors.txt")
    
    if (type == "auto") {
      if (file.exists(rds_file)) {
        type <- "rds"
        path <- rds_file
      } else if (file.exists(txt_file)) {
        type <- "txt"
        path <- txt_file
      } else {
        stop("No anchor files found in directory")
      }
    }
  }
  
  if (type == "rds" || (type == "auto" && grepl("\\.rds$", path))) {
    result <- readRDS(path)
    if (inherits(result, "anchor_selection_result")) {
      return(result$selected_anchors)
    } else {
      return(result)
    }
  } else {
    # Read text file
    as.integer(readLines(path))
  }
}

#' Print Anchor Selection Result
#' @export
print.anchor_selection_result <- function(x, ...) {
  cat("HATSA Anchor Selection Result\n")
  cat("=============================\n")
  cat(sprintf("Selected %d anchors using %s method\n", 
              length(x$selected_anchors), x$selection_method))
  cat(sprintf("Training subjects: %d\n", x$n_training_subjects))
  cat(sprintf("Components: %d\n", x$n_components))
  cat(sprintf("Timestamp: %s\n\n", x$timestamp))
  
  cat("Quality Metrics:\n")
  cat(sprintf("  Mean anchor error: %.4f\n", x$quality_metrics$mean_anchor_error))
  cat(sprintf("  Rotation dispersion: %.4f\n", x$quality_metrics$rotation_dispersion))
  cat(sprintf("  Condition number: %.1f\n", x$quality_metrics$condition_number))
  cat(sprintf("  Non-anchor reconstruction: %.4f\n", x$quality_metrics$non_anchor_reconstruction))
  
  if (!is.null(x$all_methods_comparison) && nrow(x$all_methods_comparison) > 1) {
    cat("\nMethod Comparison:\n")
    print(x$all_methods_comparison, row.names = FALSE)
  }
  
  invisible(x)
}

#' Print Anchor Validation Result  
#' @export
print.anchor_validation_result <- function(x, ...) {
  cat("HATSA Anchor Validation Result\n")
  cat("==============================\n")
  cat(sprintf("Subjects: %d\n", x$n_subjects))
  cat(sprintf("Anchors: %d\n", x$n_anchors))
  cat(sprintf("Timestamp: %s\n\n", x$timestamp))
  
  cat("Performance Metrics:\n")
  cat(sprintf("  Mean anchor error: %.4f\n", x$mean_anchor_error))
  cat(sprintf("  Rotation dispersion: %.4f\n", x$rotation_dispersion))
  cat(sprintf("  Condition number: %.1f\n", x$condition_number))
  cat(sprintf("  Non-anchor reconstruction: %.4f\n", x$non_anchor_reconstruction))
  
  if (!is.null(x$comparison)) {
    cat("\nComparison to Training:\n")
    print(x$comparison, row.names = FALSE)
    
    # Flag concerning changes
    if (any(x$comparison$ratio > 1.5)) {
      cat("\n[WARNING] Some metrics are >50% worse than training!\n")
    }
  }
  
  invisible(x)
}