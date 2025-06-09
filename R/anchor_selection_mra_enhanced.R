#' Enhanced MRA Anchor Selection with Parcel Quality Integration
#'
#' @description
#' Extensions to the MRA anchor selection algorithm that incorporate
#' parcel quality information into the selection process.
#'
#' @name anchor-selection-enhanced
NULL

#' Prepare Parcel Quality Information
#'
#' @description
#' Creates a parcel quality information structure that can be used
#' to guide anchor selection.
#'
#' @param U_pilot_list List of pilot subject sketch matrices
#' @param spectral_rank_k Number of spectral components
#' @param compute_metrics Character vector of metrics to compute
#' @param custom_quality Optional pre-computed quality scores
#'
#' @return A data frame with parcel quality information
#' @export
prepare_parcel_quality_info <- function(U_pilot_list,
                                        spectral_rank_k,
                                        compute_metrics = c("stability", "coverage", "snr"),
                                        custom_quality = NULL) {
  
  n_parcels <- nrow(U_pilot_list[[1]])
  quality_df <- data.frame(parcel_idx = 1:n_parcels)
  
  # Stability: How consistent is each parcel across subjects
  if ("stability" %in% compute_metrics) {
    stability_scores <- numeric(n_parcels)
    
    for (p in 1:n_parcels) {
      # Extract parcel p across all subjects
      parcel_rows <- lapply(U_pilot_list, function(U) {
        if (!is.null(U) && nrow(U) >= p) U[p, ] else NULL
      })
      parcel_rows <- Filter(Negate(is.null), parcel_rows)
      
      if (length(parcel_rows) > 1) {
        # Compute pairwise correlations
        parcel_mat <- do.call(rbind, parcel_rows)
        if (nrow(parcel_mat) > 1) {
          cor_mat <- cor(t(parcel_mat))
          # Average off-diagonal correlation as stability
          stability_scores[p] <- mean(cor_mat[upper.tri(cor_mat)])
        }
      }
    }
    
    quality_df$stability <- stability_scores
  }
  
  # Coverage: How much variance each parcel captures
  if ("coverage" %in% compute_metrics) {
    coverage_scores <- numeric(n_parcels)
    
    for (p in 1:n_parcels) {
      parcel_vars <- sapply(U_pilot_list, function(U) {
        if (!is.null(U) && nrow(U) >= p) {
          sum(U[p, ]^2) / spectral_rank_k
        } else NA
      })
      coverage_scores[p] <- mean(parcel_vars, na.rm = TRUE)
    }
    
    quality_df$coverage <- coverage_scores
  }
  
  # Signal-to-noise ratio estimate
  if ("snr" %in% compute_metrics) {
    snr_scores <- numeric(n_parcels)
    
    for (p in 1:n_parcels) {
      parcel_rows <- lapply(U_pilot_list, function(U) {
        if (!is.null(U) && nrow(U) >= p) U[p, ] else NULL
      })
      parcel_rows <- Filter(Negate(is.null), parcel_rows)
      
      if (length(parcel_rows) > 1) {
        parcel_mat <- do.call(rbind, parcel_rows)
        # Signal: variance of mean
        signal <- var(colMeans(parcel_mat))
        # Noise: mean within-subject variance
        noise <- mean(apply(parcel_mat, 1, var))
        snr_scores[p] <- if (noise > 0) signal / noise else 0
      }
    }
    
    quality_df$snr <- snr_scores
  }
  
  # Add custom quality if provided
  if (!is.null(custom_quality)) {
    if (length(custom_quality) == n_parcels) {
      quality_df$custom_quality <- custom_quality
    } else {
      warning("custom_quality length doesn't match number of parcels")
    }
  }
  
  # Compute composite quality score
  quality_metrics <- setdiff(names(quality_df), "parcel_idx")
  if (length(quality_metrics) > 0) {
    # Normalize each metric to 0-1
    for (metric in quality_metrics) {
      vals <- quality_df[[metric]]
      if (any(is.finite(vals))) {
        quality_df[[paste0(metric, "_norm")]] <- 
          (vals - min(vals, na.rm = TRUE)) / 
          (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
      }
    }
    
    # Composite score (average of normalized metrics)
    norm_cols <- grep("_norm$", names(quality_df), value = TRUE)
    if (length(norm_cols) > 0) {
      quality_df$composite_quality <- rowMeans(
        quality_df[norm_cols], 
        na.rm = TRUE
      )
    }
  }
  
  quality_df
}

#' Enhanced MRA Selection with Quality Filtering
#'
#' @description
#' Wrapper around select_anchors_mra that uses parcel quality information
#' to pre-filter candidates and weight selection.
#'
#' @param U_original_list_pilot Pilot subject sketch matrices
#' @param spectral_rank_k Number of spectral components
#' @param m_target Target number of anchors
#' @param parcel_quality_info Data frame with parcel quality information
#' @param quality_threshold Minimum quality score for candidate parcels (0-1)
#' @param quality_weight Weight for quality in selection (0-1)
#' @param ... Additional arguments passed to select_anchors_mra
#'
#' @return Selected anchor indices
#' @export
select_anchors_mra_quality <- function(U_original_list_pilot,
                                       spectral_rank_k,
                                       m_target,
                                       parcel_quality_info = NULL,
                                       quality_threshold = 0.3,
                                       quality_weight = 0.2,
                                       ...) {
  
  # Prepare quality info if not provided
  if (is.null(parcel_quality_info)) {
    message("Computing parcel quality metrics...")
    parcel_quality_info <- prepare_parcel_quality_info(
      U_pilot_list = U_original_list_pilot,
      spectral_rank_k = spectral_rank_k
    )
  }
  
  # Validate quality info
  if (!is.data.frame(parcel_quality_info) || 
      !"parcel_idx" %in% names(parcel_quality_info)) {
    warning("Invalid parcel_quality_info format. Proceeding without quality filtering.")
    return(select_anchors_mra(
      U_original_list_pilot = U_original_list_pilot,
      spectral_rank_k = spectral_rank_k,
      m_target = m_target,
      ...
    ))
  }
  
  # Determine quality metric to use
  quality_col <- if ("composite_quality" %in% names(parcel_quality_info)) {
    "composite_quality"
  } else if ("stability" %in% names(parcel_quality_info)) {
    "stability"
  } else {
    warning("No recognized quality metric in parcel_quality_info")
    NULL
  }
  
  # Filter candidate pool based on quality
  candidate_pool <- NULL
  if (!is.null(quality_col)) {
    quality_scores <- parcel_quality_info[[quality_col]]
    high_quality_parcels <- which(!is.na(quality_scores) & 
                                  quality_scores >= quality_threshold)
    
    message(sprintf("Found %d parcels above quality threshold %.2f",
                    length(high_quality_parcels), quality_threshold))
    
    if (length(high_quality_parcels) >= m_target) {
      candidate_pool <- high_quality_parcels
    } else {
      warning(sprintf(
        "Only %d parcels meet quality threshold, but need %d anchors. Using all parcels.",
        length(high_quality_parcels), m_target
      ))
    }
  }
  
  # TODO: Implement quality-weighted scoring within select_anchors_mra
  # This would require modifying the core function to accept and use
  # parcel quality scores in the .evaluate_anchor_set function
  
  # For now, use quality filtering via candidate_pool
  select_anchors_mra(
    U_original_list_pilot = U_original_list_pilot,
    spectral_rank_k = spectral_rank_k,
    m_target = m_target,
    candidate_pool = candidate_pool,
    parcel_quality_info = parcel_quality_info,  # Pass through for future use
    ...
  )
}

#' Network-Informed Anchor Selection
#'
#' @description
#' Select anchors that represent different functional networks or
#' brain regions.
#'
#' @param U_original_list_pilot Pilot subject sketch matrices
#' @param spectral_rank_k Number of spectral components
#' @param m_target Target number of anchors
#' @param network_labels Vector of network labels for each parcel
#' @param anchors_per_network Minimum anchors per network
#' @param ... Additional arguments for select_anchors_mra
#'
#' @return Selected anchor indices
#' @export
select_anchors_network_balanced <- function(U_original_list_pilot,
                                            spectral_rank_k,
                                            m_target,
                                            network_labels,
                                            anchors_per_network = 2,
                                            ...) {
  
  n_parcels <- nrow(U_original_list_pilot[[1]])
  
  if (length(network_labels) != n_parcels) {
    stop("network_labels must have one entry per parcel")
  }
  
  unique_networks <- unique(network_labels)
  n_networks <- length(unique_networks)
  
  # First, select minimum anchors from each network
  initial_selection <- integer(0)
  
  for (network in unique_networks) {
    network_parcels <- which(network_labels == network)
    
    if (length(network_parcels) >= anchors_per_network) {
      # Select top parcels from this network based on quality
      quality_info <- prepare_parcel_quality_info(
        U_pilot_list = U_original_list_pilot,
        spectral_rank_k = spectral_rank_k,
        compute_metrics = "coverage"
      )
      
      network_quality <- quality_info$coverage[network_parcels]
      top_idx <- order(network_quality, decreasing = TRUE)[1:anchors_per_network]
      initial_selection <- c(initial_selection, network_parcels[top_idx])
    }
  }
  
  message(sprintf("Selected %d initial anchors across %d networks",
                  length(initial_selection), n_networks))
  
  # Use MRA to select remaining anchors
  select_anchors_mra(
    U_original_list_pilot = U_original_list_pilot,
    spectral_rank_k = spectral_rank_k,
    m_target = m_target,
    initial_selection = initial_selection,
    ...
  )
}