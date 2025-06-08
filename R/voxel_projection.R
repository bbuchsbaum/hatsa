#' Compute Nyström Voxel Basis (Phi_voxel)
#'
#' Internal helper function to compute the Nyström extension basis for projecting
#' voxel-level data into a parcel-defined spectral space. This function uses
#' a k-nearest neighbors approach and a Gaussian kernel to create an affinity
#' matrix between voxels and parcels, then applies the Nyström formula.
#'
#' @param voxel_coords A numeric matrix (V_v x 3) of voxel coordinates.
#' @param parcel_coords A numeric matrix (V_p x 3) of parcel centroid coordinates.
#' @param U_orig_parcel A numeric matrix (V_p x k) of the subject's original
#'   parcel-level eigenvectors (spectral sketch).
#' @param Lambda_orig_parcel A numeric vector of length k, representing the subject's
#'   original parcel-level eigenvalues. Values should be positive.
#' @param n_nearest_parcels An integer, the number of nearest parcels (k_nn) to
#'   consider for each voxel when constructing the affinity matrix *if* `W_vox_parc`
#'   is not provided. Must be at least 1. Ignored if `W_vox_parc` is provided.
#' @param kernel_sigma A numeric scalar, the bandwidth (sigma) for the Gaussian kernel,
#'   or the string \"auto\". Used only *if* `W_vox_parc` is not provided.
#'   If \"auto\", sigma is estimated as
#'   `median(dist_to_1st_nn_parcel) / sqrt(2)`, where `dist_to_1st_nn_parcel` are
#'   the Euclidean distances from each voxel to its closest parcel centroid.
#'   A fallback value (e.g., 1.0) is used if auto-estimation is not possible.
#'   Defaults to 5.0. Ignored if `W_vox_parc` is provided.
#' @param row_normalize_W A logical. Controls whether the effective affinity matrix
#'   (`W_vox_parc`, either computed internally or provided) is row-normalized before
#'   computing Phi_voxel. The core HATSA algorithm uses an alpha-lazy random-walk
#'   normalized Laplacian (`L_rw_lazy`), and for consistency, the standard Nyström
#'   extension involves row-normalizing the voxel-parcel affinities.
#'   Defaults to `TRUE`.
#' @param eigenvalue_floor A small positive numeric value to floor near-zero eigenvalues
#'   before inversion. Defaults to 1e-8.
#' @param W_vox_parc Optional. A pre-computed sparse matrix (dgCMatrix, V_v x V_p)
#'   representing voxel-to-parcel affinities or weights. If provided, the function
#'   will skip the internal k-NN search and Gaussian kernel calculation and use
#'   this matrix directly. It must have dimensions V_v x V_p. Default is `NULL`,
#'   triggering internal calculation.
#' @param ... Additional arguments (currently unused).
#'
#' @return A dense numeric matrix \code{Phi_voxel} (V_v x k), representing the
#'   Nyström basis for voxels. Note: For large numbers of voxels (V_v) and/or
#'   components (k), this matrix can be memory-intensive (e.g., V_v=400k, k=50
#'   using 8-byte doubles is ~153MB).
#'
#' @importFrom Matrix Matrix Diagonal sparseMatrix t 
#' @keywords internal
#' @examples
#' # V_p <- 100; V_v <- 500; k <- 10; n_nearest <- 5; sigma <- 5
#' # p_coords <- matrix(rnorm(V_p*3), V_p, 3)
#' # v_coords <- matrix(rnorm(V_v*3), V_v, 3)
#' # U_p <- matrix(rnorm(V_p*k), V_p, k)
#' # L_p <- runif(k, 0.1, 1)
#' # if (requireNamespace("RANN", quietly = TRUE) && 
#' #     requireNamespace("Matrix", quietly = TRUE)) {
#' #   Phi_v <- compute_voxel_basis_nystrom(v_coords, p_coords, U_p, L_p, 
#' #                                        n_nearest, sigma)
#' #   # dim(Phi_v) # Should be V_v x k
#' # }
compute_voxel_basis_nystrom <- function(voxel_coords, parcel_coords, 
                                        U_orig_parcel, Lambda_orig_parcel,
                                        n_nearest_parcels = 10, kernel_sigma = 5.0,
                                        row_normalize_W = TRUE,
                                        eigenvalue_floor = 1e-8, 
                                        W_vox_parc = NULL, ...) {

  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("The 'RANN' package is required for compute_voxel_basis_nystrom. Please install it.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required. Please install it.") # Should be a dependency anyway
  }

  V_v <- nrow(voxel_coords)
  V_p <- nrow(parcel_coords)
  k <- ncol(U_orig_parcel)

  if (V_v == 0) return(matrix(0, nrow = 0, ncol = k))
  if (V_p == 0 || k == 0) return(matrix(0, nrow = V_v, ncol = k))
  if (nrow(U_orig_parcel) != V_p) stop("U_orig_parcel row count must match parcel_coords row count.")
  if (length(Lambda_orig_parcel) != k) stop("Length of Lambda_orig_parcel must match column count of U_orig_parcel (k).")
  if (any(Lambda_orig_parcel <= 0) && interactive()) {
      warning("Some Lambda_orig_parcel are <= 0. This might cause issues with inversion. Applying floor.")
  }

  if (is.null(W_vox_parc)) {
    # 1. Find k_nn nearest parcel centroids for each voxel
    if (n_nearest_parcels < 1) stop("`n_nearest_parcels` must be at least 1.")
    nn_results <- RANN::nn2(data = parcel_coords, query = voxel_coords, k = n_nearest_parcels, treetype = "kd")
    
    # nn_results$nn.idx gives V_v x n_nearest_parcels matrix of parcel indices
    # nn_results$nn.dists gives V_v x n_nearest_parcels matrix of SQUARED Euclidean distances

    # Determine effective kernel_sigma (Ticket V-R2)
    kernel_sigma_effective <- kernel_sigma
    if (is.character(kernel_sigma) && kernel_sigma == "auto") {
      if (V_v > 0 && ncol(nn_results$nn.dists) >= 1) {
        # Use distance to the 1st nearest neighbor for auto-tuning sigma
        # nn.dists contains SQUARED distances
        first_nn_sqrt_dists <- sqrt(nn_results$nn.dists[, 1])
        median_first_nn_dist <- median(first_nn_sqrt_dists, na.rm = TRUE)
        
        if (is.finite(median_first_nn_dist) && median_first_nn_dist > 1e-6) { # Avoid zero or tiny sigma
          kernel_sigma_effective <- median_first_nn_dist / sqrt(2)
          if (interactive()) {
            message(sprintf("Auto-tuned kernel_sigma to: %.3f (based on median 1st NN distance / sqrt(2))", kernel_sigma_effective))
          }
        } else {
          kernel_sigma_effective <- 1.0 # Fallback if median dist is zero/NA (e.g. all voxels on parcel centers)
          if (interactive()) {
            message(sprintf("Could not auto-tune kernel_sigma (median 1st NN distance was %.3f). Using fallback: %.3f", median_first_nn_dist, kernel_sigma_effective))
          }
        }
      } else {
        kernel_sigma_effective <- 1.0 # Fallback if no voxels or nn_results are unexpected
        if (interactive()) {
          message(sprintf("Could not auto-tune kernel_sigma (no voxel data or NN issue). Using fallback: %.3f", kernel_sigma_effective))
        }
      }
    } else if (!is.numeric(kernel_sigma) || length(kernel_sigma) != 1 || kernel_sigma <= 0) {
      stop("`kernel_sigma` must be a positive numeric value or the string \"auto\".")
    }

    # 2. Compute similarities using Gaussian kernel
    # W_ij = exp(-d_ij^2 / (2 * sigma_eff^2))
    similarities <- exp(-nn_results$nn.dists / (2 * kernel_sigma_effective^2))

    # affinities_flat: all affinity values, flattened
    affinities_flat <- as.vector(t(similarities))

    # Construct sparse affinity matrix W_vox_parc (V_v x V_p)
    # First, create in triplet form (TsparseMatrix). This allows for duplicate (i,j) entries.
    W_T <- Matrix::sparseMatrix(i = rep(1:V_v, each = n_nearest_parcels),
                              j = as.vector(t(nn_results$nn.idx)),
                              x = affinities_flat,
                              dims = c(V_v, V_p),
                              repr = "T", # Specify triplet representation
                              giveCsparse = FALSE) # Do not convert to Csparse yet

    # Convert to CsparseMatrix (dgCMatrix). This step automatically sums the 'x' values
    # for any duplicate (i,j) pairs present in the triplet form.
    # This handles cases where a voxel's k-NN includes the same parcel index multiple times,
    # or if multiple (voxel_idx, parcel_idx) raw pairs exist before aggregation.
    W_vox_parc <- as(W_T, "dgCMatrix")

  } else {
    message_stage("Using provided W_vox_parc matrix.", interactive_only=TRUE)
    if (!inherits(W_vox_parc, "dgCMatrix")) {
        warning("Provided `W_vox_parc` is not a dgCMatrix. Attempting to coerce.")
        W_vox_parc <- as(W_vox_parc, "dgCMatrix")
    }
    if (nrow(W_vox_parc) != V_v || ncol(W_vox_parc) != V_p) {
      stop(sprintf("Provided `W_vox_parc` has dimensions [%d x %d], expected [%d x %d] based on voxel/parcel counts.",
                   nrow(W_vox_parc), ncol(W_vox_parc), V_v, V_p))
    }
  }

  # 3. Optionally row-normalize W_vox_parc
  if (row_normalize_W) {
    # Row-normalize W_vox_parc so that the sum of each row is 1.
    # This is a common heuristic in Nyström extensions based on random walks.
    row_sums_W <- Matrix::rowSums(W_vox_parc)
    
    # Identify rows that are not all zero.
    non_zero_row_indices <- row_sums_W != 0
    
    inv_row_sums <- numeric(V_v) # Initialize with zeros
    if (any(non_zero_row_indices)) { # Proceed only if there are non-zero rows
        inv_row_sums[non_zero_row_indices] <- 1 / row_sums_W[non_zero_row_indices]
    }
    
    # Create a sparse diagonal matrix for efficient multiplication.
    # For rows that were originally all-zero (e.g., voxels with no affinity to any
    # selected parcels), their `inv_row_sums` entry is 0.
    # Multiplying by D_inv_sparse will therefore keep these all-zero rows as all-zero.
    # This means such voxels will have a zero projection, which is the chosen behavior.
    D_inv_sparse <- Matrix::Diagonal(n = V_v, x = inv_row_sums)
    W_vox_parc <- D_inv_sparse %*% W_vox_parc
  }

  # 4. Prepare inverse of eigenvalues (Lambda_orig_parcel_safe)
  Lambda_orig_parcel_safe <- pmax(Lambda_orig_parcel, eigenvalue_floor)
  Lambda_inv_diag <- Matrix::Diagonal(n = k, x = 1 / Lambda_orig_parcel_safe)

  # 5. Compute Phi_voxel = W_vox_parc %*% U_orig_parcel %*% Lambda_inv_diag
  # (V_v x V_p) %*% (V_p x k) %*% (k x k)  => (V_v x k)
  # Matrix::%*% handles sparse-dense multiplication efficiently
  Phi_voxel <- W_vox_parc %*% U_orig_parcel %*% Lambda_inv_diag
  
  if (!is.matrix(Phi_voxel)) { # e.g. if it became a Matrix object
      Phi_voxel <- as.matrix(Phi_voxel)
  }

  return(Phi_voxel)
}

#' Internal helper to validate coordinate system inputs
#'
#' Checks for gross inconsistencies between voxel and parcel coordinates, 
#' such as vastly different ranges or potential scale issues, by comparing
#' their bounding boxes.
#'
#' @param voxel_coords Matrix of voxel coordinates (V_v x 3).
#' @param parcel_coords Matrix of parcel coordinates (V_p x 3).
#' @param scale_threshold Factor by which coordinate spans can differ before warning.
#' @param absolute_range_threshold Warn if ranges don't overlap and are separated by more than this.
#' @keywords internal
.validate_coordinate_inputs <- function(voxel_coords, parcel_coords, 
                                        scale_threshold = 10, 
                                        absolute_range_threshold = 50) { # e.g. 50mm
  if (is.null(voxel_coords) || is.null(parcel_coords) || 
      !is.matrix(voxel_coords) || !is.matrix(parcel_coords) ||
      ncol(voxel_coords) != 3 || ncol(parcel_coords) != 3 ||
      nrow(voxel_coords) == 0 || nrow(parcel_coords) == 0) {
    # Basic checks already in project_voxels, but good for direct use too
    return(invisible(NULL)) 
  }

  dims <- c("X", "Y", "Z")
  warnings_found <- character(0)

  for (i in 1:3) {
    vox_range <- range(voxel_coords[, i], na.rm = TRUE)
    par_range <- range(parcel_coords[, i], na.rm = TRUE)

    vox_span <- diff(vox_range)
    par_span <- diff(par_range)

    # Check for non-overlapping ranges with significant separation
    # (e.g., max of one is far from min of other)
    ranges_disjoint_by_thresh <- (par_range[1] > vox_range[2] + absolute_range_threshold) || 
                                 (vox_range[1] > par_range[2] + absolute_range_threshold)
    
    if (ranges_disjoint_by_thresh) {
      warnings_found <- c(warnings_found, 
        sprintf("Dimension %s: Voxel range [%.1f, %.1f] and parcel range [%.1f, %.1f] are substantially separated. Possible coordinate system mismatch.", 
                dims[i], vox_range[1], vox_range[2], par_range[1], par_range[2]))
    }
    
    # Check for scale differences (e.g. mm vs cm)
    if (par_span > 0 && vox_span > 0) { # Avoid div by zero if one is a point cloud
      if (vox_span / par_span > scale_threshold || par_span / vox_span > scale_threshold) {
        warnings_found <- c(warnings_found,
          sprintf("Dimension %s: Voxel span (%.1f) and parcel span (%.1f) differ by more than %.0fx. Possible unit or scale mismatch.", 
                  dims[i], vox_span, par_span, scale_threshold))
      }
    }
  }
  
  # Check if parcel bounding box is roughly contained within voxel bounding box (allowing some leeway)
  # This checks if parcel min/max are "too far" outside voxel min/max respectively
  for (i in 1:3) {
      vox_min <- min(voxel_coords[,i], na.rm=TRUE); vox_max <- max(voxel_coords[,i], na.rm=TRUE)
      par_min <- min(parcel_coords[,i], na.rm=TRUE); par_max <- max(parcel_coords[,i], na.rm=TRUE)
      
      # Heuristic: if a parcel coord is outside the voxel range by more than, say, 10% of voxel span
      leeway <- diff(range(voxel_coords[,i], na.rm=TRUE)) * 0.1 
      
      if (par_min < vox_min - leeway || par_max > vox_max + leeway) {
          warnings_found <- c(warnings_found,
              sprintf("Dimension %s: Parcel coordinates range [%.1f, %.1f] extends notably beyond voxel coordinates range [%.1f, %.1f]. Ensure parcels are within the voxel volume.",
                      dims[i], par_min, par_max, vox_min, vox_max))
      }
  }

  if (length(warnings_found) > 0 && interactive()) {
    message("Coordinate system validation found potential issues (run non-interactively to suppress):")
    for(w in warnings_found) message(paste("  -", w))
  }
  return(invisible(NULL))
}

#' Project Voxel-Level Data
#'
#' Generic S3 method for projecting voxel-level data using a fitted model object.
#'
#' @param object A fitted model object (e.g., \code{hatsa_projector}).
#' @param voxel_timeseries_list A list of voxel time-series matrices.
#' @param voxel_coords Coordinates of the voxels.
#' @param parcel_coords Coordinates of the parcels used in the model.
#' @param ... Additional arguments specific to the method.
#'
#' @return A list of projected voxel data, specific to the method implementation.
#' @export
project_voxels <- function(object, voxel_timeseries_list, voxel_coords, parcel_coords, ...) {
  UseMethod("project_voxels")
}

#' Project Voxel-Level Data using a HATSA Projector
#'
#' Projects voxel-level time series data into the common aligned space defined
#' by a \code{hatsa_projector} object. This method uses Nyström extension.
#'
#' It is assumed that \code{voxel_coords} and \code{parcel_coords} are in the
#' same Cartesian coordinate system and units (e.g., millimeters in an MNI-aligned space).
#' The function includes heuristic checks for gross inconsistencies in these coordinate
#' systems (see \code{.\link{.validate_coordinate_inputs}} for details on checks performed),
#' issuing messages if potential issues like substantially different ranges or scales are detected.
#' The Nyström extension is formulated to be consistent with the unnormalized graph
#' Laplacian used in the core HATSA algorithm for parcel-level decomposition. See
#' \code{\link{compute_voxel_basis_nystrom}} for more details on parameters like
#' \code{row_normalize_W} that control the Nyström calculation.
#'
#' @section Methodological Details for Voxel Projection:
#' This section provides further details on the assumptions and computations
#' involved in the Nyström-based voxel projection implemented here.
#'
#' \subsection{Coordinate Systems and Validation:}
#' It is critically assumed that the \code{voxel_coords} (V_v x 3 matrix of voxel
#' coordinates) and \code{parcel_coords} (V_p x 3 matrix of parcel centroid
#' coordinates) are in the **same Cartesian coordinate system and units**.
#' Typically, this would be a standard neuroimaging space like MNI, with
#' coordinates in millimeters. The function \code{\link{.validate_coordinate_inputs}}
#' performs heuristic checks for gross inconsistencies (e.g., vastly different
#' data ranges or scales between voxel and parcel coordinates) and issues
#' messages if potential issues are detected. However, ensuring coordinate
#' system compatibility remains the user\'s responsibility.
#'
#' \subsection{Distance Metric for Voxel-Parcel Affinities:}
#' The affinities between voxels and parcels, used to construct the
#' \code{W_vox_parc} matrix in the Nyström extension (see
#' \code{\link{compute_voxel_basis_nystrom}}), are based on spatial proximity.
#' The k-nearest parcel neighbors for each voxel are found using
#' \code{RANN::nn2}, which, by default, computes **Euclidean distances**
#' in the provided coordinate space. These squared Euclidean distances are
#' then used in the Gaussian kernel.
#'
#' \subsection{Assumptions for Voxel Time-Series Data:}
#' The input \code{voxel_timeseries_list} should contain matrices of
#' pre-processed voxel BOLD fMRI time-series (T_i time points x V_v voxels).
#' While this function does not perform time-series pre-processing itself,
#' the quality and interpretation of the projected voxel coefficients
#' (\code{C_voxel_aligned_i}) will depend on the input data. It is generally
#' assumed that standard fMRI pre-processing steps such as motion correction,
#' slice-timing correction, spatial normalization (to the same space as
#' \code{voxel_coords} and \code{parcel_coords}), detrending, and potentially
#' temporal filtering and nuisance regression have been applied to the
#' voxel time-series prior to their use with this function. The projection
#' \code{C_voxel_coeffs_i = (1/T_i) * t(X_voxel_ts_i) %*% Phi_voxel_i}
#' effectively calculates the covariance of the voxel time series with the
#' Nyström-extended parcel-derived basis functions.
#'
#' \subsection{Kernel Sigma (\code{kernel_sigma}):}
#' The \code{kernel_sigma} parameter in \code{\link{compute_voxel_basis_nystrom}}
#' controls the bandwidth of the Gaussian kernel used to compute affinities
#' between voxels and their nearest parcel neighbors.
#' \itemize{
#'   \item If a numeric value is provided, it is used directly as sigma.
#'   \item If \code{kernel_sigma = "auto"} (the default in \code{project_voxels}),
#'     sigma is estimated as \code{median(dist_to_1st_nn_parcel) / sqrt(2)},
#'     where \code{dist_to_1st_nn_parcel} are the Euclidean distances from each
#'     voxel to its single closest parcel centroid. A fallback value (e.g., 1.0)
#'     is used if auto-estimation is problematic (e.g., all distances are zero).
#' }
#' The choice of sigma can influence the smoothness and reach of the
#' voxel-to-parcel affinities.
#'
#' \subsection{Nyström Formulation and Laplacian Consistency:}
#' The core HATSA algorithm, as implemented in \code{\link{run_hatsa_core}} and
#' its helper \code{compute_spectral_sketch_sparse} (using the updated
#' \code{compute_graph_laplacian_sparse}), derives the original
#' subject-specific parcel components (\code{U_original_list\\\[\\\[i\\\]\\\]} and
#' \code{Lambda_original_list\\\[\\\[i\\\]\\\]}) from an **alpha-lazy random-walk
#' normalized graph Laplacian** (\eqn{L_{rw\_lazy} = (I - \alpha D_p^{-1} W_p + (I - \alpha D_p^{-1} W_p)^T)/2},
#' where \eqn{W_p} is the parcel-parcel affinity matrix, \eqn{D_p} is the
#' corresponding degree matrix (typically based on absolute weights), and \eqn{\alpha}
#' is the laziness parameter, defaulting to 0.93). This symmetric form ensures compatibility
#' with symmetric eigensolvers.
#'
#' For mathematical consistency with this parcel-level eigensystem, the
#' Nyström extension formula used in \code{\link{compute_voxel_basis_nystrom}}
#' effectively involves row-normalizing the voxel-parcel affinity matrix \eqn{W_{vp}}
#' (making it row-stochastic, equivalent to pre-multiplying by \eqn{D_v^{-1}})
#' before applying the projection using the parcel eigenvectors and eigenvalues.
#' This formulation corresponds to setting the \code{row_normalize_W} parameter in
#' \code{\link{compute_voxel_basis_nystrom}} to \code{TRUE}, which is its default.
#' Using \code{row_normalize_W = FALSE} would generally be inconsistent with the
#' default HATSA Laplacian type.
#' (Note: The exact Nyström formula relating \eqn{L_{rw\_lazy}} eigenvectors/values
#' to the row-normalized extension requires careful derivation, but using the
#' eigenvectors of \eqn{L_{rw\_lazy}} with the row-normalized \eqn{W_{vp}} extension
#' is a common and practically effective approach consistent with random-walk principles).
#'
#' \subsection{Handling of DC Component in Parcel-Level Spectral Decomposition:}
#' The parcel-level eigenvectors \code{U_original_list\\\[\\\[i\\\]\\\]} (and their
#' corresponding eigenvalues \code{Lambda_original_list\\\[\\\[i\\\]\\\]}) used as input
#' to the Nyström extension are **non-DC components**.
#' In \code{compute_spectral_sketch_sparse}, eigenvectors of the parcel
#' Laplacian associated with eigenvalues near zero (i.e., numerically zero,
#' including the DC component for a connected graph) are explicitly filtered out
#' using a tolerance (\code{eigenvalue_tol}). The \code{spectral_rank_k}
#' parameter thus specifies the number of desired "informative" or non-DC
#' components from the parcel-level decomposition. This ensures that the
#' Nyström extension \code{Phi_voxel_i} is constructed based on these
#' non-trivial modes of parcel co-variation.
#'
#' @param object A fitted \code{hatsa_projector} object.
#' @param voxel_timeseries_list A list of voxel time-series matrices (T_i x V_v).
#'   It is assumed that the i-th element of this list corresponds to the i-th
#'   subject as stored in the \code{hatsa_projector} object (e.g., for U_original_list).
#' @param voxel_coords A numeric matrix (V_v x 3) of voxel coordinates.
#' @param parcel_coords A numeric matrix (V_p x 3) of parcel centroid coordinates
#'   corresponding to the parcellation used to fit the \code{object}.
#' @param n_nearest_parcels Integer, number of nearest parcels for Nyström extension.
#'   Passed to \code{\link{compute_voxel_basis_nystrom}}.
#' @param kernel_sigma Numeric or "auto", kernel bandwidth for Nyström extension.
#'   Passed to \code{\link{compute_voxel_basis_nystrom}}.
#' @param data_type Character string, either "timeseries" (default) or
#'   "coefficients". If "timeseries", the projection involves scaling by `1/T_i`
#'   (number of time points) to estimate covariance with the basis.
#'   If "coefficients" (e.g., for projecting beta maps or statistical maps),
#'   this scaling is omitted.
#' @param ... Additional arguments passed to \code{\link{compute_voxel_basis_nystrom}}.
#'
#' @return A list of aligned voxel coefficient matrices. Each element `[[i]]` is a
#'   matrix of dimensions T_i x k, where T_i is the number of time points (or rows)
#'   in the input `voxel_timeseries_list[[i]]`, and k is the number of HATSA components.
#'   These represent the projection coefficients of the voxel data onto the aligned
#'   HATSA basis for each subject.
#'
#' @examples
#' # This is a conceptual example. For it to run, you need a fitted hatsa_projector.
#' # First, let's set up parameters for a minimal run_hatsa_core call.
#' set.seed(456)
#' N_subj_fit <- 2
#' V_parc_fit <- 20 # Number of parcels in the fitted model
#' k_comp_fit <- 3   # Number of components in the fitted model
#' T_time_fit <- 40  # Number of time points for parcel data
#'
#' # Generate mock parcel-level data for fitting HATSA
#' subject_parcel_data <- lapply(1:N_subj_fit, function(i) {
#'   matrix(stats::rnorm(T_time_fit * V_parc_fit), ncol = V_parc_fit)
#' })
#' anchor_idx_fit <- sample(1:V_parc_fit, min(V_parc_fit, 5))
#'
#' # Fit a hatsa_projector object (requires Matrix and RSpectra)
#' fitted_hatsa_model <- NULL
#' if (requireNamespace("Matrix", quietly = TRUE) && 
#'     requireNamespace("RSpectra", quietly = TRUE) &&
#'     exists("run_hatsa_core")) {
#'   fitted_hatsa_model <- tryCatch({
#'     run_hatsa_core(
#'       subject_data_list = subject_parcel_data,
#'       anchor_indices    = anchor_idx_fit,
#'       spectral_rank_k = k_comp_fit,
#'       k_conn_pos        = min(5, V_parc_fit -1),
#'       k_conn_neg        = min(2, V_parc_fit -1),
#'       n_refine          = 1
#'     )
#'   }, error = function(e) NULL) # Return NULL on error for example
#' }
#'
#' if (!is.null(fitted_hatsa_model)) {
#'   # Now, prepare data for project_voxels
#'   N_subj_proj <- N_subj_fit # Number of subjects for voxel projection
#'   V_vox_proj  <- 50         # Number of voxels
#'   T_time_vox  <- 35         # Number of time points for voxel data
#'
#'   # Mock voxel time-series data
#'   voxel_ts_list <- lapply(1:N_subj_proj, function(i) {
#'     matrix(stats::rnorm(T_time_vox * V_vox_proj), ncol = V_vox_proj)
#'   })
#'
#'   # Mock coordinates
#'   voxel_coords_map <- matrix(stats::rnorm(V_vox_proj * 3), ncol = 3)
#'   # Parcel coords must match the V_p used when fitting the model
#'   parcel_coords_map <- matrix(stats::rnorm(V_parc_fit * 3), ncol = 3) 
#'
#'   # Project voxel data
#'   projected_vox_coeffs <- project_voxels(
#'     object = fitted_hatsa_model,
#'     voxel_timeseries_list = voxel_ts_list,
#'     voxel_coords = voxel_coords_map,
#'     parcel_coords = parcel_coords_map,
#'     n_nearest_parcels = 5, # Nystrom param
#'     kernel_sigma = "auto"    # Nystrom param
#'   )
#'
#'   # print(str(projected_vox_coeffs, max.level=1))
#'   # if (length(projected_vox_coeffs) > 0) {
#'   #   print(dim(projected_vox_coeffs[[1]])) # Should be T_time_vox x k_comp_fit
#'   # }
#' } else {
#'   if (interactive()) message("Skipping project_voxels example: fitted_hatsa_model not created.")
#' }
#' @export
#' @importFrom Matrix t
project_voxels.hatsa_projector <- function(object, 
                                           voxel_timeseries_list, 
                                           voxel_coords, 
                                           parcel_coords, 
                                           n_nearest_parcels = 10, 
                                           kernel_sigma = "auto", 
                                           data_type = c("timeseries", "coefficients"),
                                           ...) {
  
  data_type <- match.arg(data_type)
  
  if (!is.list(voxel_timeseries_list)) stop("`voxel_timeseries_list` must be a list.")
  num_new_subjects <- length(voxel_timeseries_list)
  if (num_new_subjects == 0) return(list())
  
  # Basic input validation for coordinates (critical for function operation)
  if (is.null(voxel_coords) || !is.matrix(voxel_coords) || ncol(voxel_coords) != 3) {
      stop("`voxel_coords` must be a V_v x 3 matrix.")
  }
  if (is.null(parcel_coords) || !is.matrix(parcel_coords) || ncol(parcel_coords) != 3 || nrow(parcel_coords) != object$parameters$V_p) {
      stop(sprintf("`parcel_coords` must be a V_p x 3 matrix, where V_p (%d) matches the model.", object$parameters$V_p))
  }

  # Perform heuristic coordinate system consistency checks (Ticket V-R1)
  .validate_coordinate_inputs(voxel_coords, parcel_coords)

  # Validate that the number of new subjects does not exceed N_subjects in the model
  # This simple check assumes the list corresponds to the first num_new_subjects subjects
  # A more robust solution would use named lists or explicit subject matching.
  if (num_new_subjects > object$parameters$N_subjects) {
    stop(sprintf("Number of subjects in voxel_timeseries_list (%d) exceeds N_subjects in model (%d).",
                 num_new_subjects, object$parameters$N_subjects))
  }

  V_v <- nrow(voxel_coords)
  k <- object$parameters$k

  results_list <- vector("list", num_new_subjects)

  for (i in 1:num_new_subjects) {
    X_voxel_ts_i <- voxel_timeseries_list[[i]]
    if (is.null(X_voxel_ts_i) || !is.matrix(X_voxel_ts_i) || ncol(X_voxel_ts_i) != V_v) {
      warning(sprintf("Voxel timeseries for subject %d is invalid or dimensions mismatch V_v (%d). Skipping.", i, V_v))
      results_list[[i]] <- matrix(NA, nrow = V_v, ncol = k)
      next
    }
    T_i <- nrow(X_voxel_ts_i)
    if (T_i == 0) {
        warning(sprintf("Voxel timeseries for subject %d has zero timepoints. Skipping.", i))
        results_list[[i]] <- matrix(NA, nrow = V_v, ncol = k)
        next
    }

    # Retrieve subject-specific components from the fitted object
    U_orig_parcel_i <- object$U_original_list[[i]]
    Lambda_orig_parcel_i <- object$Lambda_original_list[[i]]
    R_i <- object$R_final_list[[i]]

    if (is.null(U_orig_parcel_i) || is.null(Lambda_orig_parcel_i) || is.null(R_i)) {
      warning(sprintf("Missing original components (U, Lambda, or R) for subject %d in the model. Skipping voxel projection.", i))
      results_list[[i]] <- matrix(NA, nrow = T_i, ncol = k) # Note: T_i x k size
      next
    }
    
    # --- Rotation Matrix Validation (Ticket Z-R3) ---
    if (k > 0) { # Checks only meaningful if k > 0
        # Check for proper dimensions first
        if (nrow(R_i) != k || ncol(R_i) != k) {
            warning(sprintf("Rotation matrix R for subject %d has incorrect dimensions [%s], expected [%d x %d]. Skipping validation & projection.",
                            i, paste(dim(R_i), collapse="x"), k, k))
             results_list[[i]] <- matrix(NA, nrow = T_i, ncol = k) 
             next                       
        }
        
        # Check orthonormality (approximate)
        on_check <- crossprod(R_i)
        eye_k <- diag(k)
        if (!isTRUE(all.equal(on_check, eye_k, tolerance = 1e-6))) {
            # Calculate max absolute difference from identity
            max_diff <- max(abs(on_check - eye_k))
            warning(sprintf("Rotation matrix R for subject %d may not be orthonormal (max diff from identity: %.2e). Results may be affected.", i, max_diff))
        }
        
        # Check determinant (approximate, should be +1 for proper rotation)
        det_R <- det(R_i)
        if (!isTRUE(all.equal(det_R, 1.0, tolerance = 1e-6))) {
            warning(sprintf("Rotation matrix R for subject %d has determinant %.3f (expected ~1.0). May include reflection.", i, det_R))
        }
    }
    # --- End Validation ---

    # 1. Compute Nyström voxel basis Phi_voxel_i (V_v x k)
    Phi_voxel_i <- compute_voxel_basis_nystrom(
      voxel_coords = voxel_coords,
      parcel_coords = parcel_coords,
      U_orig_parcel = U_orig_parcel_i,
      Lambda_orig_parcel = Lambda_orig_parcel_i,
      n_nearest_parcels = n_nearest_parcels,
      kernel_sigma = kernel_sigma,
      ...
    )
    
    if (is.null(Phi_voxel_i) || nrow(Phi_voxel_i) != V_v || ncol(Phi_voxel_i) != k) {
        warning(sprintf("Nystrom basis computation failed or returned incorrect dimensions for subject %d. Skipping.", i))
        results_list[[i]] <- matrix(NA, nrow = T_i, ncol = k)
        next
    }

    # 2. Compute subject-specific voxel spectral coefficients: C = X %*% Phi
    C_voxel_coeffs_i <- X_voxel_ts_i %*% Phi_voxel_i
    
    # Apply scaling factor (if specified and applicable)
    scale_factor <- 1.0
    if (data_type == "timeseries") {
        # Apply the originally intended 1/T_i scaling, assuming it's desired for interpretation
        if (T_i > 0) { 
            scale_factor <- 1 / T_i 
        } else {
            scale_factor <- 0 # Avoid division by zero
        }
    }
    C_voxel_coeffs_i_scaled <- C_voxel_coeffs_i * scale_factor
    
    # 3. Rotate to common space: C_aligned = C_scaled %*% R_i
    # (T_i x k) %*% (k x k) -> (T_i x k)
    C_voxel_aligned_i <- C_voxel_coeffs_i_scaled %*% R_i
    
    results_list[[i]] <- as.matrix(C_voxel_aligned_i) # Ensure it's a base matrix (T_i x k)
  }

  return(results_list)
} 