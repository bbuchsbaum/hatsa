#' HATSA Projector Object
#'
#' An S3 object of class \code{hatsa_projector} that stores the results of a
#' Harmonized Tensors SVD Alignment (HATSA) analysis. This object inherits from
#' \code{multiblock_biprojector} (from the `multivarious` package) and is designed
#' to integrate HATSA outputs into a common framework for multiblock data analysis.
#'
#' @field v A numeric matrix (V_p x k) representing the mean aligned sketch,
#'   serving as the group-level template or common loadings.
#' @field s A numeric matrix ((N*V_p) x k) of stacked aligned sketches for all
#'   subjects. These are the subject-specific parcel scores in the common space.
#' @field sdev A numeric vector of length k, representing component-wise standard
#'   deviations (or scales). Currently defaults to a vector of 1s.
#' @field preproc A \code{pre_processor} object (from `multivarious`). For HATSA,
#'   this is typically `prep(pass())` as the input data to `run_hatsa_core` is
#'   already processed up to the point of raw time-series per subject.
#' @field block_indices A list defining which rows in the scores matrix `s` belong
#'   to which subject (block).
#' @field R_final_list A list of subject-specific rotation matrices (k x k) used to
#'   align each subject's original sketch to the common space.
#' @field U_original_list A list of subject-specific original (unaligned) sketch
#'   matrices (V_p x k) derived from their parcel-level graph Laplacians.
#' @field Lambda_original_list A list of numeric vectors, where each vector contains
#'   the k original eigenvalues corresponding to the eigenvectors in `U_original_list`
#'   for that subject. These are crucial for Nyström voxel projection.
#' @field Lambda_original_gaps_list A list of numeric vectors. Each vector contains
#'   the k-1 eigengap ratios `(λ_{j+1} - λ_j) / λ_j` for the corresponding subject's
#'   original eigenvalues. Useful for assessing spectral stability.
#' @field T_anchor_final A numeric matrix (V_a x k, where V_a is the number of
#'   anchors) representing the final group anchor template after Procrustes alignment.
#' @field parameters A list containing the input parameters used for the HATSA run
#'   (e.g., `k`, `V_p`, `N_subjects`, `anchor_indices`, `k_conn_pos`, `k_conn_neg`, `n_refine`).
#' @field method A character string, typically "hatsa_core", indicating the method
#'   used to generate the projector.
#' @field U_aligned_list (Internal) A list of subject-specific aligned sketch matrices (V_p x k).
#'   While `s` provides the stacked version, this list might be retained internally from the
#'   `run_hatsa_core` output passed to the constructor. For user access to aligned sketches per subject,
#'   one would typically use `project_block(object, block = i)` or segment `scores(object)`
#'   using `block_indices(object)`.
#'
#' @return A `hatsa_projector` object.
#'
#' @seealso \code{\link{run_hatsa_core}}, \code{\link{predict.hatsa_projector}}, \code{\link{project_voxels.hatsa_projector}}
#' @name hatsa_projector
NULL # This NULL is important for roxygen2 to generate class documentation

#' Constructor for hatsa_projector S3 class
#'
#' Creates a \code{hatsa_projector} object, which stores the results of the
#' HATSA algorithm and is designed to integrate with the \code{multivarious}
#' package, inheriting from \code{multiblock_biprojector}.
#'
#' @param hatsa_core_results A list containing the core outputs from the HATSA
#'   algorithm. Expected elements include:
#'   \itemize{
#'     \item{\code{U_aligned_list}: List of subject-specific aligned sketch matrices (V_p x k).}
#'     \item{\code{R_final_list}: List of subject-specific rotation matrices (k x k).}
#'     \item{\code{U_original_list}: List of subject-specific original sketch matrices (V_p x k).}
#'     \item{\code{Lambda_original_list}: List of subject-specific original eigenvalues (length k).}
#'     \item{\code{Lambda_original_gaps_list}: List of subject-specific eigengap ratios (length k-1).}
#'     \item{\code{T_anchor_final}: The final group anchor template matrix (N_anchors x k).}
#'   }
#' @param parameters A list of parameters used to run HATSA. Expected elements include:
#'   \itemize{
#'     \item{\code{k}: The number of spectral components (rank).}
#'     \item{\code{N_subjects}: The number of subjects.}
#'     \item{\code{V_p}: The number of parcels/vertices per subject.}
#'     \item{\code{method}: A string, typically \code{hatsa_core}.}
#'   }
#'
#' @return An object of class \code{c("hatsa_projector", "multiblock_biprojector", "projector", "list")}.
#'   This object contains:
#'   \itemize{
#'     \item{\code{v}: The group-level loading matrix (mean aligned sketch, V_p x k).}
#'     \item{\code{s}: The stacked scores matrix (concatenated aligned sketches, (N*V_p) x k).}
#'     \item{\code{sdev}: Component standard deviations (defaulted to 1s, length k).}
#'     \item{\code{preproc}: A preprocessing object, set to \code{multivarious::prep(multivarious::pass())}.}
#'     \item{\code{block_indices}: A list defining rows in \code{s} corresponding to each subject block.}
#'     \item{\code{R_final_list}: Stored from input.}
#'     \item{\code{U_original_list}: Stored from input.}
#'     \item{\code{Lambda_original_list}: Stored from input (crucial for voxel projection).}
#'     \item{\code{Lambda_original_gaps_list}: Stored from input.}
#'     \item{\code{T_anchor_final}: Stored from input.}
#'     \item{\code{parameters}: Stored from input.}
#'     \item{\code{method}: Stored from input parameters, typically \code{hatsa_core}.}
#'   }
#' @export
#' @examples
#' # This is a conceptual example, as real data structures are complex.
#' # Assuming hatsa_results and params are populated from run_hatsa_core:
#' # projector_obj <- hatsa_projector(
#' #   hatsa_core_results = list(
#' #     U_aligned_list = replicate(5, matrix(rnorm(100*10), 100, 10), simplify=FALSE),
#' #     R_final_list = replicate(5, diag(10), simplify=FALSE),
#' #     U_original_list = replicate(5, matrix(rnorm(100*10), 100, 10), simplify=FALSE),
#' #     Lambda_original_list = replicate(5, runif(10, 0.1, 1), simplify=FALSE), # example
#' #     Lambda_original_gaps_list = replicate(5, runif(9, 0.05, 0.5), simplify=FALSE), # example
#' #     T_anchor_final = matrix(rnorm(5*10), 5, 10)
#' #   ),
#' #   parameters = list(
#' #     k=10,
#' #     N_subjects=5,
#' #     V_p=100,
#' #     method="hatsa_core"
#' #   )
#' # )
#' # class(projector_obj)
#' # names(projector_obj)
hatsa_projector <- function(hatsa_core_results, parameters) {
  k <- parameters$k
  N_subjects <- parameters$N_subjects
  V_p <- parameters$V_p # Number of parcels/vertices per subject

  if (is.null(k) || is.null(N_subjects) || is.null(V_p)) {
    stop("Essential parameters (k, N_subjects, V_p) missing from input parameters.")
  }

  # Calculate v: Mean aligned sketch (V_p x k)
  # Robust handling of U_aligned_list, similar to task_hatsa_projector
  valid_aligned_sketches <- list()
  if (is.list(hatsa_core_results$U_aligned_list) && length(hatsa_core_results$U_aligned_list) > 0) {
      valid_aligned_sketches <- Filter(Negate(is.null), hatsa_core_results$U_aligned_list)
  }
  
  if (length(valid_aligned_sketches) > 0) {
      # Check consistency of valid sketches before reducing
      first_sketch_dim <- dim(valid_aligned_sketches[[1]])
      if (!all(sapply(valid_aligned_sketches, function(s) identical(dim(s), first_sketch_dim)))) {
          stop("Valid aligned sketches in U_aligned_list have inconsistent dimensions.")
      }
      if (first_sketch_dim[1] != V_p || first_sketch_dim[2] != k) {
          stop(sprintf("Dimensions of valid aligned sketches [%d x %d] do not match V_p=%d, k=%d.", 
                       first_sketch_dim[1], first_sketch_dim[2], V_p, k))
      }
      v_sum <- Reduce("+", valid_aligned_sketches)
      v <- v_sum / length(valid_aligned_sketches)
  } else {
      v <- matrix(0, nrow = V_p, ncol = k) # Default if no valid sketches
  }

  # Calculate s: Stacked aligned sketches ((N_subjects * V_p) x k)
  # Pre-allocate with NAs and fill, ensuring correct dimensions even with NULLs for some subjects
  s_full <- matrix(NA, nrow = N_subjects * V_p, ncol = k)
  if (N_subjects > 0 && V_p > 0) { # only proceed if dimensions make sense
      current_row_start <- 1
      if (is.list(hatsa_core_results$U_aligned_list) && length(hatsa_core_results$U_aligned_list) == N_subjects) {
          for (subj_idx in 1:N_subjects) {
              if (!is.null(hatsa_core_results$U_aligned_list[[subj_idx]])) {
                  # Additional check for individual sketch dimensions before assignment
                  current_sketch <- hatsa_core_results$U_aligned_list[[subj_idx]]
                  if (is.matrix(current_sketch) && nrow(current_sketch) == V_p && ncol(current_sketch) == k) {
                      s_full[current_row_start:(current_row_start + V_p - 1), ] <- current_sketch
                  } else {
                      warning(sprintf("Subject %d sketch has incorrect dimensions or is not a matrix. Filling with NAs in stacked scores.", subj_idx))
                  }
              }
              current_row_start <- current_row_start + V_p
          }
      } else if (N_subjects > 0) {
          warning("U_aligned_list is not a list of length N_subjects. Stacked scores 's' will be all NAs.")
      }
  }
  s <- s_full

  # sdev: Component standard deviations (length k) - default to 1s for HATSA
  sdev <- rep(1, k)

  # preproc: Set to pass() as HATSA alignment is complex and handled by predict
  preproc_obj <- multivarious::prep(multivarious::pass())

  # block_indices: List defining rows in s for each subject block
  block_indices <- if (N_subjects > 0 && V_p > 0) {
                        split(seq_len(N_subjects * V_p), rep(seq_len(N_subjects), each = V_p))
                   } else {
                        list()
                   }

  # Construct the projector object
  obj <- list(
    v = v,
    s = s,
    sdev = sdev,
    preproc = preproc_obj,
    block_indices = block_indices,
    R_final_list = hatsa_core_results$R_final_list,
    U_original_list = hatsa_core_results$U_original_list,
    Lambda_original_list = hatsa_core_results$Lambda_original_list,
    Lambda_original_gaps_list = hatsa_core_results$Lambda_original_gaps_list,
    T_anchor_final = hatsa_core_results$T_anchor_final,
    parameters = parameters,
    method = parameters$method # Store method, e.g., "hatsa_core"
  )

  # --- Initialize cache and store Fréchet mean of rotations ---
  obj$._cache <- list()
  valid_Rs_for_mean <- Filter(function(x) is.matrix(x) && !is.null(x), hatsa_core_results$R_final_list)
  if (length(valid_Rs_for_mean) > 0) {
    obj$._cache$R_frechet_mean <- tryCatch(
      frechet_mean_so_fast(valid_Rs_for_mean, refine = TRUE),
      error = function(e) if (k > 0) diag(k) else matrix(0, 0, 0)
    )
  } else {
    obj$._cache$R_frechet_mean <- if (k > 0) diag(k) else matrix(0, 0, 0)
  }

  # Set the class for S3 dispatch
  # Inherits from multiblock_biprojector, which itself inherits from projector
  class(obj) <- c("hatsa_projector", "multiblock_biprojector", "projector", "list")

  return(obj)
}

#' Print method for hatsa_projector objects
#'
#' Provides a concise summary of the hatsa_projector object.
#'
#' @param x A \code{hatsa_projector} object.
#' @param ... Additional arguments passed to print.
#' @export
print.hatsa_projector <- function(x, ...) {
  cat("HATSA Projector Object\n")
  cat("----------------------\n")
  cat("Method: ", x$parameters$method, "\n")
  cat("Number of Subjects (N): ", x$parameters$N_subjects, "\n")
  cat("Number of Parcels (V_p): ", x$parameters$V_p, "\n")
  cat("Number of Components (k): ", x$parameters$k, "\n")
  cat("\nKey components:\n")
  cat("  Mean Aligned Sketch (v): dimensions [", paste(dim(x$v), collapse = " x "), "]\n")
  cat("  Stacked Aligned Sketches (s): dimensions [", paste(dim(x$s), collapse = " x "), "]\n")
  cat("  Component Standard Deviations (sdev): length [", length(x$sdev), "]\n")
  invisible(x)
}

#' Extract coefficients (mean aligned sketch) from a hatsa_projector object
#'
#' @param object A \code{hatsa_projector} object.
#' @param ... Additional arguments (unused).
#' @return The mean aligned sketch matrix (v).
#' @export
coef.hatsa_projector <- function(object, ...) {
  return(object$v)
}

#' Extract scores (stacked aligned sketches) from a hatsa_projector object
#'
#' @param x A \code{hatsa_projector} object.
#' @param ... Additional arguments (unused).
#' @return The stacked aligned sketches matrix (s).
#' @export
scores.hatsa_projector <- function(x, ...) {
  return(x$s)
}

#' Extract component standard deviations from a hatsa_projector object
#'
#' @param x A \code{hatsa_projector} object.
#' @param ... Additional arguments (unused).
#' @return A vector of component standard deviations (sdev).
#' @export
sdev.hatsa_projector <- function(x, ...) {
  return(x$sdev)
}

#' Extract block indices from a hatsa_projector object
#'
#' These indices map rows of the scores matrix (s) to their respective blocks (subjects).
#'
#' @param x A \code{hatsa_projector} object.
#' @param ... Additional arguments (unused).
#' @return A list of block indices.
#' @export
block_indices.hatsa_projector <- function(x, ...) {
  return(x$block_indices)
}

#' Predict method for hatsa_projector objects
#'
#' Projects new subject-level parcel time-series data into the common space
#' defined by a fitted HATSA model.
#'
#' @param object A fitted \code{hatsa_projector} object.
#' @param newdata_list A list of new subject data. Each element of the list should be a
#'   numeric matrix representing parcel time-series data for one subject
#'   (T_i time points x V_p parcels). Parcel count (V_p) must match the original data.
#' @param ... Additional arguments (currently unused, for S3 compatibility).
#'
#' @return A list of aligned spectral sketch matrices (U_aligned_new, V_p x k),
#'   one for each subject in \code{newdata_list}.
#' @export
#' @importFrom stats predict
predict.hatsa_projector <- function(object, newdata_list, ...) {
  if (!is.list(newdata_list)) {
    stop("`newdata_list` must be a list of matrices.")
  }
  if (length(newdata_list) == 0) {
    return(list())
  }

  # Extract necessary parameters from the fitted object
  k <- object$parameters$k
  V_p_model <- object$parameters$V_p # V_p from the model
  anchor_indices <- object$parameters$anchor_indices
  T_anchor_final <- object$T_anchor_final
  k_conn_pos <- object$parameters$k_conn_pos
  k_conn_neg <- object$parameters$k_conn_neg
  use_dtw_model <- object$parameters$use_dtw

  aligned_sketches_new_list <- vector("list", length(newdata_list))
  for (i in seq_along(newdata_list)) {
    newdata_item <- newdata_list[[i]]
    if (!is.matrix(newdata_item)) {
      warning(sprintf("Item %d in newdata_list is not a matrix. Skipping.", i))
      aligned_sketches_new_list[i] <- list(NULL)
      next
    }
    if (ncol(newdata_item) != V_p_model) {
      warning(sprintf("Item %d in newdata_list has %d parcels, but model expects %d. Skipping.", 
                      i, ncol(newdata_item), V_p_model))
      aligned_sketches_new_list[i] <- list(NULL)
      next
    }

    # Generate parcel names if not present, consistent with run_hatsa_core logic
    current_pnames <- colnames(newdata_item)
    if (is.null(current_pnames) || length(current_pnames) != V_p_model) {
      current_pnames <- paste0("P", 1:V_p_model) 
    }

    # 1. Compute subject connectivity graph
    W_conn_new <- compute_subject_connectivity_graph_sparse(newdata_item, 
                                                            current_pnames, 
                                                            k_conn_pos, 
                                                            k_conn_neg, 
                                                            use_dtw = use_dtw_model)
    # 2. Compute graph Laplacian
    L_conn_new <- compute_graph_laplacian_sparse(W_conn_new)

    # 3. Compute original spectral sketch (eigenvectors)
    # compute_spectral_sketch_sparse returns a list(vectors=U, values=Lambda)
    sketch_result_new <- compute_spectral_sketch_sparse(L_conn_new, k)
    U_orig_new <- sketch_result_new$vectors

    if (is.null(U_orig_new) || nrow(U_orig_new) != V_p_model || ncol(U_orig_new) != k) {
        warning(sprintf("Subject %d: Original sketch computation failed or has incorrect dimensions. Skipping.", i))
        aligned_sketches_new_list[[i]] <- matrix(NA, nrow = V_p_model, ncol = k)
        next
    }

    # 4. Extract anchors from U_orig_new
    if (length(anchor_indices) == 0 && k > 0) {
        warning(sprintf("Subject %d: No anchor indices defined in the model, but k > 0. Cannot perform alignment. Returning original sketch.", i))
        aligned_sketches_new_list[[i]] <- U_orig_new # Or handle as error
        next
    } else if (k == 0) {
        aligned_sketches_new_list[[i]] <- U_orig_new # k=0 sketch
        next
    }
    A_orig_new <- U_orig_new[anchor_indices, , drop = FALSE]

    # 5. Find rotation R_new by aligning A_orig_new to stored T_anchor_final
    # Ensure T_anchor_final is correctly dimensioned
    if (is.null(T_anchor_final) || nrow(T_anchor_final) != length(anchor_indices) || ncol(T_anchor_final) != k) {
        warning(sprintf("Subject %d: T_anchor_final in the model is invalid. Cannot perform alignment. Returning original sketch.", i))
        aligned_sketches_new_list[[i]] <- U_orig_new
        next
    }
    R_new <- solve_procrustes_rotation(A_orig_new, T_anchor_final)

    # 6. Compute aligned sketch U_aligned_new
    U_aligned_new <- U_orig_new %*% R_new
    
    aligned_sketches_new_list[[i]] <- U_aligned_new
  }

  # Directly return
  return(aligned_sketches_new_list)
}

#' Project a specific block (subject) using a hatsa_projector object
#'
#' Retrieves the aligned spectral sketch for a specified subject. If new data
#' for that subject is provided, it projects that new data. Otherwise, it extracts
#' the stored aligned sketch from the original analysis.
#'
#' @param object A fitted \code{hatsa_projector} object.
#' @param newdata Optional. A single numeric matrix of parcel time-series data
#'   (T_i x V_p) for the subject specified by \code{block}. If NULL, the stored
#'   sketch for the block is returned. If provided, this new data is projected.
#' @param block An integer, the index of the subject (block) to project.
#' @param ... Additional arguments passed to \code{predict.hatsa_projector} if
#'   \code{newdata} is provided.
#'
#' @return A matrix (V_p x k) representing the aligned spectral sketch for the
#'   specified subject/block. Returns NULL or raises an error if the block is invalid
#'   or data is inappropriate.
#' @importFrom multivarious project_block
#' @export
project_block.hatsa_projector <- function(object, newdata = NULL, block, ...) {
  V_p <- object$parameters$V_p
  k <- object$parameters$k

  if (is.null(newdata)) {
    # Retrieve stored sketch for the block
    if (!is.numeric(block) || length(block) != 1 || block < 1 || block > object$parameters$N_subjects) {
      stop(sprintf("Invalid block index. Must be an integer between 1 and %d (N_subjects).", object$parameters$N_subjects))
    }
    if (is.null(object$s) || is.null(object$block_indices) || length(object$block_indices) < block) {
        stop("Stored scores 's' or 'block_indices' are missing or incomplete in the object.")
    }

    block_idx_rows <- object$block_indices[[block]]
    if (is.null(block_idx_rows) || length(block_idx_rows) != V_p) {
        stop(sprintf("Block indices for block %d are invalid or do not match V_p.", block))
    }
    
    # Ensure scores matrix has enough rows and correct number of columns
    if (nrow(object$s) < max(block_idx_rows) || ncol(object$s) != k) {
        stop("Scores matrix 's' has incorrect dimensions for the requested block or k.")
    }

    subject_sketch <- object$s[block_idx_rows, , drop = FALSE]
    
    # Reshape to V_p x k if it was extracted as a vector (e.g. k=1)
    # Although drop=FALSE should prevent this for matrices, good to be robust.
    # However, with block_indices, object$s is already (N*V_p) x k, so subsetting rows gives (V_p) x k.
    # So, explicit dim check might be more direct.
    if (nrow(subject_sketch) != V_p || ncol(subject_sketch) != k) {
        stop(sprintf("Extracted sketch for block %d has dimensions [%s] but expected [%d x %d].",
                     block, paste(dim(subject_sketch), collapse="x"), V_p, k))
    }
    return(subject_sketch)

  } else {
    # Project new data for the specified block (newdata is for this one block)
    if (!is.matrix(newdata)) {
      stop("`newdata` must be a numeric matrix (time-series data for one subject).")
    }
    if (ncol(newdata) != V_p) {
      stop(sprintf("`newdata` has %d parcels, but model expects %d.", ncol(newdata), V_p))
    }

    # predict.hatsa_projector expects a list of matrices
    # It returns a list of aligned sketches.
    aligned_sketch_list <- tryCatch(
      predict.hatsa_projector(object, newdata_list = list(newdata), ...),
      error = function(e) {
        warning(sprintf("Prediction for the new data block failed: %s", e$message))
        list(NULL)
      }
    )
    
    if (length(aligned_sketch_list) == 1 && !is.null(aligned_sketch_list[[1]])) {
      return(aligned_sketch_list[[1]])
    } else {
      warning("Prediction for the new data block failed or returned an unexpected result.")
      # Return an empty or NA matrix of appropriate size
      return(matrix(NA, nrow = V_p, ncol = k))
    }
  }
}

#' Summary method for hatsa_projector objects
#'
#' Provides a detailed summary of the hatsa_projector object, including
#' mean anchor alignment error for the original subjects.
#'
#' @param object A \code{hatsa_projector} object.
#' @param ... Additional arguments (unused).
#' @param compute_riemannian_dispersion Logical, if TRUE, computes and includes
#'   Riemannian dispersion metrics for SPD representations (default: FALSE).
#' @param riemannian_dispersion_type Character string, type of SPD representation
#'   for dispersion (e.g., "cov_coeffs"). Passed to `riemannian_dispersion_spd`.
#' @param riemannian_dispersion_options A list of additional arguments passed to
#'   `riemannian_dispersion_spd` (e.g., for `use_geometric_median`, `spd_regularize_epsilon`, `verbose`).
#' @param recompute_R_bar Logical; if TRUE the Fr\u00e9chet mean of the rotations
#'   is recomputed even if a cached value is present. Default FALSE.
#' @return A list object of class \code{summary.hatsa_projector} containing
#'   summary statistics.
#' @export
#' @importFrom stats predict sd median na.omit
#' @importFrom expm logm expm sqrtm
summary.hatsa_projector <- function(object, ...,
                                    compute_riemannian_dispersion = FALSE,
                                    riemannian_dispersion_type = "cov_coeffs",
                                    riemannian_dispersion_options = list(),
                                    recompute_R_bar = FALSE) {
  # Basic information
  summary_info <- list(
    method = object$parameters$method,
    N_subjects = object$parameters$N_subjects,
    V_p = object$parameters$V_p,
    k = object$parameters$k,
    num_anchors = length(object$parameters$anchor_indices)
  )

  # Calculate mean anchor alignment error
  if (object$parameters$N_subjects > 0 && object$parameters$k > 0 && summary_info$num_anchors > 0) {
    anchor_errors <- numeric(object$parameters$N_subjects)
    U_original_list <- object$U_original_list
    R_final_list <- object$R_final_list
    T_anchor_final <- object$T_anchor_final
    anchor_indices <- object$parameters$anchor_indices

    all_data_valid_for_error_calc <- TRUE
    if (is.null(U_original_list) || length(U_original_list) != object$parameters$N_subjects ||
        is.null(R_final_list) || length(R_final_list) != object$parameters$N_subjects ||
        is.null(T_anchor_final) || is.null(anchor_indices)) {
      all_data_valid_for_error_calc <- FALSE
    }

    if (all_data_valid_for_error_calc) {
      for (i in 1:object$parameters$N_subjects) {
        U_orig_i <- U_original_list[[i]]
        R_i <- R_final_list[[i]]

        if (is.null(U_orig_i) || is.null(R_i) || 
            nrow(U_orig_i) != object$parameters$V_p || ncol(U_orig_i) != object$parameters$k ||
            nrow(R_i) != object$parameters$k || ncol(R_i) != object$parameters$k || 
            nrow(T_anchor_final) != summary_info$num_anchors || ncol(T_anchor_final) != object$parameters$k) {
          anchor_errors[i] <- NA # Mark as NA if data is inconsistent
          next
        }
        
        A_orig_i_anchors <- U_orig_i[anchor_indices, , drop = FALSE]
        A_aligned_i_anchors <- A_orig_i_anchors %*% R_i
        
        # Frobenius norm of the difference
        # Ensure dimensions match for subtraction
        if (nrow(A_aligned_i_anchors) == nrow(T_anchor_final) && ncol(A_aligned_i_anchors) == ncol(T_anchor_final)) {
            error_matrix <- A_aligned_i_anchors - T_anchor_final
            anchor_errors[i] <- norm(error_matrix, type = "F")
        } else {
            anchor_errors[i] <- NA # Inconsistent dimensions for error calculation
        }
      }
      summary_info$mean_anchor_alignment_error <- mean(anchor_errors, na.rm = TRUE)
      summary_info$anchor_alignment_errors_per_subject <- anchor_errors
    } else {
      summary_info$mean_anchor_alignment_error <- NA
      summary_info$anchor_alignment_errors_per_subject <- rep(NA, object$parameters$N_subjects)
      warning("Could not calculate anchor alignment error due to missing or incomplete data in the object.")
    }
  } else {
    summary_info$mean_anchor_alignment_error <- NA
    summary_info$anchor_alignment_errors_per_subject <- rep(NA, object$parameters$N_subjects)
  }

  # --- HMET-001: Rotation Dispersion (sigma_R) --- 
  if (object$parameters$N_subjects > 0 && object$parameters$k > 0 && !is.null(object$R_final_list)) {
    k_dim <- object$parameters$k
    Rs <- object$R_final_list
    # Filter out any NULLs from Rs, which might occur if a subject had issues
    Rs_valid <- Filter(Negate(is.null), Rs)
    Rs_valid <- Filter(function(R_mat) is.matrix(R_mat) && all(dim(R_mat) == c(k_dim, k_dim)), Rs_valid)
    N_valid_subj_for_rot <- length(Rs_valid)

    if (N_valid_subj_for_rot > 1) {

      if (!recompute_R_bar && !is.null(object$._cache$R_frechet_mean)) {
        R_bar <- object$._cache$R_frechet_mean
      } else {
        R_bar <- tryCatch(
          frechet_mean_so_fast(Rs_valid, refine = TRUE),
          error = function(e) {
            warning("Failed to compute Fr\u00e9chet mean rotation: ", e$message, ". Using identity matrix as fallback.")
            if (k_dim > 0) diag(k_dim) else matrix(0, 0, 0)
          }
        )
        object$._cache$R_frechet_mean <- R_bar
      }
      
      summary_info$R_frechet_mean <- R_bar # Store the mean rotation
      
      geodesic_distances_to_mean <- sapply(Rs_valid, function(R_i) geodesic_dist_so_k(R_i, R_bar))
      
      summary_info$rotation_dispersion_mean_geo_dist <- mean(geodesic_distances_to_mean, na.rm = TRUE)
      summary_info$rotation_dispersion_sd_geo_dist <- sd(geodesic_distances_to_mean, na.rm = TRUE)
      summary_info$rotation_geodesic_distances <- geodesic_distances_to_mean
    } else {
      summary_info$R_frechet_mean <- if (N_valid_subj_for_rot == 1 && k_dim > 0) Rs_valid[[1]] else if (k_dim > 0) diag(k_dim) else matrix(0,0,0)
      summary_info$rotation_dispersion_mean_geo_dist <- NA
      summary_info$rotation_dispersion_sd_geo_dist <- NA
      summary_info$rotation_geodesic_distances <- numeric(0)
    }
  } else {
    summary_info$R_frechet_mean <- if(object$parameters$k > 0) diag(object$parameters$k) else matrix(0,0,0)
    summary_info$rotation_dispersion_mean_geo_dist <- NA
    summary_info$rotation_dispersion_sd_geo_dist <- NA
    summary_info$rotation_geodesic_distances <- numeric(0)
  }

  # --- HMET-001: Variance Explained by k Components (sum of selected eigenvalues) ---
  if (object$parameters$N_subjects > 0 && object$parameters$k > 0 && !is.null(object$Lambda_original_list)) {
    lambdas_orig_valid <- Filter(Negate(is.null), object$Lambda_original_list)
    lambdas_orig_valid <- Filter(function(l) is.numeric(l) && length(l) == object$parameters$k, lambdas_orig_valid)
    
    if (length(lambdas_orig_valid) > 0) {
      sum_k_lambdas_per_subject <- sapply(lambdas_orig_valid, sum)
      summary_info$mean_sum_k_selected_eigenvalues <- mean(sum_k_lambdas_per_subject, na.rm = TRUE)
      summary_info$sd_sum_k_selected_eigenvalues <- sd(sum_k_lambdas_per_subject, na.rm = TRUE)
      summary_info$sum_k_selected_eigenvalues_per_subject <- sum_k_lambdas_per_subject
    } else {
      summary_info$mean_sum_k_selected_eigenvalues <- NA
      summary_info$sd_sum_k_selected_eigenvalues <- NA
      summary_info$sum_k_selected_eigenvalues_per_subject <- numeric(0)
    }
  } else {
    summary_info$mean_sum_k_selected_eigenvalues <- NA
    summary_info$sd_sum_k_selected_eigenvalues <- NA
    summary_info$sum_k_selected_eigenvalues_per_subject <- numeric(0)
  }

  # --- HMET-001: Eigengap Ratios ---
  if (object$parameters$N_subjects > 0 && object$parameters$k > 1 && !is.null(object$Lambda_original_gaps_list)) {
    gaps_list_valid <- Filter(Negate(is.null), object$Lambda_original_gaps_list)
    # Each element of gaps_list_valid should be a numeric vector of length k-1
    gaps_list_valid <- Filter(function(g) is.numeric(g) && length(g) == (object$parameters$k - 1), gaps_list_valid)

    if (length(gaps_list_valid) > 0) {
      # Focus on the (k-1)th gap: (lambda_k - lambda_{k-1}) / lambda_{k-1}
      # This is the last element of each vector in gaps_list_valid
      k_minus_1_gaps <- sapply(gaps_list_valid, function(gaps_subj) gaps_subj[object$parameters$k - 1])
      
      summary_info$median_k_minus_1_eigengap_ratio <- median(k_minus_1_gaps, na.rm = TRUE)
      summary_info$min_k_minus_1_eigengap_ratio <- if(any(!is.na(k_minus_1_gaps))) min(k_minus_1_gaps, na.rm = TRUE) else NA
      summary_info$max_k_minus_1_eigengap_ratio <- if(any(!is.na(k_minus_1_gaps))) max(k_minus_1_gaps, na.rm = TRUE) else NA
      summary_info$k_minus_1_eigengap_ratios_per_subject <- k_minus_1_gaps
    } else {
      summary_info$median_k_minus_1_eigengap_ratio <- NA
      summary_info$min_k_minus_1_eigengap_ratio <- NA
      summary_info$max_k_minus_1_eigengap_ratio <- NA
      summary_info$k_minus_1_eigengap_ratios_per_subject <- numeric(0)
    }
  } else {
    summary_info$median_k_minus_1_eigengap_ratio <- NA
    summary_info$min_k_minus_1_eigengap_ratio <- NA
    summary_info$max_k_minus_1_eigengap_ratio <- NA
    summary_info$k_minus_1_eigengap_ratios_per_subject <- numeric(0)
  }

  # --- RGEOM-006: Optionally compute Riemannian Dispersion ---
  if (compute_riemannian_dispersion) {
    # Use tryCatch to safely attempt to compute Riemannian dispersion
    verbose <- riemannian_dispersion_options$verbose %||% FALSE
    if (verbose && interactive()) { 
        message(sprintf("Computing Riemannian dispersion for type: %s...", riemannian_dispersion_type))
    }
    
    # Prepare arguments for riemannian_dispersion_spd
    disp_args <- list(
      object = object,
      type = riemannian_dispersion_type
    )
    
    # Add options from riemannian_dispersion_options
    if (length(riemannian_dispersion_options) > 0) {
        for(opt_name in names(riemannian_dispersion_options)) {
            disp_args[[opt_name]] <- riemannian_dispersion_options[[opt_name]]
        }
    }
    
    dispersion_results <- tryCatch({
      # Try to call riemannian_dispersion_spd function
      do.call(riemannian_dispersion_spd, disp_args)
    }, error = function(e) {
      # Handle missing function or other errors
      warning("Failed to compute Riemannian dispersion: ", e$message)
      NULL
    })
    
    if (!is.null(dispersion_results)) {
      summary_info$riemannian_dispersion_type_used <- riemannian_dispersion_type
      summary_info$riemannian_mean_dispersion <- dispersion_results$mean_dispersion
      summary_info$riemannian_median_dispersion <- dispersion_results$median_dispersion
      summary_info$riemannian_num_valid_subjects_for_disp <- dispersion_results$num_valid_subjects
    } else {
      summary_info$riemannian_dispersion_type_used <- riemannian_dispersion_type
      summary_info$riemannian_mean_dispersion <- NA_real_
      summary_info$riemannian_median_dispersion <- NA_real_
      summary_info$riemannian_num_valid_subjects_for_disp <- 0
    }
  } else {
    # Ensure these fields are NULL if not computed, for consistent print checks
    summary_info$riemannian_dispersion_type_used <- NULL
    summary_info$riemannian_mean_dispersion <- NULL
    summary_info$riemannian_median_dispersion <- NULL
    summary_info$riemannian_num_valid_subjects_for_disp <- NULL
  }

  class(summary_info) <- "summary.hatsa_projector"
  return(summary_info)
}


#' Print method for summary.hatsa_projector objects
#'
#' @param x A \code{summary.hatsa_projector} object.
#' @param ... Additional arguments (unused).
#' @export
print.summary.hatsa_projector <- function(x, ...) {
  cat("HATSA Projector Summary\n")
  cat("-----------------------\n")
  cat("Method: ", x$method, "\n")
  cat("Number of Subjects (N): ", x$N_subjects, "\n")
  cat("Number of Parcels (V_p): ", x$V_p, "\n")
  cat("Number of Components (k): ", x$k, "\n")
  cat("Number of Anchors: ", x$num_anchors, "\n")
  cat("\n")
  if (!is.null(x$mean_anchor_alignment_error) && !is.na(x$mean_anchor_alignment_error)) {
    cat("Mean Anchor Alignment Error (Frobenius Norm): ", sprintf("%.4f", x$mean_anchor_alignment_error), "\n")
  } else if (x$k > 0 && x$num_anchors > 0 && x$N_subjects > 0) {
    cat("Mean Anchor Alignment Error: Not calculable (data missing or inconsistent)\n")
  }

  # Display Rotation Dispersion
  if (!is.null(x$rotation_dispersion_mean_geo_dist) && !is.na(x$rotation_dispersion_mean_geo_dist)) {
    cat("Rotation Dispersion (Mean Geodesic Distance to R_bar): ", sprintf("%.4f", x$rotation_dispersion_mean_geo_dist), "\n")
    cat("Rotation Dispersion (SD Geodesic Distance to R_bar): ", sprintf("%.4f", x$rotation_dispersion_sd_geo_dist), "\n")
  }

  # Display Variance Explained (Sum of k eigenvalues)
  if (!is.null(x$mean_sum_k_selected_eigenvalues) && !is.na(x$mean_sum_k_selected_eigenvalues)) {
    cat("Mean Sum of k Selected Eigenvalues (Energy): ", sprintf("%.4f", x$mean_sum_k_selected_eigenvalues), "\n")
    cat("SD Sum of k Selected Eigenvalues: ", sprintf("%.4f", x$sd_sum_k_selected_eigenvalues), "\n")
  }

  # Display Eigengap Ratios
  if (!is.null(x$median_k_minus_1_eigengap_ratio) && !is.na(x$median_k_minus_1_eigengap_ratio)) {
    cat(sprintf("Median (λ_k - λ_{k-1})/λ_{k-1} Ratio: %.4f (Min: %.4f, Max: %.4f)\n", 
                x$median_k_minus_1_eigengap_ratio, 
                x$min_k_minus_1_eigengap_ratio, 
                x$max_k_minus_1_eigengap_ratio))
  }

  # Display Riemannian Dispersion if computed
  if (!is.null(x$riemannian_dispersion_type_used)) {
    cat("\n--- Riemannian Dispersion ---\n")
    cat("SPD Representation Type: ", x$riemannian_dispersion_type_used, "\n")
    if (!is.na(x$riemannian_mean_dispersion)) {
      cat("  Mean Riemannian Dispersion: ", sprintf("%.4f", x$riemannian_mean_dispersion), "\n")
      cat("  Median Riemannian Dispersion: ", sprintf("%.4f", x$riemannian_median_dispersion), "\n")
      cat("  (Based on ", x$riemannian_num_valid_subjects_for_disp, " valid subjects)\n")
    } else {
      cat("  Riemannian Dispersion: Not calculable or computation failed.\n")
    }
  }

  invisible(x)
}

# --- HMET-002: Reconstruction Error S3 Method ---

#' Calculate Reconstruction Error for HATSA Projector Objects
#'
#' Computes various types of reconstruction errors based on a fitted HATSA model.
#' This helps quantify alignment quality, sanity check model components, or assess
#' how well non-anchor information is captured.
#'
#' @param object A fitted \code{hatsa_projector} object.
#' @param type A character string specifying the type of reconstruction error to compute.
#'   One of:
#'   \itemize{
#'     \item{\code{"anchors"} (Default): Error between subject-specific aligned anchor sketches
#'           and the group anchor template (`T_anchor_final`).}
#'     \item{\code{"all_parcels"}: Error between original subject sketches (`U_original_list`)
#'           and sketches reconstructed from aligned versions via inverse rotation.
#'           (Sanity check, should be near zero if rotations are orthogonal).}
#'     \item{\code{"non_anchors"}: Error in predicting non-anchor parcel sketches from
#'           anchor parcel sketches in the aligned space (Anchor Residual `ε̄`).}
#'   }
#' @param ... Additional arguments, potentially passed to internal helper functions
#'   (e.g., for non-anchor prediction methods).
#'
#' @return A named list containing reconstruction error metrics. The structure depends
#'   on the `type` requested, typically including `mean_error` and `per_subject_error`.
#'   For `type="non_anchors"`, it also includes `per_parcel_error`.
#'
#' @examples
#' # Assuming `fit_obj` is a hatsa_projector from `run_hatsa_core()`
#' # Anchor reconstruction error
#' # reconstruction_error(fit_obj, type = "anchors")
#' # Sanity check using all parcels
#' # reconstruction_error(fit_obj, type = "all_parcels")
#' # Predict non-anchor parcels from anchors
#' # reconstruction_error(fit_obj, type = "non_anchors")
#' @export
reconstruction_error <- function(object, type = "anchors", ...) {
  UseMethod("reconstruction_error")
}

#' @rdname reconstruction_error
#' @export
reconstruction_error.hatsa_projector <- function(object, type = "anchors", ...) {
  # Parameter extraction from object
  N_subjects <- object$parameters$N_subjects
  V_p <- object$parameters$V_p
  k <- object$parameters$k
  anchor_indices <- object$parameters$anchor_indices
  num_anchors <- length(anchor_indices)
  
  U_original_list <- object$U_original_list
  R_final_list <- object$R_final_list
  T_anchor_final <- object$T_anchor_final

  # Initialize results list
  results <- list(
    type = type,
    mean_error = NA,
    per_subject_error = rep(NA, N_subjects),
    notes = character(0)
  )

  # --- Helper to get U_aligned_list (reconstruct from s if not directly stored) ---
  .get_U_aligned_list <- function(obj) {
    # Check if U_aligned_list was explicitly stored by constructor (current constructor does not)
    # if (!is.null(obj$._cache$U_aligned_list)) return(obj$._cache$U_aligned_list)
    # if (!is.null(obj$U_aligned_list_stored)) return(obj$U_aligned_list_stored) # If we decide to store it

    # Reconstruct from obj$s and obj$block_indices
    if (is.null(obj$s) || is.null(obj$block_indices) || obj$parameters$N_subjects == 0) {
      warning("Cannot reconstruct U_aligned_list: object$s or object$block_indices missing, or N_subjects is 0.")
      return(vector("list", obj$parameters$N_subjects)) # List of NULLs
    }
    
    U_aligned_list_recon <- vector("list", obj$parameters$N_subjects)
    for (i in 1:obj$parameters$N_subjects) {
      if (length(obj$block_indices) >= i && !is.null(obj$block_indices[[i]])) {
        rows_subj_i <- obj$block_indices[[i]]
        if (max(rows_subj_i) <= nrow(obj$s) && ncol(obj$s) == obj$parameters$k) {
          U_aligned_list_recon[[i]] <- obj$s[rows_subj_i, , drop = FALSE]
        } else {
          U_aligned_list_recon[[i]] <- NULL # Mark as NULL if dimensions mismatch
        }
      } else {
        U_aligned_list_recon[[i]] <- NULL
      }
    }
    # Optional: Cache this reconstructed list if it were to be used multiple times within this function call
    # obj$._cache$U_aligned_list <- U_aligned_list_recon 
    return(U_aligned_list_recon)
  }

  # Type dispatch
  if (type == "anchors") {
    if (N_subjects == 0 || k == 0 || num_anchors == 0) {
      results$notes <- "Skipped: N_subjects, k, or num_anchors is 0."
      return(results)
    }
    if (is.null(U_original_list) || is.null(R_final_list) || is.null(T_anchor_final)) {
      results$notes <- "Skipped: U_original_list, R_final_list, or T_anchor_final is missing."
      return(results)
    }

    subject_errors <- numeric(N_subjects)
    for (i in 1:N_subjects) {
      U_orig_i <- U_original_list[[i]]
      R_i <- R_final_list[[i]]

      if (is.null(U_orig_i) || !is.matrix(U_orig_i) || nrow(U_orig_i) != V_p || ncol(U_orig_i) != k ||
          is.null(R_i) || !is.matrix(R_i) || any(dim(R_i) != c(k, k)) || 
          !is.matrix(T_anchor_final) || nrow(T_anchor_final) != num_anchors || ncol(T_anchor_final) != k) {
        subject_errors[i] <- NA
        next
      }
      
      A_orig_i_anchors <- U_orig_i[anchor_indices, , drop = FALSE]
      if(any(is.na(A_orig_i_anchors))) {subject_errors[i] <- NA; next} # Check after subsetting

      A_aligned_i_anchors <- A_orig_i_anchors %*% R_i
      
      if (nrow(A_aligned_i_anchors) == num_anchors && ncol(A_aligned_i_anchors) == k) {
        error_matrix <- A_aligned_i_anchors - T_anchor_final
        subject_errors[i] <- norm(error_matrix, type = "F")
      } else {
        subject_errors[i] <- NA
      }
    }
    results$per_subject_error <- subject_errors
    results$mean_error <- mean(subject_errors, na.rm = TRUE)
    
  } else if (type == "all_parcels") {
    if (N_subjects == 0 || k == 0) {
      results$notes <- "Skipped: N_subjects or k is 0."
      return(results)
    }
    if (is.null(U_original_list) || is.null(R_final_list)) {
      results$notes <- "Skipped: U_original_list or R_final_list missing."
      return(results)
    }

    U_aligned_list_local <- .get_U_aligned_list(object)
    subject_errors <- numeric(N_subjects)
    for (i in seq_len(N_subjects)) {
      U_orig_i <- U_original_list[[i]]
      R_i <- R_final_list[[i]]
      U_aligned_i <- U_aligned_list_local[[i]]

      if (is.null(U_orig_i) || is.null(U_aligned_i) || is.null(R_i) ||
          !is.matrix(U_orig_i) || !is.matrix(U_aligned_i) ||
          nrow(U_orig_i) != V_p || nrow(U_aligned_i) != V_p ||
          ncol(U_orig_i) != k || ncol(U_aligned_i) != k ||
          !is.matrix(R_i) || any(dim(R_i) != c(k, k))) {
        subject_errors[i] <- NA
        next
      }

      U_recon_i <- U_aligned_i %*% t(R_i)
      subject_errors[i] <- norm(U_orig_i - U_recon_i, type = "F")
    }

    results$per_subject_error <- subject_errors
    results$mean_error <- mean(subject_errors, na.rm = TRUE)

  } else if (type == "non_anchors") {
    if (N_subjects == 0 || k == 0 || num_anchors == 0) {
      results$notes <- "Skipped: N_subjects, k, or num_anchors is 0."
      return(results)
    }
    if (is.null(U_original_list)) {
      results$notes <- "Skipped: U_original_list missing."
      return(results)
    }

    U_aligned_list_local <- .get_U_aligned_list(object)
    non_anchor_indices <- setdiff(seq_len(V_p), anchor_indices)
    n_non <- length(non_anchor_indices)

    subject_errors <- numeric(N_subjects)
    per_parcel_err_accum <- rep(0, n_non)
    valid_counts <- 0
    lambda <- 1e-4

    for (i in seq_len(N_subjects)) {
      U_orig_i <- U_original_list[[i]]
      U_aligned_i <- U_aligned_list_local[[i]]

      if (is.null(U_orig_i) || is.null(U_aligned_i) ||
          !is.matrix(U_orig_i) || !is.matrix(U_aligned_i) ||
          nrow(U_orig_i) != V_p || ncol(U_orig_i) != k ||
          nrow(U_aligned_i) != V_p || ncol(U_aligned_i) != k) {
        subject_errors[i] <- NA
        next
      }

      A_orig <- U_orig_i[anchor_indices, , drop = FALSE]
      NA_orig <- U_orig_i[non_anchor_indices, , drop = FALSE]
      A_align <- U_aligned_i[anchor_indices, , drop = FALSE]
      NA_align <- U_aligned_i[non_anchor_indices, , drop = FALSE]

      G <- A_orig %*% t(A_orig) + lambda * diag(num_anchors)
      W <- NA_orig %*% t(A_orig) %*% solve(G)
      NA_pred <- W %*% A_align

      row_err <- sqrt(rowSums((NA_align - NA_pred)^2))
      subject_errors[i] <- mean(row_err)

      if (all(!is.na(row_err))) {
        per_parcel_err_accum <- per_parcel_err_accum + row_err
        valid_counts <- valid_counts + 1
      }
    }

    results$per_subject_error <- subject_errors
    results$mean_error <- mean(subject_errors, na.rm = TRUE)
    if (valid_counts > 0) {
      results$per_parcel_error <- per_parcel_err_accum / valid_counts
    } else {
      results$per_parcel_error <- rep(NA_real_, n_non)
    }

  } else {
    stop(sprintf("Unknown reconstruction error type: '%s'. Must be 'anchors', 'all_parcels', or 'non_anchors'.", type))
  }

  return(results)
} 