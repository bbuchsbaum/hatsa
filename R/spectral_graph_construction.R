#' Compute pairwise correlations in manageable chunks
#'
#' Helper that computes Pearson correlations between all columns of
#' `X_centered` in blocks to avoid allocating the full `V_p \times V_p`
#' dense matrix. Each block contains `block_size` columns and returns a list of
#' correlation submatrices.
#'
#' @param X_centered A centered numeric matrix (`T_i` x `V_p`).
#' @param col_sds Vector of column standard deviations for `X_centered`.
#' @param block_size Number of columns per block. Defaults to 50.
#' @return A list where each element contains:
#'   \itemize{
#'     \item `indices`: integer vector of column indices covered by the block.
#'     \item `corr`: matrix of correlations (`length(indices)` x `V_p`).
#'   }
#' @importFrom Matrix crossprod
#' @keywords internal
chunked_sparse_cor <- function(X_centered, col_sds, block_size = 50L) {
  V_p <- ncol(X_centered)
  if (V_p == 0) return(list())
  T_i <- nrow(X_centered)
  block_starts <- seq(1, V_p, by = block_size)
  res <- vector("list", length(block_starts))
  idx <- 1L
  for (start in block_starts) {
    end <- min(start + block_size - 1L, V_p)
    cov_block <- Matrix::crossprod(X_centered[, start:end, drop = FALSE], X_centered) / (T_i - 1)
    corr_block <- sweep(cov_block, 2, col_sds, "/")
    corr_block <- sweep(corr_block, 1, col_sds[start:end], "/")
    corr_block[!is.finite(corr_block)] <- 0
    res[[idx]] <- list(indices = start:end, corr = corr_block)
    idx <- idx + 1L
  }
  res
}

#' Compute subject-specific sparse connectivity graph `W_conn_i`
#'
#' Calculates the sparse connectivity graph for a single subject.
#' Steps:
#' 1. Compute pairwise Pearson correlations using a memory-efficient
#'    cross-product approach that avoids constructing a dense `V_p` \times `V_p` matrix.
#' 2. Identify and mask zero-variance parcels.
#' 3. Sparsify: For each parcel, identify indices of `k_conn_pos` strongest
#'    positive and `k_conn_neg` strongest negative correlations, using partial sort.
#'    Exclude zero-variance parcels from selection candidates.
#'    Construct a directed sparse graph `W_dir` from these.
#' 4. Symmetrize `W_dir` using `W_sym_raw = (W_dir + t(W_dir)) / 2`, then `drop0`.
#' 5. Ensure strict symmetry for z-scoring: `W_symmetric = Matrix::forceSymmetric(W_sym_raw, uplo="U")`.
#' 6. Z-score non-zero edge weights in `W_symmetric` (assumes `zscore_nonzero_sparse` is stable).
#'
#' @param X_subject A numeric matrix of time-series data for one subject
#'   (`T_i` time points x `V_p` parcels).
#' @param parcel_names A character vector of parcel names.
#' @param k_conn_pos An integer, number of positive connections to retain per node.
#' @param k_conn_neg An integer, number of negative connections to retain per node.
#' @param use_dtw Logical, defaults to `FALSE`. (Placeholder).
#' @details `X_subject` must contain at least two rows (time points).
#' @return A sparse symmetric `Matrix::dgCMatrix` of size `V_p x V_p`
#'   representing the z-scored connectivity graph `W_conn_i`.
#' @importFrom Matrix Matrix sparseMatrix drop0 t forceSymmetric
#' @importFrom stats cor sd
#' @keywords internal
compute_subject_connectivity_graph_sparse <- function(X_subject, parcel_names,
                                                      k_conn_pos, k_conn_neg,
                                                      use_dtw = FALSE) {
  V_p <- ncol(X_subject)

  if (length(parcel_names) != V_p) {
    stop(sprintf("`parcel_names` must have length %d to match number of columns in `X_subject`.", V_p))
  }
  if (!is.numeric(k_conn_pos) || length(k_conn_pos) != 1 || k_conn_pos < 0 || k_conn_pos != round(k_conn_pos)) {
    stop("`k_conn_pos` must be a non-negative integer.")
  }
  if (!is.numeric(k_conn_neg) || length(k_conn_neg) != 1 || k_conn_neg < 0 || k_conn_neg != round(k_conn_neg)) {
    stop("`k_conn_neg` must be a non-negative integer.")
  }

  if (V_p == 0) return(Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(character(0), character(0))))

  T_i <- nrow(X_subject)
  if (T_i < 2) {
    stop("`X_subject` must have at least two rows (time points).")
  }

  if (use_dtw && interactive()) {
    message("Note: DTW is not yet implemented. Proceeding with standard correlations.")
  }

  col_means <- colMeans(X_subject, na.rm = TRUE)
  col_sds <- apply(X_subject, 2, stats::sd, na.rm = TRUE)
  zero_var_indices <- which(col_sds == 0)
  if (length(zero_var_indices) > 0 && interactive()) {
    message(sprintf("Found %d parcel(s) with zero variance. These will be excluded from k-NN selection and their correlations are 0.", length(zero_var_indices)))
  }
  
  X_centered <- sweep(X_subject, 2, col_means, "-")

  row_indices_list <- vector("list", V_p)
  col_indices_list <- vector("list", V_p)
  values_list <- vector("list", V_p)
  
  # Create a full index of parcels to potentially select from, excluding zero-variance ones for candidates
  selectable_parcel_indices <- setdiff(1:V_p, zero_var_indices)

  if ((k_conn_pos > 0 || k_conn_neg > 0) && length(selectable_parcel_indices) > 0) {
    cor_blocks <- chunked_sparse_cor(X_centered, col_sds, block_size = 50L)
    for (block in cor_blocks) {
      block_indices <- block$indices
      block_corrs <- block$corr
      for (offset in seq_along(block_indices)) {
        i <- block_indices[offset]
        if (i %in% zero_var_indices) {
          row_indices_list[[i]] <- integer(0)
          col_indices_list[[i]] <- integer(0)
          values_list[[i]] <- numeric(0)
          next
        }

        node_corrs_full_row <- as.numeric(block_corrs[offset, ])
        node_corrs_full_row[i] <- 0
        if (length(zero_var_indices) > 0) {
          node_corrs_full_row[zero_var_indices] <- 0
        }

      node_corrs_candidates <- node_corrs_full_row[selectable_parcel_indices]
      
      current_selected_indices_in_selectable <- integer(0)
      current_selected_values <- numeric(0)
      
      if (k_conn_pos > 0) {
        pos_candidates_idx_in_subset <- which(node_corrs_candidates > 0)
        if (length(pos_candidates_idx_in_subset) > 0) {
          pos_subset_vals <- node_corrs_candidates[pos_candidates_idx_in_subset]
          num_to_keep_pos <- min(k_conn_pos, length(pos_subset_vals))
          # Use partial sorting to avoid ordering the full vector
          part_order <- order(
            pos_subset_vals,
            decreasing = TRUE,
            partial = seq_len(num_to_keep_pos)
          )[seq_len(num_to_keep_pos)]
          # Fully sort the selected top-k to maintain previous behaviour
          top_pos_in_subset_ordered_idx <- part_order[order(pos_subset_vals[part_order], decreasing = TRUE)]

          current_selected_indices_in_selectable <- c(
            current_selected_indices_in_selectable,
            pos_candidates_idx_in_subset[top_pos_in_subset_ordered_idx]
          )
          current_selected_values <- c(
            current_selected_values,
            pos_subset_vals[top_pos_in_subset_ordered_idx]
          )
        }
      }
      
      if (k_conn_neg > 0) {
        neg_candidates_idx_in_subset <- which(node_corrs_candidates < 0)
        if (length(neg_candidates_idx_in_subset) > 0) {
          neg_subset_vals <- node_corrs_candidates[neg_candidates_idx_in_subset]
          num_to_keep_neg <- min(k_conn_neg, length(neg_subset_vals))
          part_order <- order(
            neg_subset_vals,
            decreasing = FALSE,
            partial = seq_len(num_to_keep_neg)
          )[seq_len(num_to_keep_neg)]
          top_neg_in_subset_ordered_idx <- part_order[order(neg_subset_vals[part_order], decreasing = FALSE)]

          current_selected_indices_in_selectable <- c(
            current_selected_indices_in_selectable,
            neg_candidates_idx_in_subset[top_neg_in_subset_ordered_idx]
          )
          current_selected_values <- c(
            current_selected_values,
            neg_subset_vals[top_neg_in_subset_ordered_idx]
          )
        }
      }
      
      if(length(current_selected_indices_in_selectable) > 0) {
          # Map indices from subset back to original V_p indices
          final_selected_original_indices <- selectable_parcel_indices[current_selected_indices_in_selectable]
          row_indices_list[[i]] <- rep(i, length(final_selected_original_indices))
          col_indices_list[[i]] <- final_selected_original_indices
          values_list[[i]] <- current_selected_values
      } else {
          row_indices_list[[i]] <- integer(0) # Ensure it's not NULL
          col_indices_list[[i]] <- integer(0)
          values_list[[i]] <- numeric(0)

      }
    }
  }
}

  final_row_indices <- unlist(row_indices_list)
  final_col_indices <- unlist(col_indices_list)
  final_values <- unlist(values_list)

  if (length(final_row_indices) > 0) {
      W_dir <- Matrix::sparseMatrix(
        i = final_row_indices, j = final_col_indices, x = final_values,
        dims = c(V_p, V_p), dimnames = list(parcel_names, parcel_names)
      )
      W_dir <- Matrix::drop0(W_dir)
  } else { 
      W_dir <- Matrix::Matrix(0, nrow=V_p, ncol=V_p, sparse=TRUE, dimnames = list(parcel_names, parcel_names))
  }

  # Symmetrize W_dir 
  W_dir_t <- Matrix::t(W_dir)
  W_sum <- W_dir + W_dir_t # Sum of weights where edges exist from either direction

  # Create an empty Matrix of appropriate size
  W_symmetric <- Matrix::Matrix(0, nrow=V_p, ncol=V_p, sparse=TRUE, dimnames = list(parcel_names, parcel_names))
  
  # Only process non-zero entries if they exist
  nnz_W_sum <- tryCatch(Matrix::nnzero(W_sum), error = function(e) 0)
  if (!is.na(nnz_W_sum) && nnz_W_sum > 0) {
      # Create a denominator matrix W_den based on the non-zero pattern of W_sum
      idx <- Matrix::which(W_sum != 0, arr.ind = TRUE)
      den_values_at_idx <- pmax(1, as.numeric(W_dir[idx] != 0) + as.numeric(W_dir_t[idx] != 0))
      
      W_den_sparse <- Matrix::sparseMatrix(i = idx[,1], j = idx[,2], x = den_values_at_idx, 
                                         dims = dim(W_sum), dimnames = dimnames(W_sum))
      
      # Perform element-wise division for non-zero elements of W_sum
      W_symmetric_raw <- W_sum / W_den_sparse
      W_symmetric_raw[!is.finite(W_symmetric_raw)] <- 0
      
      # Ensure diagonal is zero
      diag(W_symmetric_raw) <- 0 
      W_symmetric_raw <- Matrix::drop0(W_symmetric_raw)
      
      # Make symmetric
      W_symmetric <- Matrix::forceSymmetric(W_symmetric_raw, uplo = "U")
      if (!is(W_symmetric, "sparseMatrix")) {
        W_symmetric <- Matrix::Matrix(W_symmetric, sparse = TRUE)
      }
  }
  
  # Ensure zero-variance parcels have zero connections
  if (length(zero_var_indices) > 0) {
      W_symmetric[zero_var_indices, ] <- 0
      W_symmetric[, zero_var_indices] <- 0
      W_symmetric <- Matrix::drop0(W_symmetric)
  }

  # Z-score non-zero elements
  W_conn_i <- zscore_nonzero_sparse(W_symmetric)
  
  return(as(W_conn_i, "generalMatrix"))
}

#' Compute sparse α-lazy random-walk normalized graph Laplacian `L = I - α D⁻¹ W`
#'
#' @param W_sparse A sparse, symmetric adjacency matrix (`Matrix::dgCMatrix`, `V_p x V_p`).
#'   Symmetry is verified using \code{Matrix::isSymmetric()} (tolerance \code{1e-8});
#'   the function stops with an error if the matrix is not symmetric.
#' @param alpha Numeric, the laziness parameter. Default is 0.93.
#'   Will be clamped to `[epsilon, 1]` range if outside `(0,1]`.
#' @param degree_type Character string, how to calculate node degrees if `W_sparse` has negative values.
#'   One of `"abs"` (default, sum of absolute weights), `"positive"` (sum of positive weights only),
#'   or `"signed"` (sum of raw weights). Documented for clarity.
#' @return A sparse, symmetric graph Laplacian matrix (`Matrix::dgCMatrix`, `V_p x V_p`).
#' @importFrom Matrix Diagonal rowSums t forceSymmetric isSymmetric
#' @keywords internal
compute_graph_laplacian_sparse <- function(W_sparse, alpha = 0.93, degree_type = "abs") {
  stopifnot(inherits(W_sparse, "Matrix"))
  if (!Matrix::isSymmetric(W_sparse, tol = 1e-8)) {
    stop("W_sparse must be symmetric. Symmetrize the adjacency matrix before calling compute_graph_laplacian_sparse().")
  }
  V_p <- nrow(W_sparse)
  if (V_p == 0) return(Matrix::Matrix(0,0,0,sparse=TRUE))
  
  # Audit suggestion: Clamp alpha
  if (alpha <= 0 || alpha > 1) {
      original_alpha <- alpha
      alpha <- min(max(alpha, .Machine$double.eps), 1)
      warning(sprintf("alpha value %.3f is outside typical (0,1] range. Clamped to %.3f.", original_alpha, alpha))
  }
  
  degree_type <- match.arg(degree_type, c("abs", "positive", "signed"))

  degree_vec <- switch(degree_type,
    abs = Matrix::rowSums(abs(W_sparse)),
    positive = Matrix::rowSums(W_sparse * (W_sparse > 0)), # or pmax(0, W_sparse)
    signed = Matrix::rowSums(W_sparse)
  )
  
  if (degree_type == "abs" && any(W_sparse@x < 0, na.rm=TRUE)) {
      # message("Note: Using 'abs' degree type with a graph containing negative weights. Interpret D^-1 W with care.")
  } else if (degree_type == "signed" && any(degree_vec <= 0, na.rm=TRUE)) {
      # message("Note: Using 'signed' degree type, and some node degrees are <=0. D^-1 will be 0 or undefined for these.")
  }
  
  zero_degree_indices <- which(degree_vec == 0)
  if (length(zero_degree_indices) > 0 && interactive()) {
      # message(sprintf("Found %d nodes with zero degree (type: %s). D^-1 will be zero for these.", length(zero_degree_indices), degree_type))
  }
  
  inv_degree_vec <- ifelse(degree_vec == 0, 0, 1 / degree_vec)
  D_inv_sparse <- Matrix::Diagonal(n = V_p, x = inv_degree_vec)
  
  L_rw_lazy <- Matrix::Diagonal(n=V_p) - alpha * (D_inv_sparse %*% W_sparse)
  
  # Instead of using forceSymmetric, use explicit symmetrization to match the test
  L_rw_lazy_sym <- (L_rw_lazy + Matrix::t(L_rw_lazy)) / 2
  L_rw_lazy_sym <- Matrix::drop0(L_rw_lazy_sym)
  
  if (!inherits(L_rw_lazy_sym, "sparseMatrix")) {
      L_rw_lazy_sym <- Matrix::Matrix(L_rw_lazy_sym, sparse = TRUE)
  }

  return(as(Matrix::drop0(L_rw_lazy_sym), "generalMatrix"))
}

#' Compute spectral sketch `U_orig_i` using `PRIMME` or `base::eigen`
#'
#' Computes the `k` eigenvectors of the sparse graph Laplacian `L_conn_i_sparse`
#' corresponding to the smallest, non-trivial eigenvalues.
#' Eigenvectors for eigenvalues numerically close to zero are discarded based on a dynamic tolerance.
#'
#' @param L_conn_i_sparse A sparse graph Laplacian matrix (`Matrix::dgCMatrix`, `V_p x V_p`).
#'   Must be symmetric.
#' @param k An integer, the desired spectral rank. Must be `k >= 0`.
#' @param eigenvalue_tol Numeric floor used when dynamically determining the
#'   tolerance for filtering near-zero eigenvalues. Defaults to `1e-8`.
#' @return A list containing two elements:
#'   - `vectors`: A dense matrix `U_orig_i` (`V_p x k_actual`) of eigenvectors.
#'     `k_actual` may be less than `k` if not enough informative eigenvectors are found.
#'   - `values`: A vector of eigenvalues corresponding to the eigenvectors.
#' @importFrom PRIMME eigs_sym
#' @keywords internal
compute_spectral_sketch_sparse <- function(L_conn_i_sparse, k,
                                           eigenvalue_tol = 1e-8) {
  V_p <- nrow(L_conn_i_sparse)
  # Use provided eigenvalue_tol as the floor for dynamic adjustment
  eigenvalue_tol_floor <- eigenvalue_tol

  if (k < 0) stop("`spectral_rank_k` (k) must be non-negative.")
  if (V_p == 0) return(list(vectors = matrix(0, 0, k), values = numeric(0)))
  if (k == 0) return(list(vectors = matrix(0, V_p, 0), values = numeric(0)))

  num_eigs_to_request <- min(V_p -1, k + 10) 
  if (V_p <= 1) num_eigs_to_request = 0

  use_dense_eigen <- FALSE
  if (V_p <= 5 || k >= V_p -1 ) { 
      use_dense_eigen <- TRUE
      if (interactive() && V_p > 5) message(sprintf("Note: k (%d) is high relative to V_p (%d) or V_p is small. Using base::eigen.", k, V_p))
  }
  if (V_p > 1 && num_eigs_to_request >= V_p) { 
      num_eigs_to_request = V_p -1 
  }
  if (num_eigs_to_request <=0 && V_p > 0) { 
      use_dense_eigen <- TRUE 
  }

  if (use_dense_eigen) {
    if (interactive() && V_p > 100) { 
        message(paste("Warning: Converting potentially large sparse Laplacian (", V_p, "x", V_p,
                      ") to dense for eigendecomposition. This may be memory intensive.", sep=""))
    }
    L_dense <- as.matrix(L_conn_i_sparse)
    if(any(!is.finite(L_dense))) {
        stop("Non-finite values in dense Laplacian before eigen decomposition.")
    }
    eigen_decomp <- eigen(L_dense, symmetric = TRUE)
    eigen_vals_raw <- eigen_decomp$values
    eigen_vecs_raw <- eigen_decomp$vectors
  } else {
    eigs_result <- PRIMME::eigs_sym(L_conn_i_sparse, 
                                      NEig = num_eigs_to_request, 
                                      which = "SM", 
                                      tol = 1e-9)
    eigen_vals_raw <- eigs_result$values
    eigen_vecs_raw <- eigs_result$vectors
  }
  
  sorted_order <- order(eigen_vals_raw)
  eigen_vals_sorted <- eigen_vals_raw[sorted_order]
  eigen_vecs_sorted <- eigen_vecs_raw[, sorted_order, drop = FALSE]
  
  # Audit suggestion: Dynamic eigenvalue tolerance
  if (length(eigen_vals_sorted) > 0) {
      median_abs_lambda <- stats::median(abs(eigen_vals_sorted), na.rm = TRUE)
      dynamic_tol_component <- 1e-4 * median_abs_lambda
      eigenvalue_tol <- max(eigenvalue_tol_floor, dynamic_tol_component, na.rm = TRUE) # Ensure NA doesn't break max if median_abs_lambda is NA
      if (is.na(eigenvalue_tol) || !is.finite(eigenvalue_tol)) eigenvalue_tol <- eigenvalue_tol_floor # Fallback
  } else {
      eigenvalue_tol <- eigenvalue_tol_floor
  }
  
  non_trivial_indices <- which(eigen_vals_sorted > eigenvalue_tol)
  
  U_informative <- eigen_vecs_sorted[, non_trivial_indices, drop = FALSE]
  Lambda_informative <- eigen_vals_sorted[non_trivial_indices]
  
  num_found_informative <- ncol(U_informative)
  
  if (num_found_informative == 0 && k > 0) {
      stop(sprintf("No informative eigenvectors found (all eigenvalues <= %.2e). Requested k=%d. Check graph structure or eigenvalue tolerance (current: %.2e).", eigenvalue_tol, k, eigenvalue_tol))
  }
  
  if (num_found_informative < k) {
    stop(sprintf("Found only %d informative eigenvectors (eigenvalues > %.2e), but k=%d was requested. Try reducing k or inspecting graph connectivity/Laplacian spectrum.", 
                    num_found_informative, eigenvalue_tol, k))
  } else {
    k_to_select <- k
  }
  
  U_orig_i <- U_informative[, 1:k_to_select, drop = FALSE]
  Lambda_orig_i <- Lambda_informative[1:k_to_select]
  
  return(list(vectors = U_orig_i, values = Lambda_orig_i))
}