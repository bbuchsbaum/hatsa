# Helper function to convert to general sparse matrix format
# This avoids the deprecated direct conversion to dgCMatrix
.to_general_sparse <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    return(mat)
  }
  return(as(as(mat, "generalMatrix"), "CsparseMatrix"))
}

#' Compute Task-Based Parcel Similarity Graph (W_task) from Activations
#'
#' Calculates a sparse, z-scored similarity graph between parcels based on their
#' activation profiles across different conditions or task features.
#'
#' @param activation_matrix A numeric matrix (`C x V_p`) where `C` is the number
#'   of conditions/features and `V_p` is the number of parcels. Each column
#'   represents the activation profile for a parcel.
#' @param parcel_names A character vector of length `V_p` specifying parcel names.
#' @param k_conn_task_pos Non-negative integer. Number of strongest positive connections to
#'   retain per parcel during sparsification.
#' @param k_conn_task_neg Non-negative integer. Number of strongest negative connections to
#'   retain per parcel during sparsification.
#' @param similarity_method Character string or function. Specifies the method to
#'   compute the initial `V_p x V_p` similarity matrix from `activation_matrix`.
#'   If "pearson" (default) or "spearman", `stats::cor` is used.
#'   If a function, it must take `activation_matrix` as input and return a
#'   `V_p x V_p` numeric matrix.
#' @param use_dtw Logical, defaults to `FALSE`. (Placeholder, currently unused but
#'   kept for potential future compatibility or signature consistency).
#'
#' @return A sparse, symmetric `Matrix::dgCMatrix` of size `V_p x V_p`
#'   representing the z-scored task-based similarity graph `W_task_i`.
#'
#' @export
#' @importFrom Matrix Matrix sparseMatrix drop0 t forceSymmetric
#' @importFrom stats cor sd
compute_W_task_from_activations <- function(activation_matrix,
                                            parcel_names,
                                            k_conn_task_pos,
                                            k_conn_task_neg,
                                            similarity_method = "pearson",
                                            use_dtw = FALSE) {

  V_p <- ncol(activation_matrix)
  C_n <- if (is.matrix(activation_matrix)) nrow(activation_matrix) else 0 # Number of conditions/rows

  if (!is.numeric(k_conn_task_pos) || length(k_conn_task_pos) != 1 ||
      k_conn_task_pos < 0 || k_conn_task_pos != round(k_conn_task_pos)) {
    stop("`k_conn_task_pos` must be a non-negative integer.")
  }
  if (!is.numeric(k_conn_task_neg) || length(k_conn_task_neg) != 1 ||
      k_conn_task_neg < 0 || k_conn_task_neg != round(k_conn_task_neg)) {
    stop("`k_conn_task_neg` must be a non-negative integer.")
  }

  if (V_p == 0) {
    # Ensure consistent return type (dgCMatrix) for empty graph
    empty_mat <- Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(character(0), character(0)))
    return(.to_general_sparse(empty_mat))
  }
  if (length(parcel_names) != V_p) {
    stop("Length of 'parcel_names' must match the number of columns (parcels) in 'activation_matrix'.")
  }
  
  # Handle C_n (number of conditions/rows) < 2, as cor() behaves unexpectedly or errors.
  if (C_n < 2) {
    warning(sprintf("compute_W_task_from_activations: activation_matrix has %d row(s) (conditions/features). Need at least 2 for meaningful correlation. Resulting graph will likely be empty or based on zero similarities.", C_n))
    # If C_n is 0 or 1, similarity matrix will be all 0s after NA replacement, leading to empty graph for kNN > 0.
    # So, can proceed, but the graph will be empty. Or return empty graph directly.
    # For safety and explicitness if C_n = 0 (which causes cor error):
    if (C_n == 0) {
        empty_mat_vp <- Matrix::Matrix(0, nrow=V_p, ncol=V_p, sparse=TRUE, dimnames = list(parcel_names, parcel_names))
        return(.to_general_sparse(empty_mat_vp))
    }
    # If C_n = 1, cor returns NAs, which are then set to 0. Resulting graph is empty.
    # This path will be handled by existing logic turning NAs to 0.
  }

  if (is.character(similarity_method) && similarity_method %in% c("pearson", "spearman")) {
    return(.build_sparse_graph_crossprod(
      activation_matrix, parcel_names,
      k_conn_task_pos, k_conn_task_neg,
      "Found %d parcel(s) with zero variance in activation profiles. Similarities involving these will be affected (likely 0 or NA)."
    ))
  } else if (is.function(similarity_method)) {
    sim_matrix_dense <- tryCatch({
      similarity_method(activation_matrix)
    }, error = function(e) {
      stop(paste("The provided similarity_method function failed:", e$message))
    })
    if (!is.matrix(sim_matrix_dense) || nrow(sim_matrix_dense) != V_p || ncol(sim_matrix_dense) != V_p) {
      stop("The custom similarity_method function must return a V_p x V_p matrix.")
    }
    sim_matrix_dense[is.na(sim_matrix_dense)] <- 0
    return(.sparsify_from_dense(sim_matrix_dense, parcel_names,
                                k_conn_task_pos, k_conn_task_neg))
  } else {
    stop("'similarity_method' must be 'pearson', 'spearman', or a function.")
  }
}

#' Compute Task-Based Parcel Similarity Graph (W_task) from Encoding Weights
#'
#' Calculates a sparse, z-scored similarity graph between parcels based on their
#' encoding weight profiles for a set of features.
#'
#' @param encoding_weights_matrix A numeric matrix (`V_p x N_features`) where `V_p`
#'   is the number of parcels and `N_features` is the number of encoding features.
#'   Each row represents the encoding weight profile for a parcel.
#' @param parcel_names A character vector of length `V_p` specifying parcel names.
#' @param k_conn_task_pos Non-negative integer. Number of strongest positive connections to
#'   retain per parcel during sparsification.
#' @param k_conn_task_neg Non-negative integer. Number of strongest negative connections to
#'   retain per parcel during sparsification.
#' @param similarity_method Character string or function. Specifies the method to
#'   compute the initial `V_p x V_p` similarity matrix.
#'   If "pearson" (default) or "spearman", `stats::cor` is used on the
#'   transposed input (to compare rows/parcels).
#'   If a function, it must take `encoding_weights_matrix` (V_p x N_features)
#'   as input and return a `V_p x V_p` numeric matrix.
#'
#' @return A sparse, symmetric `Matrix::dgCMatrix` of size `V_p x V_p`
#'   representing the z-scored task-based similarity graph `W_task_i`.
#'
#' @export
#' @importFrom Matrix Matrix sparseMatrix drop0 t forceSymmetric
#' @importFrom stats cor sd
compute_W_task_from_encoding <- function(encoding_weights_matrix,
                                         parcel_names,
                                         k_conn_task_pos,
                                         k_conn_task_neg,
                                         similarity_method = "pearson") {

  V_p <- nrow(encoding_weights_matrix) # Parcels are rows
  N_features <- ncol(encoding_weights_matrix)

  if (!is.numeric(k_conn_task_pos) || length(k_conn_task_pos) != 1 ||
      k_conn_task_pos < 0 || k_conn_task_pos != round(k_conn_task_pos)) {
    stop("`k_conn_task_pos` must be a non-negative integer.")
  }
  if (!is.numeric(k_conn_task_neg) || length(k_conn_task_neg) != 1 ||
      k_conn_task_neg < 0 || k_conn_task_neg != round(k_conn_task_neg)) {
    stop("`k_conn_task_neg` must be a non-negative integer.")
  }

  if (V_p == 0) {
    # Ensure consistent return type (dgCMatrix) for empty graph
    empty_mat <- Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(character(0), character(0)))
    return(.to_general_sparse(empty_mat))
  }
  if (length(parcel_names) != V_p) {
    stop("Length of 'parcel_names' must match the number of rows (parcels) in 'encoding_weights_matrix'.")
  }
  if (N_features == 0) {
     warning("compute_W_task_from_encoding: encoding_weights_matrix has zero features. Resulting graph will be empty.")
     # Ensure consistent return type (dgCMatrix) for empty graph
     empty_mat <- Matrix::Matrix(0, nrow=V_p, ncol=V_p, sparse=TRUE, dimnames = list(parcel_names, parcel_names))
     return(.to_general_sparse(empty_mat))
  }

  if (is.character(similarity_method) && similarity_method %in% c("pearson", "spearman")) {
    return(.build_sparse_graph_crossprod(
      t(encoding_weights_matrix), parcel_names,
      k_conn_task_pos, k_conn_task_neg,
      "Found %d parcel(s) with zero variance in encoding weights. Similarities involving these will be affected."
    ))
  } else if (is.function(similarity_method)) {
    sim_matrix_dense <- tryCatch({
      similarity_method(encoding_weights_matrix)
    }, error = function(e) {
      stop(paste("The provided similarity_method function failed:", e$message))
    })
    if (!is.matrix(sim_matrix_dense) || nrow(sim_matrix_dense) != V_p || ncol(sim_matrix_dense) != V_p) {
      stop("The custom similarity_method function must return a V_p x V_p matrix.")
    }
    sim_matrix_dense[is.na(sim_matrix_dense)] <- 0
    return(.sparsify_from_dense(sim_matrix_dense, parcel_names,
                                k_conn_task_pos, k_conn_task_neg))
  } else {
    stop("'similarity_method' must be 'pearson', 'spearman', or a function.")
  }
}

#' Compute Correlation Between Two Sparse Graphs
#'
#' Calculates the Spearman correlation between the edge weights of two sparse graphs,
#' considering the union of their non-zero edges and filling missing edges with 0.
#' Only the upper triangle of the matrices is considered.
#'
#' @param W_graph1 A sparse, symmetric matrix (`Matrix::dgCMatrix`, `V_p x V_p`).
#'   Assumed to have z-scored non-zero entries.
#' @param W_graph2 A sparse, symmetric matrix (`Matrix::dgCMatrix`, `V_p x V_p`).
#'   Assumed to have z-scored non-zero entries.
#' @param max_edges An integer or `Inf`. If the number of unique edges in the
#'   union of the upper triangles exceeds `max_edges`, a random sample of
#'   `max_edges` edges will be used for the correlation calculation. Defaults to 2,000,000.
#'
#' @return The Spearman correlation coefficient (`rho`) between the edge weights
#'   of the two graphs based on the union of their upper triangle edges. Returns `NA`
#'   if correlation cannot be computed (e.g., too few edges, zero variance).
#'
#' @export
#' @importFrom Matrix summary
#' @importFrom methods is
#' @importFrom stats cor sd
compute_graph_correlation <- function(W_graph1, W_graph2, max_edges = 2000000) {

  if (!is(W_graph1, "sparseMatrix") || !is(W_graph2, "sparseMatrix")) {
    stop("Inputs must be sparse matrices (e.g., dgCMatrix).")
  }
  if (nrow(W_graph1) != ncol(W_graph1) || nrow(W_graph2) != ncol(W_graph2) || nrow(W_graph1) != nrow(W_graph2)) {
    stop("Input graphs must be square matrices of the same dimension.")
  }

  # Extract triplets (i, j, x) from summaries
  summary1 <- Matrix::summary(W_graph1)
  summary2 <- Matrix::summary(W_graph2)

  # Handle case where summary returns empty or simple structure
  if (length(summary1) == 0 || !is.data.frame(summary1)) {
    edges1 <- data.frame(i = integer(0), j = integer(0), x = numeric(0))
  } else {
    # Filter for upper triangle (i < j)
    edges1 <- summary1[summary1$i < summary1$j, , drop = FALSE]
  }
  
  if (length(summary2) == 0 || !is.data.frame(summary2)) {
    edges2 <- data.frame(i = integer(0), j = integer(0), x = numeric(0))
  } else {
    # Filter for upper triangle (i < j)
    edges2 <- summary2[summary2$i < summary2$j, , drop = FALSE]
  }

  if (nrow(edges1) == 0 && nrow(edges2) == 0) {
    return(NA_real_)
  }

  # Vectorised union using matching on edge indices
  idx1 <- paste(edges1$i, edges1$j, sep = "-")
  idx2 <- paste(edges2$i, edges2$j, sep = "-")
  all_idx <- union(idx1, idx2)
  pos1 <- match(all_idx, idx1)
  pos2 <- match(all_idx, idx2)
  w1 <- ifelse(is.na(pos1), 0, edges1$x[pos1])
  w2 <- ifelse(is.na(pos2), 0, edges2$x[pos2])

  # Sample if needed
  num_edges <- length(all_idx)
  if (num_edges > max_edges && is.finite(max_edges) && max_edges > 0) {
    if (interactive()) {
      message(sprintf("Sampling %d edges from union of %d for correlation calculation.", max_edges, num_edges))
    }
    sample_indices <- sample.int(num_edges, size = max_edges, replace = FALSE)
    w1 <- w1[sample_indices]
    w2 <- w2[sample_indices]
    num_edges <- max_edges
  }

  if (num_edges < 2) {
    warning("Too few edges (< 2) in the sampled union to compute correlation. Returning NA.")
    return(NA_real_)
  }

  # Check for zero variance before attempting correlation
  sd1 <- stats::sd(w1)
  sd2 <- stats::sd(w2)

  if (is.na(sd1) || sd1 == 0 || is.na(sd2) || sd2 == 0) {
    warning("Zero variance in edge weights for at least one graph in the sampled union set. Correlation is undefined (NA).")
  }

  rho <- stats::cor(w1, w2, method = "spearman", use = "complete.obs")

  return(rho)
}

# Internal helper: build sparse correlation graph using crossprod and partial sorting
.build_sparse_graph_crossprod <- function(data_matrix, parcel_names,
                                          k_pos, k_neg, zero_var_message) {
  V_p <- ncol(data_matrix)
  N_obs <- nrow(data_matrix)

  col_means <- colMeans(data_matrix, na.rm = TRUE)
  col_sds <- apply(data_matrix, 2, stats::sd, na.rm = TRUE)
  zero_var_indices <- which(col_sds == 0)
  if (length(zero_var_indices) > 0) {
    message(sprintf(zero_var_message, length(zero_var_indices)))
  }

  data_centered <- sweep(data_matrix, 2, col_means, "-")

  row_indices_list <- vector("list", V_p)
  col_indices_list <- vector("list", V_p)
  values_list <- vector("list", V_p)

  selectable_parcel_indices <- setdiff(seq_len(V_p), zero_var_indices)

  denom <- max(1, N_obs - 1)

  if ((k_pos > 0 || k_neg > 0) && length(selectable_parcel_indices) > 0) {
    for (i in seq_len(V_p)) {
      if (i %in% zero_var_indices) {
        row_indices_list[[i]] <- integer(0)
        col_indices_list[[i]] <- integer(0)
        values_list[[i]] <- numeric(0)
        next
      }

      cov_vec <- as.numeric(crossprod(data_centered, data_centered[, i])) / denom
      node_corrs_full <- cov_vec / (col_sds * col_sds[i])
      node_corrs_full[!is.finite(node_corrs_full)] <- 0
      node_corrs_full[i] <- 0
      if (length(zero_var_indices) > 0) node_corrs_full[zero_var_indices] <- 0

      node_corrs_candidates <- node_corrs_full[selectable_parcel_indices]

      sel_idx_subset <- integer(0)
      sel_vals <- numeric(0)

      if (k_pos > 0) {
        pos_idx <- which(node_corrs_candidates > 0)
        if (length(pos_idx) > 0) {
          pos_vals <- node_corrs_candidates[pos_idx]
          num_keep_pos <- min(k_pos, length(pos_vals))
          ord <- order(pos_vals, decreasing = TRUE)
          ord <- head(ord, num_keep_pos)
          sel_idx_subset <- c(sel_idx_subset, pos_idx[ord])
          sel_vals <- c(sel_vals, pos_vals[ord])
        }
      }

      if (k_neg > 0) {
        neg_idx <- which(node_corrs_candidates < 0)
        if (length(neg_idx) > 0) {
          neg_vals <- node_corrs_candidates[neg_idx]
          num_keep_neg <- min(k_neg, length(neg_vals))
          ord <- order(neg_vals, decreasing = FALSE)
          ord <- head(ord, num_keep_neg)
          sel_idx_subset <- c(sel_idx_subset, neg_idx[ord])
          sel_vals <- c(sel_vals, neg_vals[ord])
        }
      }

      if (length(sel_idx_subset) > 0) {
        final_idx <- selectable_parcel_indices[sel_idx_subset]
        row_indices_list[[i]] <- rep(i, length(final_idx))
        col_indices_list[[i]] <- final_idx
        values_list[[i]] <- sel_vals
      } else {
        row_indices_list[[i]] <- integer(0)
        col_indices_list[[i]] <- integer(0)
        values_list[[i]] <- numeric(0)
      }
    }
  }

  final_row_indices <- unlist(row_indices_list)
  final_col_indices <- unlist(col_indices_list)
  final_values <- unlist(values_list)

  W_dir <- if (length(final_row_indices) > 0) {
    Matrix::sparseMatrix(
      i = final_row_indices, j = final_col_indices, x = final_values,
      dims = c(V_p, V_p), dimnames = list(parcel_names, parcel_names)
    )
  } else {
    Matrix::Matrix(0, nrow = V_p, ncol = V_p, sparse = TRUE,
                  dimnames = list(parcel_names, parcel_names))
  }
  W_dir <- Matrix::drop0(W_dir)

  W_dir_t <- Matrix::t(W_dir)
  W_sum <- W_dir + W_dir_t

  nnz_sum <- tryCatch(Matrix::nnzero(W_sum), error = function(e) 0)
  if (!is.na(nnz_sum) && nnz_sum > 0) {
    idx <- Matrix::which(W_sum != 0, arr.ind = TRUE)
    den_vals <- pmax(1, as.numeric(W_dir[idx] != 0) + as.numeric(W_dir_t[idx] != 0))
    W_den <- Matrix::sparseMatrix(i = idx[,1], j = idx[,2], x = den_vals,
                                  dims = dim(W_sum), dimnames = dimnames(W_sum))
    W_sym_raw <- W_sum / W_den
    W_sym_raw[!is.finite(W_sym_raw)] <- 0
    diag(W_sym_raw) <- 0
    W_sym_raw <- Matrix::drop0(W_sym_raw)
    W_sym <- Matrix::forceSymmetric(W_sym_raw, uplo = "U")
    if (!is(W_sym, "sparseMatrix")) {
      W_sym <- Matrix::Matrix(W_sym, sparse = TRUE)
    }
  } else {
    W_sym <- Matrix::Matrix(0, nrow = V_p, ncol = V_p, sparse = TRUE,
                            dimnames = list(parcel_names, parcel_names))
  }

  if (length(zero_var_indices) > 0) {
    W_sym[zero_var_indices, ] <- 0
    W_sym[, zero_var_indices] <- 0
    W_sym <- Matrix::drop0(W_sym)
  }

  W_final <- zscore_nonzero_sparse(W_sym)
  return(.to_general_sparse(W_final))
}

# Internal helper: sparsify a dense similarity matrix and z-score edges
.sparsify_from_dense <- function(sim_matrix_dense, parcel_names, k_pos, k_neg) {
  V_p <- nrow(sim_matrix_dense)
  diag(sim_matrix_dense) <- 0

  row_indices_list <- vector("list", V_p)
  col_indices_list <- vector("list", V_p)
  values_list <- vector("list", V_p)

  if (k_pos > 0 || k_neg > 0) {
    for (i in seq_len(V_p)) {
      node_vals <- sim_matrix_dense[i, ]
      sel_idx <- integer(0)
      sel_vals <- numeric(0)

      if (k_pos > 0) {
        pos_idx <- which(node_vals > 1e-9)
        if (length(pos_idx) > 0) {
          pos_vals <- node_vals[pos_idx]
          num_keep <- min(k_pos, length(pos_vals))
          ord <- order(pos_vals, decreasing = TRUE)
          ord <- head(ord, num_keep)
          sel_idx <- c(sel_idx, pos_idx[ord])
          sel_vals <- c(sel_vals, pos_vals[ord])
        }
      }

      if (k_neg > 0) {
        neg_idx <- which(node_vals < -1e-9)
        if (length(neg_idx) > 0) {
          neg_vals <- node_vals[neg_idx]
          num_keep <- min(k_neg, length(neg_vals))
          ord <- order(neg_vals, decreasing = FALSE)
          ord <- head(ord, num_keep)
          sel_idx <- c(sel_idx, neg_idx[ord])
          sel_vals <- c(sel_vals, neg_vals[ord])
        }
      }

      if (length(sel_idx) > 0) {
        row_indices_list[[i]] <- rep(i, length(sel_idx))
        col_indices_list[[i]] <- sel_idx
        values_list[[i]] <- sel_vals
      }
    }
  }

  final_row_indices <- unlist(row_indices_list)
  final_col_indices <- unlist(col_indices_list)
  final_values <- unlist(values_list)

  W_dir <- if (length(final_row_indices) > 0) {
    Matrix::sparseMatrix(
      i = final_row_indices, j = final_col_indices, x = final_values,
      dims = c(V_p, V_p), dimnames = list(parcel_names, parcel_names)
    )
  } else {
    Matrix::Matrix(0, nrow = V_p, ncol = V_p, sparse = TRUE,
                  dimnames = list(parcel_names, parcel_names))
  }
  W_dir <- Matrix::drop0(W_dir)

  W_dir_t <- Matrix::t(W_dir)
  W_sum <- W_dir + W_dir_t
  W_den_val <- as.numeric((W_dir != 0) + (W_dir_t != 0))
  W_den <- Matrix::Matrix(pmax(1, W_den_val), nrow = V_p, ncol = V_p, sparse = TRUE)
  W_sym_raw <- W_sum / W_den
  W_sym_raw[is.nan(W_sym_raw)] <- 0
  W_sym_raw[is.infinite(W_sym_raw)] <- 0
  W_sym_raw <- Matrix::drop0(W_sym_raw)
  W_sym <- Matrix::forceSymmetric(W_sym_raw, uplo = "U")

  W_final <- zscore_nonzero_sparse(W_sym)
  return(.to_general_sparse(W_final))
}

# Internal helper to sparsify and z-score a symmetric matrix
# Takes a symmetric matrix (potentially dense) and applies k-NN
# sparsification using the same k for positive and negative edges.
# Then z-scores the non-zero elements.
.sparsify_symmetric_matrix <- function(input_matrix, k_nn, parcel_names) {
  V_p <- nrow(input_matrix)
  if (V_p == 0) {
    return(Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(character(0), character(0))))
  }

  # Work with a sparse matrix directly to avoid unnecessary densification
  input_matrix <- .to_general_sparse(input_matrix)
  diag(input_matrix) <- 0 # Ensure diagonal is zero

  row_indices_list <- vector("list", V_p)
  col_indices_list <- vector("list", V_p)
  values_list <- vector("list", V_p)

  if (k_nn > 0) {
    for (i in 1:V_p) {
      node_vals <- as.numeric(input_matrix[i, ])
      selected_indices <- integer(0)
      selected_values <- numeric(0)

      # Keep k largest positive values
      pos_idx <- which(node_vals > 1e-9)
      if (length(pos_idx) > 0) {
        pos_vals <- node_vals[pos_idx]
        num_keep <- min(k_nn, length(pos_vals))
        ordered_idx <- order(pos_vals, decreasing = TRUE)
        top_idx_local <- head(ordered_idx, num_keep)
        selected_indices <- c(selected_indices, pos_idx[top_idx_local])
        selected_values <- c(selected_values, pos_vals[top_idx_local])
      }

      # Keep k most negative values (smallest values)
      neg_idx <- which(node_vals < -1e-9)
      if (length(neg_idx) > 0) {
        neg_vals <- node_vals[neg_idx]
        num_keep <- min(k_nn, length(neg_vals))
        ordered_idx <- order(neg_vals, decreasing = FALSE)
        top_idx_local <- head(ordered_idx, num_keep)
        selected_indices <- c(selected_indices, neg_idx[top_idx_local])
        selected_values <- c(selected_values, neg_vals[top_idx_local])
      }

      if(length(selected_indices) > 0) {
        row_indices_list[[i]] <- rep(i, length(selected_indices))
        col_indices_list[[i]] <- selected_indices
        values_list[[i]] <- selected_values
      }
    }
  }

  final_row_indices <- unlist(row_indices_list[!sapply(row_indices_list, is.null)])
  final_col_indices <- unlist(col_indices_list[!sapply(col_indices_list, is.null)])
  final_values <- unlist(values_list[!sapply(values_list, is.null)])

  W_dir <- if (length(final_row_indices) > 0) {
    Matrix::sparseMatrix(
      i = final_row_indices, j = final_col_indices, x = final_values,
      dims = c(V_p, V_p), dimnames = list(parcel_names, parcel_names)
    )
  } else {
    Matrix::Matrix(0, nrow=V_p, ncol=V_p, sparse=TRUE, dimnames = list(parcel_names, parcel_names))
  }
  W_dir <- Matrix::drop0(W_dir)

  # Symmetrize
  W_dir_t <- Matrix::t(W_dir)
  W_sum <- W_dir + W_dir_t
  W_den_val <- as.numeric((W_dir != 0) + (W_dir_t != 0))
  W_den <- Matrix::Matrix(pmax(1, W_den_val), nrow=V_p, ncol=V_p, sparse=TRUE)
  W_sym_raw <- W_sum / W_den
  if (inherits(W_sym_raw, "sparseMatrix")) {
    if (any(is.nan(W_sym_raw@x))) W_sym_raw@x[is.nan(W_sym_raw@x)] <- 0
    if (any(is.infinite(W_sym_raw@x))) W_sym_raw@x[is.infinite(W_sym_raw@x)] <- 0
  } else {
    W_sym_raw[is.nan(W_sym_raw)] <- 0
    W_sym_raw[is.infinite(W_sym_raw)] <- 0
  }
  W_sym_raw <- Matrix::drop0(W_sym_raw)
  W_sym <- Matrix::forceSymmetric(W_sym_raw, uplo = "U")

  # Z-score
  W_sym <- zscore_nonzero_sparse(W_sym)
  W_final <- Matrix::drop0(W_sym)

  return(.to_general_sparse(W_final))
}

#' Residualize Graph B based on Subspace from Graph A's Laplacian
#'
#' Projects `W_graph_to_residualize` onto the subspace spanned by the first
#' `k_eigenvectors_to_remove` smallest eigenvectors of
#' `L_graph_for_projection` and subtracts this projection.
#' The residual graph is then re-sparsified and re-z-scored.
#'
#' @param W_graph_to_residualize The sparse graph matrix (`dgCMatrix`) to be
#'   residualized (e.g., `W_task`).
#' @param L_graph_for_projection The sparse Laplacian matrix (`dgCMatrix`) from
#'   which the projection subspace is derived (e.g., `L_conn`). Must be symmetric.
#' @param k_eigenvectors_to_remove Integer, the number of smallest (by magnitude)
#'   eigenvectors of `L_graph_for_projection` to define the subspace for projection.
#'   Defaults to 64.
#' @param k_nn_resparsify Integer, the k value for k-NN sparsification applied
#'   to the residual graph (same k used for positive and negative edges).
#' @param eigenvalue_tol Numeric, tolerance for eigenvalue decomposition convergence
#'   and for identifying near-zero eigenvalues if needed (though projection uses the space).
#'   Default 1e-9.
#'
#' @return A sparse, symmetric, z-scored `dgCMatrix` representing the residualized graph.
#'
#' @importFrom Matrix t crossprod tcrossprod forceSymmetric
#' @importFrom PRIMME eigs_sym
#' @importFrom methods is
#' @keywords internal
residualize_graph_on_subspace <- function(W_graph_to_residualize,
                                          L_graph_for_projection,
                                          k_eigenvectors_to_remove = 64,
                                          k_nn_resparsify,
                                          eigenvalue_tol = 1e-9) {

  if (!is(W_graph_to_residualize, "sparseMatrix")) stop("W_graph_to_residualize must be a sparse Matrix.")
  if (!is(L_graph_for_projection, "sparseMatrix")) stop("L_graph_for_projection must be a sparse Matrix.")
  V_p <- nrow(W_graph_to_residualize)
  if (V_p == 0) return(W_graph_to_residualize) # Return empty graph if input is empty
  if (nrow(L_graph_for_projection) != V_p || ncol(L_graph_for_projection) != V_p || ncol(W_graph_to_residualize) != V_p) {
    stop("Input matrices must be square and of the same dimension.")
  }
  if (missing(k_nn_resparsify) || !is.numeric(k_nn_resparsify) || k_nn_resparsify <= 0) {
     stop("`k_nn_resparsify` must be a positive integer.")
  }

  k_proj <- min(k_eigenvectors_to_remove, V_p - 1) # Cannot request more than V_p-1
  if (k_proj <= 0) {
    warning("k_eigenvectors_to_remove is too small or V_p <= 1. Returning original graph.")
    return(W_graph_to_residualize)
  }

  # 1. Get eigenvectors U of L_graph_for_projection (L_A)
  message_stage(sprintf("Computing %d eigenvectors of projection Laplacian...", k_proj), interactive_only = TRUE)
  eigs_result <- tryCatch({
    PRIMME::eigs_sym(L_graph_for_projection,
                       NEig = k_proj,
                       which = "SM", # Smallest Magnitude
                       tol = eigenvalue_tol)
  }, error = function(e) {
    warning(paste("Eigen decomposition failed:", e$message, ". Returning original graph."))
    return(NULL)
  })

  if (is.null(eigs_result) || is.null(eigs_result$vectors) || ncol(eigs_result$vectors) < k_proj) {
      warning(sprintf("Eigen decomposition did not return the requested %d vectors. Returning original graph.", k_proj))
      return(W_graph_to_residualize)
  }
  U <- eigs_result$vectors # V_p x k_proj dense matrix
  
  # Orthonormalize U (important for projection formula)
  # Using SVD: U_ortho = svd(U)$u is robust
  # Note: svd(U)$u will have min(V_p, k_proj) columns. Should be k_proj if k_proj < V_p.
  U_ortho <- svd(U, nu = k_proj, nv = 0)$u
  if (ncol(U_ortho) != k_proj) {
      warning("Orthonormalization of eigenvectors failed to produce correct dimensions. Returning original graph.")
      return(W_graph_to_residualize)
  }

  # 2. Compute residual using low-rank formula
  # W_res = W_B - P W_B - W_B P + P W_B P, where P = U U^T
  # W_res = W_B - U(U^T W_B) - (W_B U)U^T + U(U^T W_B U)U^T
  # Let W_B be the graph to residualize
  message_stage("Computing residual graph projection...", interactive_only = TRUE)
  W_B <- W_graph_to_residualize
  U <- U_ortho # Use orthonormalized version

  dense_mem_gb <- (V_p * V_p * 8) / (1024^3)
  if (dense_mem_gb > 0.5) {
    warning(sprintf("Residualization will create ~%.2f GB dense matrices; ensure sufficient memory.", dense_mem_gb))
  }

  # Calculate intermediate terms (sparse-dense and dense-dense products)
  # UT_WB = t(U) %*% W_B  (k_proj x V_p, result is dense)
  UT_WB <- Matrix::crossprod(U, W_B)
  # WB_U = W_B %*% U      (V_p x k_proj, result is dense)
  WB_U <- W_B %*% U
  # UT_WB_U = t(U) %*% W_B %*% U = UT_WB %*% U (k_proj x k_proj, dense)
  UT_WB_U <- UT_WB %*% U

  # Calculate terms of the residual formula
  # Term1 = U %*% UT_WB   (V_p x V_p, dense)
  Term1 <- U %*% UT_WB
  # Term2 = WB_U %*% t(U) (V_p x V_p, dense)
  Term2 <- Matrix::tcrossprod(WB_U, U)
  # Term3 = U %*% UT_WB_U %*% t(U) (V_p x V_p, dense)
  Term3 <- U %*% UT_WB_U %*% Matrix::t(U)
  
  # Compute residual W_res (dense). Sparse input will be promoted automatically.
  W_res_dense <- -Term1 - Term2 + Term3
  W_res_dense <- W_res_dense + W_B
  message_stage("Symmetrizing residual graph...", interactive_only = TRUE)

  # 3. Symmetrize W_res
  W_res_sym <- (W_res_dense + t(W_res_dense)) / 2
  
  # Get parcel names for the helper function
  pnames <- rownames(W_graph_to_residualize) # Assume dimnames exist and match
  if (is.null(pnames)) pnames <- paste0("P", 1:V_p)

  # 4. Re-sparsify and Re-z-score using the helper
  message_stage(sprintf("Re-sparsifying residual graph (k=%d)...", k_nn_resparsify), interactive_only = TRUE)
  W_res_final <- .sparsify_symmetric_matrix(W_res_sym, k_nn = k_nn_resparsify, parcel_names = pnames)

  message_stage("Residualization complete.", interactive_only = TRUE)
  return(W_res_final)
}

#' Blend Connectivity and Task Laplacians
#'
#' Combines a connectivity-based Laplacian and a task-based Laplacian using
#' a specified blending method. Currently, only linear blending is supported.
#'
#' @param L_conn A sparse Laplacian matrix derived from connectivity (`dgCMatrix`, `V_p x V_p`).
#' @param L_task A sparse Laplacian matrix derived from task activations/encodings
#'   (`dgCMatrix`, `V_p x V_p`). Must have the same dimensions as `L_conn`.
#' @param lambda_blend_value Numeric, the blending factor (`λ`). Must be between 0 and 1.
#'   `L_hybrid = (1 - λ) * L_conn + λ * L_task`.
#' @param method Character string, the blending method. Currently only "linear" is supported.
#'
#' @return The blended sparse Laplacian matrix `L_hybrid` (`dgCMatrix`, `V_p x V_p`).
#'
#' @importFrom methods is
#' @keywords internal
blend_laplacians <- function(L_conn, L_task, lambda_blend_value, method = "linear") {

  if (!is(L_conn, "sparseMatrix") || !is(L_task, "sparseMatrix")) {
    stop("Inputs L_conn and L_task must be sparse matrices.")
  }
  if (nrow(L_conn) != ncol(L_conn) || nrow(L_task) != ncol(L_task) || nrow(L_conn) != nrow(L_task)) {
    stop("Input Laplacians must be square matrices of the same dimension.")
  }
  if (!is.numeric(lambda_blend_value) || length(lambda_blend_value) != 1 || lambda_blend_value < 0 || lambda_blend_value > 1) {
    stop("lambda_blend_value must be a single numeric value between 0 and 1.")
  }

  if (tolower(method) == "linear") {
    if (lambda_blend_value == 0) return(L_conn)
    if (lambda_blend_value == 1) return(L_task)

    # Sparse matrix arithmetic handles the weighted sum efficiently
    L_hybrid <- (1 - lambda_blend_value) * L_conn + lambda_blend_value * L_task

    # Ensure output is dgCMatrix
    return(.to_general_sparse(L_hybrid))

  } else {
    # Placeholder for future methods like "geo"
    stop(sprintf("Blending method '%s' is not currently supported. Only 'linear' is implemented.", method))
  }
} 