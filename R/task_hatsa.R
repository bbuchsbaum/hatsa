#' Compute Task-Based Parcel Similarity Graph (W_task) from Activations
#'
#' Calculates a sparse, z-scored similarity graph between parcels based on their
#' activation profiles across different conditions or task features.
#'
#' @param activation_matrix A numeric matrix (`C x V_p`) where `C` is the number
#'   of conditions/features and `V_p` is the number of parcels. Each column
#'   represents the activation profile for a parcel.
#' @param parcel_names A character vector of length `V_p` specifying parcel names.
#' @param k_conn_task_pos Integer, number of strongest positive connections to
#'   retain per parcel during sparsification.
#' @param k_conn_task_neg Integer, number of strongest negative connections to
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
#' @importFrom Matrix Matrix sparseMatrix drop0 t forceSymmetric
#' @importFrom stats cor sd
#' @keywords internal
#compute_W_task_from_activations <- function(activation_matrix,
#                                            parcel_names,
#                                            k_conn_task_pos,
#                                            k_conn_task_neg,
#                                            similarity_method = "pearson",
#                                            use_dtw = FALSE) {
#
#  V_p <- ncol(activation_matrix)
#  if (V_p == 0) {
#    return(Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(character(0), character(0))))
#  }
#  if (length(parcel_names) != V_p) {
#    stop("Length of 'parcel_names' must match the number of columns (parcels) in 'activation_matrix'.")
#  }
#
#  # 1. Compute V_p x V_p dense similarity matrix
#  sim_matrix_dense <- NULL
#  if (is.character(similarity_method) && similarity_method %in% c("pearson", "spearman")) {
#    if (nrow(activation_matrix) > 1) {
#        col_sds <- apply(activation_matrix, 2, stats::sd, na.rm = TRUE)
#        if (any(col_sds == 0, na.rm = TRUE) && interactive()) {
#            message(sprintf("compute_W_task_from_activations: Found %d parcel(s) with zero variance in activation profiles. Similarities involving these will be affected (likely 0 or NA).", sum(col_sds==0, na.rm=TRUE)))
#        }
#    }
#    sim_matrix_dense <- stats::cor(activation_matrix, method = similarity_method, use = "pairwise.complete.obs")
#  } else if (is.function(similarity_method)) {
#    sim_matrix_dense <- tryCatch({
#      similarity_method(activation_matrix)
#    }, error = function(e) {
#      stop(paste("The provided similarity_method function failed:", e$message))
#    })
#    if (!is.matrix(sim_matrix_dense) || nrow(sim_matrix_dense) != V_p || ncol(sim_matrix_dense) != V_p) {
#      stop("The custom similarity_method function must return a V_p x V_p matrix.")
#    }
#  } else {
#    stop("'similarity_method' must be 'pearson', 'spearman', or a function.")
#  }
#
#  sim_matrix_dense[is.na(sim_matrix_dense)] <- 0
#  diag(sim_matrix_dense) <- 0
#
#  # 2. Sparsify
#  row_indices_list <- vector("list", V_p)
#  col_indices_list <- vector("list", V_p)
#  values_list <- vector("list", V_p)
#
#  if (k_conn_task_pos > 0 || k_conn_task_neg > 0) {
#    for (i in 1:V_p) {
#      node_sims <- sim_matrix_dense[i, ]
#      current_selected_indices <- integer(0)
#      current_selected_values <- numeric(0)
#
#      if (k_conn_task_pos > 0) {
#        pos_candidates_idx <- which(node_sims > 1e-9) 
#        if (length(pos_candidates_idx) > 0) {
#          pos_candidates_vals <- node_sims[pos_candidates_idx]
#          num_to_keep_pos <- min(k_conn_task_pos, length(pos_candidates_vals))
#          ordered_in_pos_group <- order(pos_candidates_vals, decreasing = TRUE)
#          top_pos_in_group_idx <- head(ordered_in_pos_group, num_to_keep_pos)
#
#          current_selected_indices <- c(current_selected_indices, pos_candidates_idx[top_pos_in_group_idx])
#          current_selected_values <- c(current_selected_values, pos_candidates_vals[top_pos_in_group_idx])
#        }
#      }
#
#      if (k_conn_task_neg > 0) {
#        neg_candidates_idx <- which(node_sims < -1e-9) 
#        if (length(neg_candidates_idx) > 0) {
#          neg_candidates_vals <- node_sims[neg_candidates_idx]
#          num_to_keep_neg <- min(k_conn_task_neg, length(neg_candidates_vals))
#          ordered_in_neg_group <- order(neg_candidates_vals, decreasing = FALSE)
#          top_neg_in_group_idx <- head(ordered_in_neg_group, num_to_keep_neg)
#
#          current_selected_indices <- c(current_selected_indices, neg_candidates_idx[top_neg_in_group_idx])
#          current_selected_values <- c(current_selected_values, neg_candidates_vals[top_neg_in_group_idx])
#        }
#      }
#
#      if(length(current_selected_indices) > 0) {
#          row_indices_list[[i]] <- rep(i, length(current_selected_indices))
#          col_indices_list[[i]] <- current_selected_indices
#          values_list[[i]] <- current_selected_values
#      }
#    }
#  }
#
#  final_row_indices <- unlist(row_indices_list[!sapply(row_indices_list, is.null)])
#  final_col_indices <- unlist(col_indices_list[!sapply(col_indices_list, is.null)])
#  final_values <- unlist(values_list[!sapply(values_list, is.null)])
#
#  W_dir_task <- if (length(final_row_indices) > 0) {
#    Matrix::sparseMatrix(
#      i = final_row_indices, j = final_col_indices, x = final_values,
#      dims = c(V_p, V_p), dimnames = list(parcel_names, parcel_names)
#    )
#  } else {
#    Matrix::Matrix(0, nrow=V_p, ncol=V_p, sparse=TRUE, dimnames = list(parcel_names, parcel_names))
#  }
#  W_dir_task <- Matrix::drop0(W_dir_task)
#
#  # 3. Symmetrize
#  W_dir_task_t <- Matrix::t(W_dir_task)
#  W_sum_task <- W_dir_task + W_dir_task_t
#  
#  W_den_task_val <- as.numeric((W_dir_task != 0) + (W_dir_task_t != 0))
#  W_den_task <- Matrix::Matrix(pmax(1, W_den_task_val), nrow=V_p, ncol=V_p, sparse=TRUE)
#
#  W_symmetric_raw_task <- W_sum_task / W_den_task
#  if (inherits(W_symmetric_raw_task, "sparseMatrix")) {
#    if (any(is.nan(W_symmetric_raw_task@x))) W_symmetric_raw_task@x[is.nan(W_symmetric_raw_task@x)] <- 0
#    if (any(is.infinite(W_symmetric_raw_task@x))) W_symmetric_raw_task@x[is.infinite(W_symmetric_raw_task@x)] <- 0
#  } else {
#    W_symmetric_raw_task[is.nan(W_symmetric_raw_task)] <- 0
#    W_symmetric_raw_task[is.infinite(W_symmetric_raw_task)] <- 0
#  }
#  W_symmetric_raw_task <- Matrix::drop0(W_symmetric_raw_task)
#  W_symmetric_task <- Matrix::forceSymmetric(W_symmetric_raw_task, uplo = "U")
#
#  # 4. Z-score non-zero edge weights
#  if (length(W_symmetric_task@x) > 0) { 
#    non_zero_vals <- W_symmetric_task@x
#    mean_val <- mean(non_zero_vals)
#    sd_val <- stats::sd(non_zero_vals)
#    if (is.na(sd_val) || sd_val == 0) { 
#      W_symmetric_task@x <- rep(0, length(non_zero_vals))
#    } else {
#      W_symmetric_task@x <- (non_zero_vals - mean_val) / sd_val
#    }
#    W_task_i <- Matrix::drop0(W_symmetric_task)
#  } else {
#    W_task_i <- W_symmetric_task 
#  }
#  
#  return(as(W_task_i, "dgCMatrix"))
#}

#' Compute Task-Based Parcel Similarity Graph (W_task) from Encoding Weights
#'
#' Calculates a sparse, z-scored similarity graph between parcels based on their
#' activation profiles across different conditions or task features.
#'
#' @param activation_matrix A numeric matrix (`C x V_p`) where `C` is the number
#'   of conditions/features and `V_p` is the number of parcels. Each column
#'   represents the activation profile for a parcel.
#' @param parcel_names A character vector of length `V_p` specifying parcel names.
#' @param k_conn_task_pos Integer, number of strongest positive connections to
#'   retain per parcel during sparsification.
#' @param k_conn_task_neg Integer, number of strongest negative connections to
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
#' @importFrom Matrix Matrix sparseMatrix drop0 t forceSymmetric
#' @importFrom stats cor sd
#' @keywords internal
#compute_W_task_from_activations <- function(activation_matrix,
#                                            parcel_names,
#                                            k_conn_task_pos,
#                                            k_conn_task_neg,
#                                            similarity_method = "pearson",
#                                            use_dtw = FALSE) {
#
#  V_p <- ncol(activation_matrix)
#  if (V_p == 0) {
#    return(Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(character(0), character(0))))
#  }
#  if (length(parcel_names) != V_p) {
#    stop("Length of 'parcel_names' must match the number of columns (parcels) in 'activation_matrix'.")
#  }
#
#  # 1. Compute V_p x V_p dense similarity matrix
#  sim_matrix_dense <- NULL
#  if (is.character(similarity_method) && similarity_method %in% c("pearson", "spearman")) {
#    if (nrow(activation_matrix) > 1) {
#        col_sds <- apply(activation_matrix, 2, stats::sd, na.rm = TRUE)
#        if (any(col_sds == 0, na.rm = TRUE) && interactive()) {
#            message(sprintf("compute_W_task_from_activations: Found %d parcel(s) with zero variance in activation profiles. Similarities involving these will be affected (likely 0 or NA).", sum(col_sds==0, na.rm=TRUE)))
#        }
#    }
#    sim_matrix_dense <- stats::cor(activation_matrix, method = similarity_method, use = "pairwise.complete.obs")
#  } else if (is.function(similarity_method)) {
#    sim_matrix_dense <- tryCatch({
#      similarity_method(activation_matrix)
#    }, error = function(e) {
#      stop(paste("The provided similarity_method function failed:", e$message))
#    })
#    if (!is.matrix(sim_matrix_dense) || nrow(sim_matrix_dense) != V_p || ncol(sim_matrix_dense) != V_p) {
#      stop("The custom similarity_method function must return a V_p x V_p matrix.")
#    }
#  } else {
#    stop("'similarity_method' must be 'pearson', 'spearman', or a function.")
#  }
#
#  sim_matrix_dense[is.na(sim_matrix_dense)] <- 0
#  diag(sim_matrix_dense) <- 0
#
#  # 2. Sparsify
#  row_indices_list <- vector("list", V_p)
#  col_indices_list <- vector("list", V_p)
#  values_list <- vector("list", V_p)
#
#  if (k_conn_task_pos > 0 || k_conn_task_neg > 0) {
#    for (i in 1:V_p) {
#      node_sims <- sim_matrix_dense[i, ]
#      current_selected_indices <- integer(0)
#      current_selected_values <- numeric(0)
#
#      if (k_conn_task_pos > 0) {
#        pos_candidates_idx <- which(node_sims > 1e-9) 
#        if (length(pos_candidates_idx) > 0) {
#          pos_candidates_vals <- node_sims[pos_candidates_idx]
#          num_to_keep_pos <- min(k_conn_task_pos, length(pos_candidates_vals))
#          ordered_in_pos_group <- order(pos_candidates_vals, decreasing = TRUE)
#          top_pos_in_group_idx <- head(ordered_in_pos_group, num_to_keep_pos)
#
#          current_selected_indices <- c(current_selected_indices, pos_candidates_idx[top_pos_in_group_idx])
#          current_selected_values <- c(current_selected_values, pos_candidates_vals[top_pos_in_group_idx])
#        }
#      }
#
#      if (k_conn_task_neg > 0) {
#        neg_candidates_idx <- which(node_sims < -1e-9) 
#        if (length(neg_candidates_idx) > 0) {
#          neg_candidates_vals <- node_sims[neg_candidates_idx]
#          num_to_keep_neg <- min(k_conn_task_neg, length(neg_candidates_vals))
#          ordered_in_neg_group <- order(neg_candidates_vals, decreasing = FALSE)
#          top_neg_in_group_idx <- head(ordered_in_neg_group, num_to_keep_neg)
#
#          current_selected_indices <- c(current_selected_indices, neg_candidates_idx[top_neg_in_group_idx])
#          current_selected_values <- c(current_selected_values, neg_candidates_vals[top_neg_in_group_idx])
#        }
#      }
#
#      if(length(current_selected_indices) > 0) {
#          row_indices_list[[i]] <- rep(i, length(current_selected_indices))
#          col_indices_list[[i]] <- current_selected_indices
#          values_list[[i]] <- current_selected_values
#      }
#    }
#  }
#
#  final_row_indices <- unlist(row_indices_list[!sapply(row_indices_list, is.null)])
#  final_col_indices <- unlist(col_indices_list[!sapply(col_indices_list, is.null)])
#  final_values <- unlist(values_list[!sapply(values_list, is.null)])
#
#  W_dir_task <- if (length(final_row_indices) > 0) {
#    Matrix::sparseMatrix(
#      i = final_row_indices, j = final_col_indices, x = final_values,
#      dims = c(V_p, V_p), dimnames = list(parcel_names, parcel_names)
#    )
#  } else {
#    Matrix::Matrix(0, nrow=V_p, ncol=V_p, sparse=TRUE, dimnames = list(parcel_names, parcel_names))
#  }
#  W_dir_task <- Matrix::drop0(W_dir_task)
#
#  # 3. Symmetrize
#  W_dir_task_t <- Matrix::t(W_dir_task)
#  W_sum_task <- W_dir_task + W_dir_task_t
#  
#  W_den_task_val <- as.numeric((W_dir_task != 0) + (W_dir_task_t != 0))
#  W_den_task <- Matrix::Matrix(pmax(1, W_den_task_val), nrow=V_p, ncol=V_p, sparse=TRUE)
#
#  W_symmetric_raw_task <- W_sum_task / W_den_task
#  if (inherits(W_symmetric_raw_task, "sparseMatrix")) {
#    if (any(is.nan(W_symmetric_raw_task@x))) W_symmetric_raw_task@x[is.nan(W_symmetric_raw_task@x)] <- 0
#    if (any(is.infinite(W_symmetric_raw_task@x))) W_symmetric_raw_task@x[is.infinite(W_symmetric_raw_task@x)] <- 0
#  } else {
#    W_symmetric_raw_task[is.nan(W_symmetric_raw_task)] <- 0
#    W_symmetric_raw_task[is.infinite(W_symmetric_raw_task)] <- 0
#  }
#  W_symmetric_raw_task <- Matrix::drop0(W_symmetric_raw_task)
#  W_symmetric_task <- Matrix::forceSymmetric(W_symmetric_raw_task, uplo = "U")
#
#  # 4. Z-score non-zero edge weights
#  if (length(W_symmetric_task@x) > 0) { 
#    non_zero_vals <- W_symmetric_task@x
#    mean_val <- mean(non_zero_vals)
#    sd_val <- stats::sd(non_zero_vals)
#    if (is.na(sd_val) || sd_val == 0) { 
#      W_symmetric_task@x <- rep(0, length(non_zero_vals))
#    } else {
#      W_symmetric_task@x <- (non_zero_vals - mean_val) / sd_val
#    }
#    W_task_i <- Matrix::drop0(W_symmetric_task)
#  } else {
#    W_task_i <- W_symmetric_task 
#  }
#  
#  return(as(W_task_i, "dgCMatrix"))
#}
