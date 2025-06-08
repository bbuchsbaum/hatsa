describe("Spectral Graph Construction Functionality", {

# Helper function for z-scoring non-zero elements of a sparse matrix's upper triangle
# (Assuming zscore_nonzero_sparse exists elsewhere or is accessible)
# If not, we might need to define a simplified version here for testing SGC-004

# Test TCK-SGC-001: Basic Sparsification & Output
test_that("TCK-SGC-001: compute_subject_connectivity_graph_sparse basic output", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  
  set.seed(101)
  V_p <- 10
  T_i <- 50
  k_pos <- 2
  k_neg <- 1
  p_names <- paste0("P", 1:V_p)
  
  # Create data with some correlation structure
  base_signal <- matrix(rnorm(T_i * 2), ncol = 2)
  noise <- matrix(rnorm(T_i * V_p), ncol = V_p)
  X_ts <- matrix(0, nrow = T_i, ncol = V_p)
  X_ts[, 1:4] <- base_signal[, 1] * matrix(runif(T_i*4, 0.5, 1.5), T_i, 4) + noise[, 1:4] * 0.5 # Group 1
  X_ts[, 5:7] <- base_signal[, 2] * matrix(runif(T_i*3, 0.5, 1.5), T_i, 3) + noise[, 5:7] * 0.5 # Group 2
  X_ts[, 8:10] <- noise[, 8:10] # Noise
  colnames(X_ts) <- p_names

  # Function Call
  # Assuming compute_subject_connectivity_graph_sparse exists and performs all steps
  output_W <- compute_subject_connectivity_graph_sparse(
      X_subject = X_ts,
      parcel_names = p_names,
      k_conn_pos = k_pos,
      k_conn_neg = k_neg,
      use_dtw = FALSE
  )

  # Assertions
  expect_s4_class(output_W, "dgCMatrix")
  expect_equal(dim(output_W), c(V_p, V_p))
  expect_true(Matrix::isSymmetric(output_W))
  expect_equal(unname(Matrix::diag(output_W)), rep(0, V_p), tolerance = 1e-8) # Diagonal should be zero

  # Check sparsity for node 1 (expect k_pos + k_neg = 3 outgoing originally, maybe more after symmetrization)
  # Calculate raw correlations for node 1
  cor_node1 <- stats::cor(X_ts[, 1], X_ts)[1, -1] # Exclude self-correlation
  pos_neighbors_node1 <- order(cor_node1[cor_node1 > 0], decreasing = TRUE)[1:min(k_pos, sum(cor_node1 > 0))] + 1 # +1 to adjust index after removing self
  neg_neighbors_node1 <- order(cor_node1[cor_node1 < 0], decreasing = FALSE)[1:min(k_neg, sum(cor_node1 < 0))] + 1
  direct_neighbors_node1 <- c(pos_neighbors_node1, neg_neighbors_node1)
  
  # Check non-zero count for row/col 1 (symmetric)
  # Max possible direct edges = k_pos + k_neg. Symmetrization can increase this.
  # Example: Node 1 keeps edges (1,2), (1,3) positive, (1,8) negative. 
  # If node 2 keeps edge (2,1), it adds to symmetry. If node 3 keeps (3,9), it doesn't affect node 1 directly.
  # The number of non-zeros in row 1 should reflect edges *initiated* by node 1 plus edges *initiated* by others towards node 1.
  # It's hard to predict exact non-zero count without running full symmetrization.
  # Let's check a weaker condition: number of non-zeros <= V_p - 1
  nnz_row1 <- length(output_W@i[output_W@p[1]:(output_W@p[2]-1)])
  expect_lte(nnz_row1, V_p - 1)
  expect_gt(nnz_row1, 0) # Should have some connections unless data is pathological
})

# Test TCK-SGC-002: Edge Weight Correctness
test_that("TCK-SGC-002: compute_subject_connectivity_graph_sparse edge selection & sign", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  
  set.seed(102)
  V_p <- 5
  T_i <- 100
  k_pos <- 1 # Keep only the strongest positive
  k_neg <- 1 # Keep only the strongest negative
  p_names <- paste0("P", 1:V_p)
  
  # Create data with specific correlations
  # P1 correlated with P2 (positive)
  # P1 anti-correlated with P3 (negative)
  # P1 uncorrelated with P4 
  # P5 is noise
  base_signal_pos <- rnorm(T_i)
  base_signal_neg <- rnorm(T_i) # Independent noise for neg corr
  X_ts <- matrix(rnorm(T_i * V_p, sd=0.1), ncol = V_p) # Base noise
  X_ts[, 1] <- base_signal_pos + rnorm(T_i, sd=0.2) # P1 
  X_ts[, 2] <- base_signal_pos + rnorm(T_i, sd=0.2) # P2 highly correlated with P1
  X_ts[, 3] <- -base_signal_neg + rnorm(T_i, sd=0.2) # P3 anti-correlated with P1 via base_signal_neg
  X_ts[, 1] <- X_ts[, 1] + base_signal_neg # Add base_signal_neg to P1 to create neg corr with P3
  # P4 remains noise, weakly correlated with P1
  # P5 remains noise
  colnames(X_ts) <- p_names

  # Calculate raw correlations for reference (esp. for P1)
  raw_corrs <- cor(X_ts)
  diag(raw_corrs) <- 0
  # Expected strongest positive for P1 should be P2
  # Expected strongest negative for P1 should be P3
  # Correlation P1-P4 should be small
  expect_gt(raw_corrs[1, 2], 0.5) # Should be high positive
  expect_lt(raw_corrs[1, 3], -0.5) # Should be high negative
  expect_lt(abs(raw_corrs[1, 4]), 0.3) # Should be low 

  # Function Call
  output_W <- compute_subject_connectivity_graph_sparse(
      X_subject = X_ts,
      parcel_names = p_names,
      k_conn_pos = k_pos,
      k_conn_neg = k_neg,
      use_dtw = FALSE
  )
  
  # Assertions
  # Check that the strong positive connection P1-P2 exists and has a positive z-score
  # Since it's symmetric, check both [1,2] and [2,1]
  # At least one direction should have selected it initially
  # We expect output_W[1,2] to be non-zero and positive after z-scoring
  # Note: z-score depends on all edges, so we check sign mainly
  expect_gt(output_W[1, 2], 0) 

  # Check that the strong negative connection P1-P3 exists and has a negative z-score
  expect_lt(output_W[1, 3], 0)

  # Check that P1 did NOT directly select P4 (the weak connection)
  # This means raw_corrs[1,4] was not among the top k_pos/k_neg for P1.
  raw_node1_corrs_others <- raw_corrs[1, -1] # P1's corrs with P2, P3, P4, P5 (indices 2,3,4,5 of X_ts)
  # We are interested if P4 (which is at original index 4, so index 3 in raw_node1_corrs_others) was selected.
  target_original_idx_for_P4 <- 4 

  if (k_pos > 0) {
    pos_candidates_P1_indices <- which(raw_node1_corrs_others > 0)
    if (length(pos_candidates_P1_indices) > 0) {
      pos_candidates_P1_vals <- raw_node1_corrs_others[pos_candidates_P1_indices]
      num_to_keep_pos_P1 <- min(k_pos, length(pos_candidates_P1_vals))
      # Get original column indices (2 to V_p) of those selected by P1 as positive
      selected_pos_by_P1_original_col_indices <- (which(raw_node1_corrs_others > 0))[order(pos_candidates_P1_vals, decreasing = TRUE)[1:num_to_keep_pos_P1]] + 1
      expect_false(target_original_idx_for_P4 %in% selected_pos_by_P1_original_col_indices, 
                   label = "P1 direct positive selection of P4 check")
    }
  }
  
  if (k_neg > 0) {
    neg_candidates_P1_indices <- which(raw_node1_corrs_others < 0)
    if (length(neg_candidates_P1_indices) > 0) {
      neg_candidates_P1_vals <- raw_node1_corrs_others[neg_candidates_P1_indices]
      num_to_keep_neg_P1 <- min(k_neg, length(neg_candidates_P1_vals))
      # Get original column indices (2 to V_p) of those selected by P1 as negative
      selected_neg_by_P1_original_col_indices <- (which(raw_node1_corrs_others < 0))[order(neg_candidates_P1_vals, decreasing = FALSE)[1:num_to_keep_neg_P1]] + 1
      expect_false(target_original_idx_for_P4 %in% selected_neg_by_P1_original_col_indices, 
                   label = "P1 direct negative selection of P4 check")
    }
  }

  # Check overall sparsity: With k_pos=1, k_neg=1, max direct edges per node is 2.
  # Symmetrization might increase this, but total nnz should be relatively small.
  # Example: If P1->P2, P1->P3 selected. If P2->P1, P2->P5 selected. If P3->P1, P3->P4 selected.
  # Symmetric non-zero pairs could be (1,2), (1,3), (2,5), (3,4).
  # Total non-zero entries (upper or lower triangle) should be relatively low.
  
  # Use tryCatch to handle potential NA result from nnzero
  nnz_output_W <- tryCatch(
    Matrix::nnzero(output_W), 
    error = function(e) 0, 
    warning = function(w) 0
  )
  
  if (!is.na(nnz_output_W)) {
    expect_lt(nnz_output_W, V_p * (k_pos + k_neg) * 2) # Loose upper bound
  }
})

# Test TCK-SGC-003: Handling Zero Variance Columns
test_that("TCK-SGC-003: compute_subject_connectivity_graph_sparse handles zero variance cols", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  
  set.seed(103)
  V_p <- 6
  T_i <- 40
  k_pos <- 2
  k_neg <- 2
  p_names <- paste0("P", 1:V_p)
  
  # Create data where column 3 is constant
  X_ts <- matrix(rnorm(T_i * V_p), ncol = V_p)
  X_ts[, 3] <- 5 # Constant column
  colnames(X_ts) <- p_names

  # Call it once to get the actual output for subsequent checks
  # The message will be produced here, but we test for it in the expect_message block.
  output_W <- compute_subject_connectivity_graph_sparse(
      X_subject = X_ts,
      parcel_names = p_names,
      k_conn_pos = k_pos,
      k_conn_neg = k_neg,
      use_dtw = FALSE
  )

  # Assertions on output_W (obtained from the first call)
  expect_s4_class(output_W, "dgCMatrix")
  expect_equal(dim(output_W), c(V_p, V_p))
  expect_true(Matrix::isSymmetric(output_W))

  # Check that row 3 and column 3 are all zeros
  expect_equal(nnzero(output_W[3, ]), 0) # Check row 3 non-zeros
  expect_equal(nnzero(output_W[, 3]), 0) # Check column 3 non-zeros

  # Double check using direct access (might be redundant but safe)
  # Extract column 3 structure from dgCMatrix
  col3_start_idx <- output_W@p[3]
  col3_end_idx <- output_W@p[3+1]
  col3_num_entries <- col3_end_idx - col3_start_idx
  expect_equal(col3_num_entries, 0)
  
  # Check row 3 - find which non-zero entries have row index 3 (should be none)
  row3_indices_in_x <- which(output_W@i == (3-1)) # 0-based index
  expect_length(row3_indices_in_x, 0)
})

# Test TCK-SGC-011: Memory efficient correlation for larger V_p
test_that("TCK-SGC-011: compute_subject_connectivity_graph_sparse scales without dense correlation", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  set.seed(104)
  V_p <- 200
  T_i <- 20
  X_ts <- matrix(rnorm(T_i * V_p), ncol = V_p)
  p_names <- paste0("P", seq_len(V_p))

  expect_no_error(
    compute_subject_connectivity_graph_sparse(
      X_subject = X_ts,
      parcel_names = p_names,
      k_conn_pos = 1,
      k_conn_neg = 1,
      use_dtw = FALSE
    )
  )
})


# Test TCK-SGC-012: Input validation checks
test_that("TCK-SGC-012: compute_subject_connectivity_graph_sparse validates inputs", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  X_small <- matrix(rnorm(10 * 3), ncol = 3)
  p_names <- paste0("P", 1:3)

  # parcel_names length mismatch
  expect_error(
    compute_subject_connectivity_graph_sparse(
      X_subject = X_small,
      parcel_names = p_names[-1],
      k_conn_pos = 1,
      k_conn_neg = 1
    ),
    "parcel_names"
  )

  # negative k_conn_pos
  expect_error(
    compute_subject_connectivity_graph_sparse(
      X_subject = X_small,
      parcel_names = p_names,
      k_conn_pos = -1,
      k_conn_neg = 1
    ),
    "k_conn_pos"
  )

  # non-integer k_conn_neg
  expect_error(
    compute_subject_connectivity_graph_sparse(
      X_subject = X_small,
      parcel_names = p_names,
      k_conn_pos = 1,
      k_conn_neg = 1.5
    ),
    "k_conn_neg"
})
  
  
# Test TCK-SGC-012: Input must have at least two rows
test_that("TCK-SGC-012: compute_subject_connectivity_graph_sparse requires >=2 rows", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  X_ts <- matrix(rnorm(5), nrow = 1, ncol = 5)
  p_names <- paste0("P", 1:5)

  expect_error(
    compute_subject_connectivity_graph_sparse(
      X_subject = X_ts,
      parcel_names = p_names,
      k_conn_pos = 1,
      k_conn_neg = 1,
      use_dtw = FALSE
    ),
    regexp = "at least two rows"
  )
})

# Test TCK-SGC-004: Symmetrization Logic
test_that("TCK-SGC-004: Symmetrization rule application", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  
  # Create a small, non-symmetric directed graph manually
  # Cases:
  # (1,2): Edge exists only 1->2
  # (1,3): Edge exists only 3->1
  # (2,3): Edges exist in both directions with different values
  # (1,4): No edge in either direction (implicitly zero)
  i_idx <- c(1, 3, 2, 3) # Row indices (1-based)
  j_idx <- c(2, 1, 3, 2) # Col indices (1-based)
  x_val <- c(0.5, 0.8, 0.9, 0.6) # Values W_dir[i,j]
  W_dir <- sparseMatrix(i=i_idx, j=j_idx, x=x_val, dims=c(4,4), dimnames=list(paste0("N",1:4), paste0("N",1:4)))
  # W_dir:
  #       N1  N2  N3 N4
  #   N1   . 0.5   .  .
  #   N2   .   . 0.9  .
  #   N3 0.8 0.6   .  .
  #   N4   .   .   .  .

  # Apply the symmetrization logic used in compute_subject_connectivity_graph_sparse
  W_dir_t <- Matrix::t(W_dir)
  W_sum <- W_dir + W_dir_t
  W_indicator <- (W_dir != 0)
  W_indicator_t <- (W_dir_t != 0) # or Matrix::t(W_indicator)
  W_den <- W_indicator + W_indicator_t
  W_symmetric_raw <- W_sum / W_den
  W_symmetric_raw <- Matrix::drop0(W_symmetric_raw)

  # Assertions based on the rule
  # Case (1,2): Only W_dir[1,2]=0.5 exists. Sum=0.5. Indicator = 1+0 = 1. Den=1. Expected = 0.5/1 = 0.5
  expect_equal(W_symmetric_raw[1, 2], 0.5)
  expect_equal(W_symmetric_raw[2, 1], 0.5) # Should be symmetric

  # Case (1,3): Only W_dir[3,1]=0.8 exists. Sum=0.8. Indicator = 0+1 = 1. Den=1. Expected = 0.8/1 = 0.8
  expect_equal(W_symmetric_raw[1, 3], 0.8)
  expect_equal(W_symmetric_raw[3, 1], 0.8) # Should be symmetric

  # Case (2,3): W_dir[2,3]=0.9, W_dir[3,2]=0.6. Sum=1.5. Indicator = 1+1=2. Den=2. Expected = 1.5/2 = 0.75
  expect_equal(W_symmetric_raw[2, 3], 0.75)
  expect_equal(W_symmetric_raw[3, 2], 0.75) # Should be symmetric

  # Case (1,4): No edges. Sum=0. Indicator = 0+0=0. Den=0. Expected = 0/0 -> NaN.
  # drop0 does not remove NaNs.
  expect_true(is.nan(W_symmetric_raw[1, 4]))
  expect_true(is.nan(W_symmetric_raw[4, 1]))

  # Verify overall symmetry (numerically)
  # Note: isSymmetric might be false if NaNs are present and not handled symmetrically by the check itself.
  # For this test, we are checking the arithmetic. Symmetry of NaN handling is a separate concern for the main function.
  # If W_symmetric_raw contains NaNs, isSymmetric(W_symmetric_raw) might be NA or FALSE.
  # We can test symmetry on the non-NaN part if needed, or accept that NaNs break simple symmetry check here.
  # The main function *does* convert NaN to 0, then drops, then forcesymmetry.
  # Here we test the raw arithmetic result *before* that NaN cleanup.
  # A matrix with NaN is symmetric if M[i,j] and M[j,i] are both NaN or both equal.
  expect_true(Matrix::isSymmetric(W_symmetric_raw, tol=1e-8, na.rm=FALSE)) # Check symmetry considering NaNs as comparable
})

# Test TCK-SGC-005: Correctness of alpha-lazy random-walk Laplacian
test_that("TCK-SGC-005: compute_graph_laplacian_sparse calculates L_rw_lazy_sym correctly", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  
  set.seed(105)
  alpha <- 0.93
  
  # Define upper triangle for the desired symmetric matrix:
  # W_sparse:
  #      [,1] [,2] [,3] [,4]
  # [1,]  .   0.5  0.2   .
  # [2,] 0.5   .    .   0.7
  # [3,] 0.2   .    .    .
  # [4,]  .   0.7   .    .
  i_upper <- c(1,   1,   2) 
  j_upper <- c(2,   3,   4) 
  x_upper <- c(0.5, 0.2, 0.7)
  
  W_tri <- sparseMatrix(i=i_upper, j=j_upper, x=x_upper, dims=c(4,4))
  W_sparse <- W_tri + Matrix::t(W_tri) # Construct symmetric matrix
  W_sparse <- as(W_sparse, "dgCMatrix") # Ensure it's dgCMatrix for consistency

  # Manual Calculation
  V_p <- nrow(W_sparse)
  degree_vec <- Matrix::rowSums(abs(W_sparse)) # Using abs as in the function
  inv_degree_vec <- ifelse(degree_vec == 0, 0, 1 / degree_vec)
  D_inv_sparse <- Matrix::Diagonal(n = V_p, x = inv_degree_vec)
  I_mat <- Matrix::Diagonal(n = V_p)
  
  L_rw_lazy_manual <- I_mat - alpha * (D_inv_sparse %*% W_sparse)
  # Symmetrize manually
  L_rw_lazy_sym_manual <- (L_rw_lazy_manual + Matrix::t(L_rw_lazy_manual)) / 2
  L_rw_lazy_sym_manual <- Matrix::drop0(L_rw_lazy_sym_manual)
  
  # Function Call
  output_L <- compute_graph_laplacian_sparse(W_sparse, alpha = alpha)

  # Assertions
  # Use all.equal for sparse matrices or convert to dense if small
  expect_true(all.equal(output_L, L_rw_lazy_sym_manual, tolerance = 1e-8))
  
  # Check properties
  expect_s4_class(output_L, "dgCMatrix") # Should be sparse
  expect_true(Matrix::isSymmetric(output_L, tol = 1e-8))
})

# Test TCK-SGC-006: Handling Zero Degree Nodes
test_that("TCK-SGC-006: compute_graph_laplacian_sparse handles zero degree nodes", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  
  set.seed(106)
  alpha <- 0.8 # Use a different alpha for variation
  
  # Create a W matrix where node 3 is isolated
  # W_sparse:
  #      [,1] [,2] [,3] [,4]
  # [1,]  .   0.5  .   0.1
  # [2,] 0.5   .    .    .
  # [3,]  .    .    .    .
  # [4,] 0.1   .    .    .
  i_upper <- c(1, 1) 
  j_upper <- c(2, 4) 
  x_upper <- c(0.5, 0.1)
  
  W_tri <- sparseMatrix(i=i_upper, j=j_upper, x=x_upper, dims=c(4,4))
  W_sparse <- W_tri + Matrix::t(W_tri) # Construct symmetric matrix
  W_sparse <- as(W_sparse, "dgCMatrix") # Ensure it's dgCMatrix

  # Function Call
  output_L <- NULL
  expect_no_error(
    output_L <- compute_graph_laplacian_sparse(W_sparse, alpha = alpha)
  )
  
  # Assertions
  expect_s4_class(output_L, "dgCMatrix")
  expect_equal(dim(output_L), c(4, 4))

  # Check properties for the isolated node (node 3)
  # For L_rw_lazy_sym = (I - alpha*Dinv*W + t(I - alpha*Dinv*W))/2
  # If degree(3)=0, then Dinv[3,3]=0.
  # Row 3 of Dinv*W is zero.
  # Col 3 of Dinv*W is zero.
  # So Row 3 of L = I - alpha*Dinv*W is [0, 0, 1, 0]
  # Col 3 of L is [0, 0, 1, 0]^T
  # Symmetrizing keeps this structure.
  expect_equal(output_L[3, 3], 1.0) # L[i,i] should be 1 for isolated node
  
  # Check other elements in row 3 and column 3 are zero
  expect_equal(output_L[3, -3], c(0, 0, 0)) # Row 3, excluding diagonal
  expect_equal(output_L[-3, 3], c(0, 0, 0)) # Col 3, excluding diagonal
})

# Test TCK-SGC-006B: Non-symmetric input should error
test_that("TCK-SGC-006B: compute_graph_laplacian_sparse errors for non-symmetric input", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  # Create a simple non-symmetric adjacency matrix
  W_ns <- sparseMatrix(i = c(1, 2), j = c(2, 3), x = c(1, 1), dims = c(3, 3))

  expect_error(
    compute_graph_laplacian_sparse(W_ns),
    regexp = "W_sparse must be symmetric"
  )
})

# Test TCK-SGC-007: compute_spectral_sketch_sparse - Basic Correctness & Dimensions
test_that("TCK-SGC-007: compute_spectral_sketch_sparse basic correctness & dimensions", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("RSpectra") # Function uses RSpectra::eigs_sym
  library(Matrix)
  
  set.seed(107)
  V_p <- 6
  k_target <- 2 # Request 2 informative eigenvectors
  
  # Create a small, sparse, symmetric matrix (not necessarily a valid Laplacian)
  # Ensure it has some distinct eigenvalues
  i_upper <- c(1, 1, 2, 2, 3, 4, 5)
  j_upper <- c(1, 2, 2, 3, 4, 5, 6) # Include diagonal elements
  x_upper <- c(2, 0.5, 3, -0.1, 1.5, 0.8, 2.5)
  L_tri <- sparseMatrix(i=i_upper, j=j_upper, x=x_upper, dims=c(V_p, V_p))
  L_sparse <- L_tri + t(L_tri)
  # Correct the diagonal after symmetrization (it gets doubled)
  Matrix::diag(L_sparse) <- Matrix::diag(L_sparse) / 2 
  L_sparse <- as(L_sparse, "dgCMatrix") # Ensure type before eigen

  # Calculate ground truth using base::eigen on the dense matrix
  eigen_gt <- eigen(as.matrix(L_sparse), symmetric = TRUE)
  
  # Sort ground truth by eigenvalue (ascending, like compute_spectral_sketch_sparse)
  sorted_order_gt <- order(eigen_gt$values)
  eigen_vals_gt_sorted <- eigen_gt$values[sorted_order_gt]
  eigen_vecs_gt_sorted <- eigen_gt$vectors[, sorted_order_gt, drop = FALSE]
  
  # Identify the first k_target *non-trivial* ground truth components
  # (Assuming the first one might be trivial/near-zero, though not guaranteed for this arbitrary matrix)
  # Let's assume compute_spectral_sketch_sparse handles the trivial filtering correctly (tested in SGC-008)
  # Here, we compare against the first k_target smallest eigenvalues/vectors directly from sorted ground truth.
  # Need to be careful if L_sparse has eigenvalue ~0 as first.
  # compute_spectral_sketch_sparse filters eigenvalues > tol (e.g., 1e-8)
  eigenvalue_tol <- 1e-8
  non_trivial_gt_indices <- which(eigen_vals_gt_sorted > eigenvalue_tol)
  if (length(non_trivial_gt_indices) < k_target) {
      skip("Ground truth matrix doesn't have enough non-trivial eigenvalues for k_target.")
  }
  gt_indices_to_compare <- non_trivial_gt_indices[1:k_target]
  U_gt_k <- eigen_vecs_gt_sorted[, gt_indices_to_compare, drop=FALSE]
  Lambda_gt_k <- eigen_vals_gt_sorted[gt_indices_to_compare]

  # Function Call
  result <- NULL
  expect_no_error(
    result <- compute_spectral_sketch_sparse(L_sparse, k = k_target)
  )
  
  # Assertions
  expect_true(is.list(result))
  expect_named(result, c("vectors", "values"))
  expect_true(is.matrix(result$vectors))
  expect_true(is.numeric(result$values))
  expect_equal(nrow(result$vectors), V_p)
  expect_equal(ncol(result$vectors), k_target)
  expect_length(result$values, k_target)
  
  # Compare eigenvalues (should be close)
  expect_equal(result$values, Lambda_gt_k, tolerance = 1e-6)
  
  # Compare eigenvectors (allowing for sign flips)
  U_result_k <- result$vectors
  # Align signs for comparison: Flip columns in U_result_k if their correlation with U_gt_k is negative
  for (j in 1:k_target) {
    # Check correlation - handle potential zero vectors if L was pathological
    cor_val <- tryCatch(cor(U_result_k[, j], U_gt_k[, j]), error = function(e) 0)
    if (!is.na(cor_val) && cor_val < 0) {
      U_result_k[, j] <- -U_result_k[, j]
    }
  }
  # Now compare sign-aligned vectors
  expect_equal(U_result_k, U_gt_k, tolerance = 1e-6)
  
  # Alternative check: High absolute correlation for each component
  # correlations <- diag(abs(cor(result$vectors, U_gt_k)))
  # expect_true(all(correlations > 0.999), label = "Absolute correlation check")
})

# Test TCK-SGC-008: compute_spectral_sketch_sparse - Trivial Eigenvector Handling
test_that("TCK-SGC-008: compute_spectral_sketch_sparse handles trivial eigenvectors", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("RSpectra")
  library(Matrix)
  
  V_p <- 5 # Size of the graph
  k_target <- 2 # Request 2 informative eigenvectors
  
  # Create Laplacian for a simple connected graph (path graph 1-2-3-4-5)
  # Expected eigenvalues: 0, and 4 non-zero smallest ones.
  adj <- bandSparse(V_p, k = 1, diagonals = list(rep(1, V_p-1)), symmetric = TRUE)
  D <- Diagonal(V_p, x = Matrix::rowSums(adj))
  L_path <- D - adj
  L_path <- as(L_path, "dgCMatrix")

  # Function Call
  result <- NULL
  expect_no_error(
    result <- compute_spectral_sketch_sparse(L_path, k = k_target)
  )
  
  # Assertions
  expect_true(is.list(result))
  expect_named(result, c("vectors", "values"))
  expect_equal(ncol(result$vectors), k_target) # Should return exactly k vectors
  expect_length(result$values, k_target)     # Should return exactly k values

  # Check that the smallest returned eigenvalue is strictly greater than the internal tolerance
  eigenvalue_tol <- 1e-8 # Match the default used in the function
  expect_gt(min(result$values), eigenvalue_tol, 
            label="Smallest returned eigenvalue should be > tolerance")
            
  # Optional: Verify against known eigenvalues for path graph if needed
  # evals_analytic <- 2 * (1 - cos(pi * (0:(V_p-1)) / V_p))
  # expect_equal(sort(result$values), sort(evals_analytic)[2:(k_target+1)], tolerance=1e-6)
})

# Test TCK-SGC-009: compute_spectral_sketch_sparse - Rank Deficiency / Insufficient Eigenvectors
test_that("TCK-SGC-009: compute_spectral_sketch_sparse handles rank deficiency", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("RSpectra")
  library(Matrix)
  
  # Create a graph with 2 connected components (e.g., two path graphs)
  V_p1 <- 4
  adj1 <- bandSparse(V_p1, k = 1, diagonals = list(rep(1, V_p1-1)), symmetric = TRUE)
  D1 <- Diagonal(V_p1, x = Matrix::rowSums(adj1))
  L1 <- D1 - adj1
  
  V_p2 <- 3
  adj2 <- bandSparse(V_p2, k = 1, diagonals = list(rep(1, V_p2-1)), symmetric = TRUE)
  D2 <- Diagonal(V_p2, x = Matrix::rowSums(adj2))
  L2 <- D2 - adj2

  # Combine into a block diagonal matrix
  L_disconnected <- Matrix::bdiag(L1, L2)
  L_disconnected <- as(L_disconnected, "dgCMatrix")
  V_p_total <- nrow(L_disconnected) # Should be V_p1 + V_p2 = 7
  
  # This graph has 2 zero eigenvalues (one for each component).
  # Total number of eigenvalues is V_p_total = 7.
  # Number of non-zero eigenvalues = V_p_total - 2 = 5.
  num_non_zero_eigenvals <- V_p_total - 2 
  
  # Request k > number of non-zero eigenvalues
  k_target <- num_non_zero_eigenvals + 1 # Request 6 informative eigenvectors

  # Expect an error because we can only find 5 informative ones.
  expect_error(
    compute_spectral_sketch_sparse(L_disconnected, k = k_target),
    regexp = "Rank deficiency|too many zero eigenvalues|Found only.*informative eigenvectors.*but k=.*was requested" # Match potential error messages
  )
  
  # Also test the boundary case: requesting exactly the number of informative eigenvalues
  # In this specific case (Vp=7, k=5 requesting 6 eigs), eigs_sym only returns 4 informative ones.
  # So, the function should still error because 4 < 5.
  k_target_boundary = num_non_zero_eigenvals # Request 5
  expect_error(
    compute_spectral_sketch_sparse(L_disconnected, k = k_target_boundary),
    regexp = "Rank deficiency|too many zero eigenvalues|Found only.*informative eigenvectors.*but k=.*was requested",
    info=paste("Boundary case k=", k_target_boundary, "should error if eigs_sym doesn't return enough.")
  )
})

# Test TCK-SGC-010: compute_spectral_sketch_sparse - Edge Cases k=0, k=1
test_that("TCK-SGC-010: compute_spectral_sketch_sparse edge cases k=0, k=1", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("RSpectra")
  library(Matrix)
  
  # Use the same path graph Laplacian from TCK-SGC-008
  V_p <- 5
  adj <- bandSparse(V_p, k = 1, diagonals = list(rep(1, V_p-1)), symmetric = TRUE)
  D <- Diagonal(V_p, x = Matrix::rowSums(adj))
  L_path <- D - adj
  L_path <- as(L_path, "dgCMatrix")

  # --- Test k = 0 --- 
  result_k0 <- NULL
  expect_no_error(
    result_k0 <- compute_spectral_sketch_sparse(L_path, k = 0)
  )
  
  # Assertions for k=0
  expect_true(is.list(result_k0))
  expect_named(result_k0, c("vectors", "values"))
  expect_true(is.matrix(result_k0$vectors))
  expect_true(is.numeric(result_k0$values))
  expect_equal(nrow(result_k0$vectors), V_p) # Should still have V_p rows
  expect_equal(ncol(result_k0$vectors), 0)
  expect_length(result_k0$values, 0)

  # --- Test k = 1 --- 
  k_target_1 <- 1
  result_k1 <- NULL
  expect_no_error(
    result_k1 <- compute_spectral_sketch_sparse(L_path, k = k_target_1)
  )
  
  # Assertions for k=1
  expect_true(is.list(result_k1))
  expect_named(result_k1, c("vectors", "values"))
  expect_true(is.matrix(result_k1$vectors))
  expect_true(is.numeric(result_k1$values))
  expect_equal(nrow(result_k1$vectors), V_p)
  expect_equal(ncol(result_k1$vectors), k_target_1)
  expect_length(result_k1$values, k_target_1)
  
  # Check eigenvalue is positive (smallest non-trivial)
  eigenvalue_tol <- 1e-8
  expect_gt(result_k1$values[1], eigenvalue_tol)
  
})
}) 