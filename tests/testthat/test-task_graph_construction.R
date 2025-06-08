# tests/testthat/test-task_graph_construction.R

library(testthat)
library(Matrix)
# Assuming 'hatsa' package functions are available during testing, or use devtools::load_all()

describe("compute_W_task_from_activations", {
  V_p <- 10 # Number of parcels
  C <- 5    # Number of conditions/features
  parcel_names <- paste0("P", 1:V_p)

  # Helper to check for z-scoring (approximate)
  is_z_scored <- function(x) {
    if (length(x) < 2) return(TRUE) # Cannot assess sd for < 2 values
    abs(mean(x)) < 1e-9 && abs(stats::sd(x) - 1) < 1e-9
  }

  it("produces a valid sparse graph with pearson correlation", {
    activation_matrix <- matrix(rnorm(C * V_p), nrow = C, ncol = V_p)
    colnames(activation_matrix) <- parcel_names
    
    W_task <- compute_W_task_from_activations(
      activation_matrix = activation_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 2,
      k_conn_task_neg = 1,
      similarity_method = "pearson"
    )
    
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    expect_equal(rownames(W_task), parcel_names)
    expect_equal(colnames(W_task), parcel_names)
    expect_true(Matrix::isSymmetric(W_task))
    if (length(W_task@x) > 1) { # Only check z-score if there are edges
      # This is an approximate check as k-NN and symmetrization affect distribution
      # expect_true(is_z_scored(W_task@x)) 
      # A looser check might be that values are not wildly out of z-score range
      expect_true(all(abs(W_task@x) < 5)) # Arbitrary sanity check for z-scored values
    }
  })

  it("produces a valid sparse graph with spearman correlation", {
    activation_matrix <- matrix(rnorm(C * V_p), nrow = C, ncol = V_p)
    colnames(activation_matrix) <- parcel_names
    
    W_task <- compute_W_task_from_activations(
      activation_matrix = activation_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 2,
      k_conn_task_neg = 1,
      similarity_method = "spearman"
    )
    
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    # ... other common checks ...
    expect_true(Matrix::isSymmetric(W_task))
  })
  
  it("handles custom similarity function", {
    activation_matrix <- matrix(rnorm(C * V_p), nrow = C, ncol = V_p)
    custom_sim_func <- function(act_mat) {
      # Example: simple cosine similarity (not as robust as cor for this structure)
      # For testing, ensure it returns V_p x V_p
      m <- t(act_mat) # V_p x C
      sim <- m %*% t(m) / (sqrt(rowSums(m^2)) %*% t(sqrt(rowSums(m^2))))
      diag(sim) <- 0 # Ensure diagonal is 0 before kNN
      return(sim)
    }
    
    W_task <- compute_W_task_from_activations(
      activation_matrix = activation_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 1,
      k_conn_task_neg = 1,
      similarity_method = custom_sim_func
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    expect_true(Matrix::isSymmetric(W_task))
  })
  
  it("handles zero-variance columns in activation_matrix", {
    activation_matrix <- matrix(rnorm(C * V_p), nrow = C, ncol = V_p)
    activation_matrix[, 1] <- rep(1, C) # First parcel has zero variance
    colnames(activation_matrix) <- parcel_names
    
    # Expect the message from our function, suppress potential underlying cor() warning for this specific test
    expect_message( 
      suppressWarnings(W_task <- compute_W_task_from_activations(
        activation_matrix = activation_matrix,
        parcel_names = parcel_names,
        k_conn_task_pos = 2,
        k_conn_task_neg = 1,
        similarity_method = "pearson"
      )),
      regexp = "Found 1 parcel\\(s\\) with zero variance"
    )
    expect_s4_class(W_task, "dgCMatrix")
    # Further checks could be added here if W_task is not NA/NaN
  })
  
  it("handles k_conn_pos=0 and k_conn_task_neg=0", {
    activation_matrix <- matrix(rnorm(C * V_p), nrow = C, ncol = V_p)
    colnames(activation_matrix) <- parcel_names
    W_task <- compute_W_task_from_activations(
      activation_matrix = activation_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 0,
      k_conn_task_neg = 0,
      similarity_method = "pearson"
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    expect_equal(length(W_task@x), 0) # Should be an empty graph (no non-zero elements)
  })
  
  it("handles input with 0 parcels (V_p=0)", {
    activation_matrix_empty_vp <- matrix(nrow = C, ncol = 0)
    W_task <- compute_W_task_from_activations(
        activation_matrix = activation_matrix_empty_vp,
        parcel_names = character(0),
        k_conn_task_pos = 1, k_conn_task_neg = 1
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(0,0))
  })
  
  it("handles input with 0 conditions (C=0)", {
    activation_matrix_empty_c <- matrix(nrow = 0, ncol = V_p)
    colnames(activation_matrix_empty_c) <- parcel_names
    
    expect_warning(
      W_task_c0 <- compute_W_task_from_activations(
          activation_matrix = activation_matrix_empty_c,
          parcel_names = parcel_names,
          k_conn_task_pos = 1, k_conn_task_neg = 1
      ),
      regexp = "activation_matrix has 0 row\\(s\\)"
    )
    expect_s4_class(W_task_c0, "dgCMatrix")
    expect_equal(dim(W_task_c0), c(V_p,V_p))
    expect_equal(length(W_task_c0@x), 0)

    activation_matrix_c1 <- matrix(rnorm(V_p), nrow = 1, ncol = V_p)
    colnames(activation_matrix_c1) <- parcel_names
    expect_warning(
      W_task_c1 <- compute_W_task_from_activations(
          activation_matrix = activation_matrix_c1,
          parcel_names = parcel_names,
          k_conn_task_pos = 1, k_conn_task_neg = 1
      ),
      regexp = "activation_matrix has 1 row\\(s\\)"
    )
    expect_s4_class(W_task_c1, "dgCMatrix")
    expect_equal(dim(W_task_c1), c(V_p,V_p))
    expect_equal(length(W_task_c1@x), 0)
  })

  it("errors for invalid k_conn_task_pos or k_conn_task_neg", {
    activation_matrix <- matrix(rnorm(C * V_p), nrow = C, ncol = V_p)
    expect_error(
      compute_W_task_from_activations(
        activation_matrix = activation_matrix,
        parcel_names = parcel_names,
        k_conn_task_pos = -1,
        k_conn_task_neg = 1
      ),
      "k_conn_task_pos"
    )
    expect_error(
      compute_W_task_from_activations(
        activation_matrix = activation_matrix,
        parcel_names = parcel_names,
        k_conn_task_pos = 1,
        k_conn_task_neg = -1
      ),
      "k_conn_task_neg"
    )
  })

}) # end describe compute_W_task_from_activations


describe("compute_W_task_from_encoding", {
  V_p <- 10 # Number of parcels
  N_features <- 5 # Number of encoding features
  parcel_names <- paste0("P", 1:V_p)

  is_z_scored <- function(x) { # Helper from above, can be defined globally in test file
    if (length(x) < 2) return(TRUE)
    abs(mean(x)) < 1e-9 && abs(stats::sd(x) - 1) < 1e-9
  }

  it("produces a valid sparse graph with pearson correlation", {
    encoding_matrix <- matrix(rnorm(V_p * N_features), nrow = V_p, ncol = N_features)
    rownames(encoding_matrix) <- parcel_names

    W_task <- compute_W_task_from_encoding(
      encoding_weights_matrix = encoding_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 2,
      k_conn_task_neg = 1,
      similarity_method = "pearson"
    )

    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    expect_equal(rownames(W_task), parcel_names)
    expect_equal(colnames(W_task), parcel_names)
    expect_true(Matrix::isSymmetric(W_task))
    if (length(W_task@x) > 1) {
        expect_true(all(abs(W_task@x) < 5)) # Sanity check for z-scored values
    }
  })
  
  it("produces a valid sparse graph with spearman correlation", {
    encoding_matrix <- matrix(runif(V_p * N_features), nrow = V_p, ncol = N_features) # Use runif for spearman
    rownames(encoding_matrix) <- parcel_names
    W_task <- compute_W_task_from_encoding(
      encoding_weights_matrix = encoding_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 2,
      k_conn_task_neg = 1,
      similarity_method = "spearman"
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    expect_true(Matrix::isSymmetric(W_task))
  })
  
  it("handles custom similarity function for encoding weights", {
    encoding_matrix <- matrix(rnorm(V_p * N_features), nrow = V_p, ncol = N_features)
    custom_sim_func_encoding <- function(enc_mat) { # V_p x N_features
      # Example: simple dot product similarity between parcel encoding profiles
      sim <- enc_mat %*% t(enc_mat) 
      diag(sim) <- 0
      return(sim) # Should return V_p x V_p
    }
    
    W_task <- compute_W_task_from_encoding(
      encoding_weights_matrix = encoding_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 1,
      k_conn_task_neg = 1,
      similarity_method = custom_sim_func_encoding
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    expect_true(Matrix::isSymmetric(W_task))
  })
  
  it("handles zero-variance rows in encoding_matrix", {
    encoding_matrix <- matrix(rnorm(V_p * N_features), nrow = V_p, ncol = N_features)
    encoding_matrix[1, ] <- rep(0.5, N_features) # First parcel has zero variance
    rownames(encoding_matrix) <- parcel_names
    
    expect_message(
      suppressWarnings(W_task <- compute_W_task_from_encoding(
        encoding_weights_matrix = encoding_matrix,
        parcel_names = parcel_names,
        k_conn_task_pos = 2,
        k_conn_task_neg = 1,
        similarity_method = "pearson"
      )),
      regexp = "Found 1 parcel\\(s\\) with zero variance"
    )
    expect_s4_class(W_task, "dgCMatrix")
  })
  
  it("handles k_conn_pos=0 and k_conn_neg=0 for encoding", {
    encoding_matrix <- matrix(rnorm(V_p * N_features), nrow = V_p, ncol = N_features)
    rownames(encoding_matrix) <- parcel_names
    W_task <- compute_W_task_from_encoding(
      encoding_weights_matrix = encoding_matrix,
      parcel_names = parcel_names,
      k_conn_task_pos = 0,
      k_conn_task_neg = 0,
      similarity_method = "pearson"
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p, V_p))
    expect_equal(length(W_task@x), 0) 
  })
  
  it("handles input with 0 parcels (V_p=0) for encoding", {
    encoding_matrix_empty_vp <- matrix(nrow = 0, ncol = N_features)
    W_task <- compute_W_task_from_encoding(
        encoding_weights_matrix = encoding_matrix_empty_vp,
        parcel_names = character(0),
        k_conn_task_pos = 1, k_conn_task_neg = 1
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(0,0))
  })
  
  it("handles input with 0 features (N_features=0) for encoding", {
    encoding_matrix_empty_nf <- matrix(nrow = V_p, ncol = 0)
    rownames(encoding_matrix_empty_nf) <- parcel_names
    
    expect_warning(
      W_task <- compute_W_task_from_encoding(
          encoding_weights_matrix = encoding_matrix_empty_nf,
          parcel_names = parcel_names,
          k_conn_task_pos = 1, k_conn_task_neg = 1
      ),
      regexp = "encoding_weights_matrix has zero features"
    )
    expect_s4_class(W_task, "dgCMatrix")
    expect_equal(dim(W_task), c(V_p,V_p))
    expect_equal(length(W_task@x), 0) # Graph should be empty
  })

  it("errors for invalid k_conn_task_pos or k_conn_task_neg for encoding", {
    encoding_matrix <- matrix(rnorm(V_p * N_features), nrow = V_p, ncol = N_features)
    rownames(encoding_matrix) <- parcel_names
    expect_error(
      compute_W_task_from_encoding(
        encoding_weights_matrix = encoding_matrix,
        parcel_names = parcel_names,
        k_conn_task_pos = -1,
        k_conn_task_neg = 1
      ),
      "k_conn_task_pos"
    )
    expect_error(
      compute_W_task_from_encoding(
        encoding_weights_matrix = encoding_matrix,
        parcel_names = parcel_names,
        k_conn_task_pos = 1,
        k_conn_task_neg = -1
      ),
      "k_conn_task_neg"
    )
  })

}) # end describe compute_W_task_from_encoding 