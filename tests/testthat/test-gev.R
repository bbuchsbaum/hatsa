library(testthat)
library(Matrix)

# Assuming solve_gev_laplacian_primme is available from hatsa package (e.g. via devtools::load_all())

describe("GEV Helper Functions: solve_gev_laplacian_primme", {

test_that("solve_gev_laplacian_primme handles block-diagonal Laplacians and filters correctly", {
  skip_if_not_installed("PRIMME")

  # 1. Construct L1 (3x3 path graph Laplacian)
  W1 <- Matrix::sparseMatrix(i = c(1, 2, 2, 3), j = c(2, 1, 3, 2), x = 1, dims = c(3, 3))
  D1 <- Matrix::Diagonal(3, x = Matrix::rowSums(W1))
  L1 <- D1 - W1

  # 2. Construct L2 (2x2 path graph Laplacian)
  W2 <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = 1, dims = c(2, 2))
  D2 <- Matrix::Diagonal(2, x = Matrix::rowSums(W2))
  L2 <- D2 - W2

  # 3. Combine into a block-diagonal A
  A_mat <- Matrix::bdiag(L1, L2)
  V_p <- nrow(A_mat) # Should be 5

  # 4. B matrix (Identity)
  B_mat <- Matrix::Diagonal(V_p)

  # 5. Parameters for the test
  k_request_test <- V_p  # Request all eigenpairs
  lambda_max_thresh_test <- 1.5
  epsilon_reg_B_test <- 1e-6 # Matches default in solve_gev_laplacian_primme

  # Eigenvalues of A_mat are {0, 1, 3} from L1 and {0, 2} from L2.
  # Combined and sorted: {0, 0, 1, 2, 3}
  # PRIMME solves A v = lambda_primme B_reg v, where B_reg = (1 + epsilon_reg_B_test) I
  # So, lambda_primme = lambda_A / (1 + epsilon_reg_B_test)
  # Expected PRIMME values (approx): 0, 0, 1/(1+eps), 2/(1+eps), 3/(1+eps)
  # Approx: 0, 0, 0.999999, 1.999998, 2.999997

  # Filtering abs(lambda_primme) < 1.5 should keep the first three: ~0, ~0, ~1.
  expected_n_filtered <- 3
  expected_values_approx <- c(0, 0, 1) / (1 + epsilon_reg_B_test)

  result <- suppressMessages(solve_gev_laplacian_primme(
    A = A_mat,
    B = B_mat,
    k_request = k_request_test,
    lambda_max_thresh = lambda_max_thresh_test,
    epsilon_reg_B = epsilon_reg_B_test,
    tol = 1e-8 # Use a reasonable tolerance for PRIMME
  ))

  # Assertions
  expect_true(is.list(result))
  expect_named(result, c("vectors", "values", "n_converged", "n_filtered", "primme_stats"))

  # Check number of converged eigenvalues (PRIMME might struggle with exact multiplicity)
  # Allow for slight variation if k_request_test = V_p
  expect_gte(result$n_converged, k_request_test - 2) # e.g. at least 3 for V_p=5 if zeros are tricky
  expect_lte(result$n_converged, k_request_test)

  # Check filtered results
  expect_equal(result$n_filtered, expected_n_filtered)
  expect_equal(ncol(result$vectors), expected_n_filtered)
  expect_length(result$values, expected_n_filtered)

  # Verify that all returned values satisfy the threshold
  if (result$n_filtered > 0) {
    expect_true(all(abs(result$values) < lambda_max_thresh_test))
  }

  # Verify the actual filtered eigenvalues (sorted by absolute magnitude in the function)
  # The function sorts by abs value, so for {0,0, ~1}, order is preserved.
  if (result$n_filtered == expected_n_filtered) {
    expect_equal(result$values, expected_values_approx, tolerance = 1e-5)
  }

  # Test with a different threshold that should filter more
  lambda_max_thresh_tight <- 0.5
  expected_n_filtered_tight <- 2 # Should only keep the two ~0 eigenvalues
  expected_values_tight_approx <- c(0,0) / (1+epsilon_reg_B_test)
  
  result_tight <- suppressMessages(solve_gev_laplacian_primme(
    A = A_mat,
    B = B_mat,
    k_request = k_request_test,
    lambda_max_thresh = lambda_max_thresh_tight,
    epsilon_reg_B = epsilon_reg_B_test,
    tol = 1e-8
  ))
  
  expect_equal(result_tight$n_filtered, expected_n_filtered_tight)
  if (result_tight$n_filtered == expected_n_filtered_tight) {
      expect_equal(result_tight$values, expected_values_tight_approx, tolerance = 1e-5)
      if (result_tight$n_filtered > 0) {
        expect_true(all(abs(result_tight$values) < lambda_max_thresh_tight))
      }
  }
})

})
