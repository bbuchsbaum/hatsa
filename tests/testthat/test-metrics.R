describe("metrics - misalign_deg", {

test_that("misalign_deg basic functionality and known angles", {
  # Identity matrices
  R_ident2 <- diag(2)
  expect_equal(misalign_deg(R_ident2, R_ident2), 0)

  R_ident3 <- diag(3)
  expect_equal(misalign_deg(R_ident3, R_ident3), 0)

  # Known 2D rotations
  theta_30_rad <- pi/6
  R_2d_30deg <- matrix(c(cos(theta_30_rad), -sin(theta_30_rad),
                         sin(theta_30_rad),  cos(theta_30_rad)), 2, 2)
  expect_equal(misalign_deg(R_2d_30deg, R_ident2), 30, tolerance = 1e-7)
  expect_equal(misalign_deg(R_ident2, R_2d_30deg), 30, tolerance = 1e-7) # Symmetry

  theta_45_rad <- pi/4
  R_2d_45deg <- matrix(c(cos(theta_45_rad), -sin(theta_45_rad),
                         sin(theta_45_rad),  cos(theta_45_rad)), 2, 2)
  expect_equal(misalign_deg(R_2d_45deg, R_ident2), 45, tolerance = 1e-7)

  theta_90_rad <- pi/2
  R_2d_90deg <- matrix(c(cos(theta_90_rad), -sin(theta_90_rad),
                         sin(theta_90_rad),  cos(theta_90_rad)), 2, 2)
  expect_equal(misalign_deg(R_2d_90deg, R_ident2), 90, tolerance = 1e-7)

  # Known 3D rotation (around Z-axis)
  theta_60_rad <- pi/3
  R_3d_z_60deg <- matrix(c(cos(theta_60_rad), -sin(theta_60_rad), 0,
                           sin(theta_60_rad),  cos(theta_60_rad), 0,
                           0,                  0,                 1), 3, 3, byrow = TRUE)
  expect_equal(misalign_deg(R_3d_z_60deg, R_ident3), 60, tolerance = 1e-7)
  expect_equal(misalign_deg(R_ident3, R_3d_z_60deg), 60, tolerance = 1e-7) # Symmetry
  
  # Known 3D rotation (around Y-axis)
  theta_90_rad_y <- pi/2
  R_3d_y_90deg <- matrix(c(cos(theta_90_rad_y), 0, sin(theta_90_rad_y),
                           0,                   1, 0,
                           -sin(theta_90_rad_y),0, cos(theta_90_rad_y)), 3, 3, byrow = TRUE)
  expect_equal(misalign_deg(R_3d_y_90deg, R_ident3), 90, tolerance = 1e-7)

  # Test composition: R_AB, R_BC, expect dist(R_AC, I) related to dist(R_AB,I) + dist(R_BC,I)
  # Not strictly additive for large rotations, but for small ones it's close.
  # For SO(3), angle(R1*R2) is not angle(R1)+angle(R2) in general.
  # Geodesic from R_A to R_C via R_B: d(A,C) <= d(A,B) + d(B,C)
  R_A <- R_2d_30deg
  R_B <- R_2d_45deg %*% R_A # R_B is R_A rotated by 45 deg more
  # Misalignment between R_B and R_A should be 45 deg
  expect_equal(misalign_deg(R_B, R_A), 45, tolerance = 1e-7)
  
})

test_that("misalign_deg input validation and error handling", {
  R_valid2 <- diag(2)
  R_valid3 <- diag(3)

  # Non-matrix inputs
  expect_warning(res <- misalign_deg(as.vector(R_valid2), R_valid2), "Inputs must be matrices.")
  expect_true(is.na(res))
  expect_warning(res <- misalign_deg(R_valid2, "not_a_matrix"), "Inputs must be matrices.")
  expect_true(is.na(res))

  # Mismatched dimensions
  expect_warning(res <- misalign_deg(R_valid3, R_valid2), "Rotation matrices must have the same dimensions.")
  expect_true(is.na(res))

  # Non-square matrices
  R_nonsquare1 <- matrix(1:6, 2, 3)
  R_nonsquare2 <- matrix(1:6, 2, 3)
  expect_warning(res <- misalign_deg(R_nonsquare1, R_nonsquare2), "Matrices must be square.")
  expect_true(is.na(res))
  
  # Case where M_rel is identity (already covered by basic tests but good to be explicit)
  R1 <- diag(3)
  expect_equal(misalign_deg(R1,R1), 0)
  
})

# It's hard to reliably trigger the expm::logm failure for a valid SO(k) matrix in a simple test
# without creating a matrix that is numerically very problematic (e.g., extremely close to singular,
# or far from orthogonal, which misalign_deg doesn't strictly check before crossprod).
# The safe_logm_internal is designed to catch errors from expm::logm.
# We can test the fallback if `expm` is not available by temporarily mocking `requireNamespace`.

test_that("misalign_deg fallback mechanism when expm is not available", {
  # Mock requireNamespace to simulate expm not being installed
  mockery::stub(misalign_deg, 'requireNamespace', function(pkg, ...) if(pkg=='expm') FALSE else TRUE)
  
  R1 <- diag(3)
  theta <- pi/4 # 45 degrees
  R2 <- matrix(c(cos(theta), -sin(theta), 0,
                 sin(theta),  cos(theta), 0,
                 0,           0,          1), nrow=3, byrow=TRUE)
  
  # For 3x3, the fallback (acos((tr(M_rel)-1)/2)) should be accurate
  expect_warning(
    angle_fallback <- misalign_deg(R2, R1, method="geodesic"), 
    "Package 'expm' not available for geodesic method. Using fallback trace-based calculation."
  )
  expect_equal(angle_fallback, 45, tolerance = 1e-7)
  
  # For 2x2, the fallback trace is trace(R)/2 for cos(theta)
  R1_2d <- diag(2)
  R2_2d <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2,2)
  
  # Capture all warnings for 2D case
  warnings_2d <- capture_warnings(
    angle_fallback_2d <- misalign_deg(R2_2d, R1_2d, method="geodesic")
  )
  expect_true(any(grepl("Package 'expm' not available", warnings_2d)))
  expect_true(any(grepl("Fallback trace-based angle is most accurate for 3x3 rotations", warnings_2d)))
  expect_equal(angle_fallback_2d, 45, tolerance = 1e-7)
})

test_that("misalign_deg specific k_dim warnings for fallback", {
  # Mock requireNamespace to force fallback
  mockery::stub(misalign_deg, 'requireNamespace', function(pkg, ...) if(pkg=='expm') FALSE else TRUE)
  
  # k=3 should not warn about k_dim != 3 for fallback accuracy
  R1_3d <- diag(3)
  theta <- pi/6
  R2_3d <- matrix(c(cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1), 3,3, byrow=TRUE)
  expect_warning(
    misalign_deg(R2_3d, R1_3d, method="geodesic"), # This warns about expm missing
    "Package 'expm' not available for geodesic method" 
  )
  
  # k=2 should warn about k_dim != 3 for fallback accuracy
  R1_2d <- diag(2)
  R2_2d <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2,2)
  
  # Check that we get the k_dim warning when using the fallback
  # We need to capture all warnings and check that the k_dim warning is present
  warnings_captured <- capture_warnings(
    misalign_deg(R2_2d, R1_2d, method="geodesic")
  )
  expect_true(any(grepl("Package 'expm' not available", warnings_captured)))
  expect_true(any(grepl("Fallback trace-based angle is most accurate for 3x3 rotations", warnings_captured)))
  
  # k=4 should warn about k_dim != 3 for fallback accuracy
  R1_4d <- diag(4)
  # Simple rotation in first 2 dims
  R2_4d_block <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2,2)
  R2_4d <- diag(4)
  R2_4d[1:2,1:2] <- R2_4d_block
  
  # Check that we get the k_dim warning when using the fallback
  warnings_captured_4d <- capture_warnings(
    misalign_deg(R2_4d, R1_4d, method="geodesic")
  )
  expect_true(any(grepl("Package 'expm' not available", warnings_captured_4d)))
  expect_true(any(grepl("Fallback trace-based angle is most accurate for 3x3 rotations", warnings_captured_4d)))
}) })
