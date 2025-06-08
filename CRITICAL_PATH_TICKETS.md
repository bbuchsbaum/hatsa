# CRITICAL PATH IMPLEMENTATION TICKETS

**Priority: URGENT** - These 5 functions must be implemented before any other HATSA functionality can work.

---

## **TICKET CP-001: `compute_subject_connectivity_graph_sparse()`**

**File**: `R/spectral_graph_construction.R`  
**Status**: ❌ Function signature exists, no implementation  
**Estimated Time**: 1-2 days  

### **Function Signature**
```r
compute_subject_connectivity_graph_sparse <- function(X_subject, parcel_names, k_conn_pos, k_conn_neg, use_dtw = FALSE)
```

### **Purpose**
Build sparse correlation graphs from subject time series data. This is the foundational function that creates the connectivity matrices used throughout HATSA.

### **Input**
- `X_subject`: Numeric matrix (T x V_p) - time series for one subject
- `parcel_names`: Character vector of length V_p - parcel identifiers  
- `k_conn_pos`: Integer - number of strongest positive connections to retain per parcel
- `k_conn_neg`: Integer - number of strongest negative connections to retain per parcel
- `use_dtw`: Logical - whether to use Dynamic Time Warping (not implemented yet)

### **Expected Output**
Sparse symmetric matrix (dgCMatrix, V_p x V_p) with:
- Sparsified correlation values
- Z-scored edge weights (non-zero elements)
- Symmetric structure
- Zero diagonal

### **Algorithm**
1. Compute full correlation matrix: `cor(X_subject)` (V_p x V_p)
2. Set diagonal to 0
3. Apply k-NN sparsification separately for positive and negative correlations
4. Symmetrize: `(W + t(W)) / 2` 
5. Z-score non-zero edge weights using existing `zscore_nonzero_sparse()` function
6. Return as sparse dgCMatrix

### **Dependencies**
- `Matrix` package (already imported)
- `zscore_nonzero_sparse()` function (already implemented)

### **Tests Ready**
✅ `tests/testthat/test-spectral_graph_construction.R` lines 15-200

---

## **TICKET CP-002: `compute_graph_laplacian_sparse()`**

**File**: `R/spectral_graph_construction.R`  
**Status**: ❌ Function signature exists, no implementation  
**Estimated Time**: 1 day  

### **Function Signature**
```r
compute_graph_laplacian_sparse <- function(W_sparse, alpha = 0.93, degree_type = "abs")
```

### **Purpose**
Compute alpha-lazy random-walk normalized Laplacian from sparse adjacency matrix.

### **Input**
- `W_sparse`: Sparse symmetric matrix (dgCMatrix) - adjacency/weight matrix
- `alpha`: Numeric in [0,1] - laziness parameter (default 0.93)
- `degree_type`: Character - how to compute degrees ("abs", "positive", "signed")

### **Mathematical Formula**
```
L = I - alpha * D^(-1) * W
where D is degree matrix based on degree_type
```

### **Expected Output**
Sparse symmetric matrix (dgCMatrix, V_p x V_p):
- Alpha-lazy random-walk normalized Laplacian
- Symmetric: `(L + t(L))/2`
- Proper handling of zero-degree nodes

### **Algorithm**
1. Compute degree vector based on `degree_type`:
   - "abs": `rowSums(abs(W_sparse))`
   - "positive": `rowSums(pmax(W_sparse, 0))`
   - "signed": `rowSums(W_sparse)`
2. Handle zero degrees: set to small epsilon or special case
3. Create degree matrix: `D_inv = Diagonal(x = 1/degrees)`
4. Compute: `L_rw = I - alpha * D_inv %*% W_sparse`
5. Symmetrize: `L = (L_rw + t(L_rw)) / 2`
6. Return as sparse dgCMatrix

### **Dependencies**
- `Matrix` package (already imported)

### **Tests Ready**
✅ `tests/testthat/test-spectral_graph_construction.R` lines 295-380

---

## **TICKET CP-003: `compute_spectral_sketch_sparse()`**

**File**: `R/spectral_graph_construction.R`  
**Status**: ❌ Function signature exists, no implementation  
**Estimated Time**: 1-2 days  

### **Function Signature**
```r
compute_spectral_sketch_sparse <- function(L_sparse, k, eigenvalue_tol = 1e-8)
```

### **Purpose**
Compute sparse eigendecomposition with DC component filtering for HATSA spectral sketches.

### **Input**
- `L_sparse`: Sparse symmetric matrix (dgCMatrix) - Laplacian matrix
- `k`: Integer - number of eigenvectors to return
- `eigenvalue_tol`: Numeric - threshold for filtering trivial eigenvectors

### **Expected Output**
List with two elements:
- `vectors`: Numeric matrix (V_p x k) - eigenvectors  
- `values`: Numeric vector (length k) - eigenvalues

### **Algorithm**
1. Determine how many eigenvectors to request (need k + buffer for filtering)
2. Use `RSpectra::eigs_sym()` to compute smallest eigenvalues/vectors
3. Filter out trivial eigenvectors where `eigenvalue < eigenvalue_tol`
4. Take first `k` non-trivial eigenvectors
5. Return as list with `vectors` and `values`

### **Edge Cases**
- Handle `k = 0`: return empty matrices
- Handle insufficient non-trivial eigenvectors: throw informative error
- Handle numerical issues in eigendecomposition

### **Dependencies**
- `RSpectra` package (already imported)

### **Tests Ready**
✅ `tests/testthat/test-spectral_graph_construction.R` lines 400-590

---

## **TICKET CP-004: `misalign_deg()`**

**File**: `R/metrics.R`  
**Status**: ❌ Exported in NAMESPACE but no implementation  
**Estimated Time**: 0.5 days  

### **Function Signature**
```r
misalign_deg <- function(R1, R2, method = "geodesic")
```

### **Purpose**
Compute geodesic distance between two rotation matrices on SO(k) manifold, returned in degrees.

### **Input**
- `R1`, `R2`: Numeric matrices (k x k) - rotation matrices in SO(k)
- `method`: Character - distance computation method ("geodesic")

### **Mathematical Formula**
```
d = acos((trace(t(R1) %*% R2) - 1) / 2) * 180/π
```

### **Expected Output**
Numeric scalar - geodesic distance in degrees [0, 180]

### **Algorithm**
1. Validate inputs are square matrices of same dimension
2. Compute trace: `tr = sum(diag(t(R1) %*% R2))`
3. Handle numerical edge cases:
   - If `tr > k + 1`: clamp to `k + 1`
   - If `tr < k - 1`: clamp to `k - 1`
4. Compute: `d_rad = acos((tr - 1) / 2)`
5. Convert to degrees: `d_deg = d_rad * 180 / pi`
6. Return scalar value

### **Error Handling**
- Check matrix dimensions match
- Handle numerical issues in `acos()`
- Validate rotation matrix properties (optional)

### **Dependencies**
- Base R only

### **Tests Ready**
✅ `tests/testthat/test-metrics.R` lines 1-150

---

## **TICKET CP-005: `solve_gev_laplacian_primme()`**

**File**: `R/gev_helpers.R`  
**Status**: ❌ Exported in NAMESPACE but no implementation  
**Estimated Time**: 1 day  

### **Function Signature**
```r
solve_gev_laplacian_primme <- function(A, B, k_request, lambda_max_thresh = 0.8, epsilon_reg_B = 1e-6, tol = 1e-8)
```

### **Purpose**
Solve generalized eigenvalue problem A*v = λ*B*v for task HATSA GEV patches.

### **Input**
- `A`: Sparse symmetric matrix (dgCMatrix) - typically task Laplacian
- `B`: Sparse symmetric matrix (dgCMatrix) - typically connectivity Laplacian  
- `k_request`: Integer - number of eigenpairs to compute
- `lambda_max_thresh`: Numeric - maximum eigenvalue to retain
- `epsilon_reg_B`: Numeric - regularization for B matrix
- `tol`: Numeric - convergence tolerance

### **Expected Output**
List with:
- `vectors`: Numeric matrix (V_p x k_filtered) - filtered eigenvectors
- `values`: Numeric vector (length k_filtered) - filtered eigenvalues  
- `n_converged`: Integer - number of eigenvalues that converged
- `n_filtered`: Integer - number retained after filtering

### **Algorithm**
1. Regularize B matrix: `B_reg = B + epsilon_reg_B * I`
2. Use `PRIMME::eigs_sym(A = A, B = B_reg, NEig = k_request, which = "SM")`
3. Filter results: keep only eigenvalues with `abs(lambda) < lambda_max_thresh`
4. Return list with filtered vectors, values, and diagnostic counts

### **Error Handling**
- Check matrix dimensions match
- Handle convergence failures
- Validate eigenvalue filtering results

### **Dependencies**
- `PRIMME` package (already imported)
- `Matrix` package (already imported)

### **Tests Ready**
✅ `tests/testthat/test-gev.R` lines 30-100

---

## **🚀 IMPLEMENTATION ORDER**

### **Day 1-2**: Foundation
1. **CP-001** (connectivity graphs) - Most critical, everything depends on this
2. **CP-002** (Laplacian) - Direct dependency of spectral sketches

### **Day 3-4**: Core Analysis  
3. **CP-003** (spectral sketches) - Core eigendecomposition functionality
4. **CP-004** (misalign_deg) - Needed for validation metrics

### **Day 5**: Advanced Features
5. **CP-005** (GEV solver) - Enables task HATSA GEV patches

---

## **✅ VALIDATION WORKFLOW**

After implementing each function:

```bash
# Test specific function
R -e "devtools::test_file('tests/testthat/test-spectral_graph_construction.R')"

# Test integration  
R -e "devtools::test_file('tests/testthat/test-hatsa_core_functionality.R')"

# Full test suite (after all 5 implemented)
R -e "devtools::test()"
```

**Success Criteria**: 
- All function-specific tests pass
- `run_hatsa_core()` executes without errors
- Integration tests pass
- Task HATSA and voxel projection work end-to-end

---

**🎯 Start with CP-001 and work through the list sequentially. Each function builds on the previous ones.** 