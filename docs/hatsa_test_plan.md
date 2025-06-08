
**Granular Testing Tickets for R `hatsa` Package**

**Module: `spectral_graph_construction.R`**

*   **TCK-SGC-001: `sparse_correlation_r` - Basic Sparsification & Output**
    *   **Synthetic Data:** Small matrix `X_ts` (e.g., 50x10) with known correlation structure.
    *   **Test:** Call `sparse_correlation_r` with small `k_conn_pos=2`, `k_conn_neg=1`.
    *   **Assertions:**
        *   `expect_s3_class(output, "dgCMatrix")`.
        *   `expect_equal(dim(output), c(10, 10))`.
        *   Check manually calculated sparsity for a few nodes matches output non-zero counts.
        *   `expect_true(Matrix::isSymmetric(output))` (after z-scoring step).
        *   `expect_equal(Matrix::diag(output), rep(0, 10))`.

*   **TCK-SGC-002: `sparse_correlation_r` - Edge Weight Correctness**
    *   **Synthetic Data:** `X_ts` where specific pairs have known high positive, high negative, and near-zero correlations.
    *   **Test:** Call `sparse_correlation_r` with `k_conn` chosen to select/exclude specific edges.
    *   **Assertions:**
        *   Check `output@x` (non-zero values) match the *z-scored* expected correlation values for selected edges.
        *   Verify edges near zero are correctly excluded.

*   **TCK-SGC-003: `sparse_correlation_r` - Handling Zero Variance Columns**
    *   **Synthetic Data:** `X_ts` where one or more columns are constant.
    *   **Test:** Call `sparse_correlation_r`.
    *   **Assertions:**
        *   `expect_message()` containing "zero variance".
        *   Verify correlations involving the constant column are 0 in the output sparse matrix.

*   **TCK-SGC-004: `sparse_correlation_r` - Symmetrization Logic**
    *   **Synthetic Data:** Create a small *directed* sparse matrix `W_dir` manually where `W_dir[i,j]` exists but `W_dir[j,i]` doesn't, or where they have different values.
    *   **Test:** Manually apply the internal symmetrization logic `(W_dir + t(W_dir)) / pmax(1, (W_dir != 0) + (t(W_dir) != 0))` followed by `forceSymmetric` and `zscore_nonzero_sparse`.
    *   **Assertions:** Check specific `[i,j]` and `[j,i]` entries match expected symmetric, z-scored value based on the averaging rule.

*   **TCK-SGC-005: `compute_graph_laplacian_sparse` - Correctness (α-lazy L_rw)**
    *   **Synthetic Data:** Small, known sparse symmetric adjacency `W_sparse` (e.g., 5x5). Manually calculate degrees `D`, `D_inv`, and `L = I - alpha * D_inv %*% W_sparse`. Manually symmetrize `(L + t(L))/2`.
    *   **Test:** Call `compute_graph_laplacian_sparse(W_sparse, alpha=0.93)`.
    *   **Assertions:** `expect_equal(output_matrix, manual_L_sym, tolerance=1e-8)`. Check output is sparse and symmetric.

*   **TCK-SGC-006: `compute_graph_laplacian_sparse` - Handling Zero Degrees**
    *   **Synthetic Data:** `W_sparse` containing one or more isolated nodes (zero degree).
    *   **Test:** Call `compute_graph_laplacian_sparse`.
    *   **Assertions:**
        *   No errors/warnings about division by zero.
        *   Verify corresponding rows/columns in the output `L` are zero (or identity on diagonal if `I` dominates).

*   **TCK-SGC-007: `compute_spectral_sketch_sparse` - Basic Correctness & Dimensions**
    *   **Synthetic Data:** Known graph Laplacian `L` (e.g., from SBM in T-1) with analytically known first few eigenvectors/values.
    *   **Test:** Call `compute_spectral_sketch_sparse(L, k=2)`.
    *   **Assertions:**
        *   `expect_true(is.list(output))` with names "vectors" and "values".
        *   `expect_equal(dim(output$vectors), c(nrow(L), 2))`.
        *   Compare `output$vectors` to analytical `U_analytic[,1:2]` (allow for sign flips, `expect_fro_error(abs(output$vectors), abs(U_analytic), tol=1e-5)`).
        *   Compare `output$values` to analytical eigenvalues `lambda_analytic[1:2]`.

*   **TCK-SGC-008: `compute_spectral_sketch_sparse` - Trivial Eigenvector Handling**
    *   **Synthetic Data:** Laplacian `L` from a connected graph (guaranteed zero eigenvalue).
    *   **Test:** Call `compute_spectral_sketch_sparse(L, k=3)`. Check eigenvalues returned.
    *   **Assertions:**
        *   Verify the *smallest* eigenvalue returned is `> eigenvalue_tol` (default 1e-8).
        *   Verify the eigenvector corresponding to the true zero eigenvalue was correctly discarded.

*   **TCK-SGC-009: `compute_spectral_sketch_sparse` - Rank Deficiency / Insufficient Eigenvectors**
    *   **Synthetic Data:** Graph `L` designed to have fewer than `k+1` non-zero eigenvalues (e.g., multiple connected components).
    *   **Test:** Call `compute_spectral_sketch_sparse(L, k=5)` where only, say, 3 informative eigenvectors exist.
    *   **Assertions:** `expect_error()` with informative message about rank deficiency.

*   **TCK-SGC-010: `compute_spectral_sketch_sparse` - Edge Cases `k=0`, `k=1`**
    *   **Synthetic Data:** Standard `L`.
    *   **Test:** Call with `k=0`.
    *   **Assertions:** `expect_equal(ncol(output$vectors), 0)`, `expect_length(output$values, 0)`.
    *   **Test:** Call with `k=1`.
    *   **Assertions:** `expect_equal(ncol(output$vectors), 1)`, `expect_length(output$values, 1)`.

**Module: `procrustes_alignment.R`**

*   **TCK-GPA-001: `solve_procrustes_rotation` - Correctness & Properties**
    *   **Synthetic Data:** Create `A_orig` (e.g., 10x5 random matrix) and known rotation `R_true` (e.g., using `pracma::randortho(5)`). Calculate `T_target = A_orig %*% R_true`.
    *   **Test:** Call `R_est = solve_procrustes_rotation(A_orig, T_target)`.
    *   **Assertions:**
        *   `expect_equal(R_est, R_true, tolerance=1e-7)`.
        *   `expect_equal(crossprod(R_est), diag(5), tolerance=1e-7)` (orthogonality).
        *   `expect_equal(det(R_est), 1.0, tolerance=1e-7)` (proper rotation).

*   **TCK-GPA-002: `solve_procrustes_rotation` - Reflection Handling**
    *   **Synthetic Data:** As above, but force `R_true` to be a reflection (`det(R_true)=-1`). Calculate `T_target = A_orig %*% R_true`.
    *   **Test:** Call `R_est = solve_procrustes_rotation(A_orig, T_target)`.
    *   **Assertions:** Check `det(R_est)` is close to `+1`, confirming reflection correction logic worked. `expect_false(isTRUE(all.equal(R_est, R_true)))` but `expect_true(norm(A_orig %*% R_est - T_target, 'F')` is minimized).

*   **TCK-GPA-003: `solve_procrustes_rotation` - Edge Cases `k=0`, `k=1`**
    *   **Test:** Call with `k=0` inputs (`A`: m x 0, `T`: m x 0).
    *   **Assertions:** `expect_equal(output, matrix(0,0,0))`.
    *   **Test:** Call with `k=1` inputs (`A`: m x 1, `T`: m x 1).
    *   **Assertions:** `expect_equal(dim(output), c(1,1))`. Output should be `matrix(1)` or `matrix(-1)` depending on sign correlation, ensure `det(output)==1` (trivial for 1x1).

*   **TCK-GPA-004: `perform_gpa_refinement` - Basic Convergence**
    *   **Synthetic Data (T-2):** Generate `U_analytic` (e.g., 60x40 SBM eigenvectors). Create 10 subjects: `U_orig_i` = `U_analytic`. Apply random `R_true_i ∈ SO(40)` to get `A_target_i = U_orig_i[1:10,] %*% R_true_i`. (Here target is known). For test, use `A_originals_list` derived from `U_analytic[1:10,]`.
    *   **Test:** Call `perform_gpa_refinement(A_originals_list, n_refine=5, k=40)`. Let `R_final = output$R_final_list`.
    *   **Assertions:**
        *   Verify rotations converge (e.g., `norm(R_iter_n - R_iter_{n-1})` decreases).
        *   Check anchor residual `mean(sapply(1:10, function(i) norm(A_originals_list[[i]] %*% R_final[[i]] - output$T_anchor_final, 'F')))` is small (< 1e-4).
        *   *Note: This doesn't compare to `R_true_i` directly, just checks internal consistency.*

*   **TCK-GPA-005: `perform_gpa_refinement` - Handling `n_refine=0`**
    *   **Synthetic Data:** `A_originals_list` as above.
    *   **Test:** Call `perform_gpa_refinement(A_originals_list, n_refine=0, k=40)`.
    *   **Assertions:** `expect_true(all(sapply(output$R_final_list, function(R) isTRUE(all.equal(R, diag(40))))))` (rotations are identity). `expect_true(norm(output$T_anchor_final - Reduce('+', A_originals_list)/10) < 1e-8)` (template is simple mean).

**Module: `run_hatsa_core_algorithm.R` & `hatsa_projector.R` (End-to-End & Object)**

*   **TCK-CORE-001: `run_hatsa_core` - End-to-End Recovery (T-2 Combined)**
    *   **Synthetic Data (T-2):** Generate `U_analytic` (e.g., 60x40 SBM). Create `X_ts_list` for 10 subjects such that their `U_orig_i` *would be* `U_analytic`. Apply known random `R_true_i` to `U_analytic` to get target aligned `U_aligned_target_i`. *Challenge: Creating `X_ts` that yields exact `U_analytic` is hard. Alternative: Use `U_analytic` directly in a mocked `compute_spectral_sketch` OR use T-2 Procrustes test setup.*
    *   **Test (Using T-2 Setup):** Generate `U_analytic`, apply `R_true_i` -> `U_rotated_i`. Use `A_orig = U_analytic[anchors,]` and `T_target = mean(U_rotated_i[anchors,])`. Run Procrustes `perform_gpa_refinement` starting with `A_orig` and target `T_target`.
    *   **Assertions:** Mean angular error `mean(sapply(1:10, function(i) acos((sum(diag(t(R_est_i) %*% R_true_i)) - 1)/2)))` < 0.05 rad. *Self-correction: This tests GPA part well, less so graph->sketch part.*
    *   **Test (Full `run_hatsa_core`):** Use SBM to generate *graph `W`*. Create `X_ts` list where `cor(X_ts)` approximates `W`. Run `run_hatsa_core`. *Challenge: Ground truth `R_i` is unknown.* Assert low final anchor residual, reasonable `κ`.

*   **TCK-CORE-002: `run_hatsa_core` - Input Validation Passthrough**
    *   **Test:** Call `run_hatsa_core` with invalid inputs (e.g., `k > m`, non-integer `k`, `anchor_indices` out of bounds, mismatched `V_p`).
    *   **Assertions:** `expect_error()` with messages matching those from `validate_hatsa_inputs`.

*   **TCK-CORE-003: `hatsa_projector` Object Structure**
    *   **Test:** Run `run_hatsa_core` successfully.
    *   **Assertions:** Check output object contains all documented fields (`v`, `s`, `sdev`, `preproc`, `block_indices`, `R_final_list`, `U_original_list`, `Lambda_original_list`, `T_anchor_final`, `parameters`, `method`) with correct types and basic dimension consistency.

*   **TCK-CORE-004: S3 Methods (`print`, `summary`, `coef`, `scores`, etc.)**
    *   **Test:** Call each S3 method on a valid `hatsa_projector` object.
    *   **Assertions:**
        *   `expect_output()` for `print`.
        *   `expect_s3_class(summary(obj), "summary.hatsa_projector")`. `expect_output(print(summary(obj)))`.
        *   Check dimensions/types returned by `coef`, `scores`, `sdev`, `block_indices` match expectations.

*   **TCK-CORE-005: `predict.hatsa_projector` - Basic Functionality**
    *   **Test:** Create `hatsa_obj`. Generate `newdata_list`. Call `predict(hatsa_obj, newdata_list)`.
    *   **Assertions:** Check output is list of matrices with correct dimensions (`V_p x k`). Test error handling for mismatched `V_p`.

*   **TCK-CORE-006: `project_block.hatsa_projector` - Functionality**
    *   **Test:** Call `project_block` with `newdata=NULL` (retrieve stored) and with `newdata` provided (predict new).
    *   **Assertions:** Check output dimensions (`V_p x k`). Test error handling for invalid `block` index.

**Module: `voxel_projection.R` (Nyström)**

*   **TCK-NYSTROM-001: `compute_voxel_basis_nystrom` - Identity Check (T-6)**
    *   **Synthetic Data:** `V_p=50`, `k=10`. `parcel_coords` = random. `voxel_coords` = exactly `parcel_coords`. `U_orig_parcel` = random `V_p x k`. `Lambda_orig_parcel` = random `k`.
    *   **Test:** Call `compute_voxel_basis_nystrom` with `n_nearest_parcels=1`, `kernel_sigma=1e-6` (approx delta function), `row_normalize_W=FALSE` (or adjust formula if needed).
    *   **Assertions:** `expect_equal(output_Phi_voxel, U_orig_parcel %*% diag(1/Lambda_orig_parcel), tolerance=1e-6)`.

*   **TCK-NYSTROM-002: `compute_voxel_basis_nystrom` - Smoothness/Locality (T-7)**
    *   **Synthetic Data:** Points on sphere. Parcels at fixed vertices, voxels nearby. Calculate analytical expected decay `exp(-dist^2 / 2*sigma^2)`.
    *   **Test:** Call `compute_voxel_basis_nystrom`. Project a simple spatial pattern (e.g., linear gradient) defined on parcels into voxel space using `output_Phi_voxel`.
    *   **Assertions:** Check `cor(projected_voxel_pattern, expected_smooth_pattern) > 0.97`.

*   **TCK-NYSTROM-003: `compute_voxel_basis_nystrom` - Coordinate Checks**
    *   **Synthetic Data:** `voxel_coords` and `parcel_coords` with deliberately mismatched scales (e.g., mm vs m) or disjoint ranges.
    *   **Test:** Call `compute_voxel_basis_nystrom` (or just `.validate_coordinate_inputs`).
    *   **Assertions:** `expect_message()` containing warnings about scale or range mismatch (when run interactively).

*   **TCK-NYSTROM-004: `project_voxels.hatsa_projector` - End-to-End**
    *   **Test:** Create `hatsa_obj`. Create mock `voxel_timeseries_list`, `voxel_coords`, `parcel_coords`. Call `project_voxels`.
    *   **Assertions:** Check output list length and element dimensions (`V_v x k`). Check basic properties (finite values).

*   **TCK-NYSTROM-005: Back-projection Consistency (T-10)**
    *   **Synthetic Data:** Random `k`-space vector `w_k`.
    *   **Test:**
        1. Project to parcels: `w_parcel = w_k %*% t(hatsa_obj$v)` (using mean template `v`).
        2. Project `w_parcel` back to `k`: `w_k_recon = w_parcel %*% hatsa_obj$v` (assuming `v` is approx orthogonal basis).
        3. (Optional) Project `w_k` to voxels `C_vox = Phi_vox %*% w_k`, project `C_vox` back `w_k_recon_vox = pinv(Phi_vox) %*% C_vox`.
    *   **Assertions:** `expect_equal(cor(w_k, w_k_recon), 1, tolerance=1e-6)`. (Voxel round trip R² > 0.995). *Self-correction: Back-projection needs careful implementation.*

**Module: Advanced Features (Simulated/Conceptual)**

*   **TCK-ENSEMBLE-001: Ensemble Aggregation & Uncertainty (T-8)**
    *   **Test:** Simulate 32 noisy `R_i_e`. Implement/call geometric mean function. Calculate dispersion `σ_R`.
    *   **Assertions:** Check `R_i_agg` is SO(k). Check `σ_R` calculation. Check `R_i_agg` properties relative to `R_i_e` (e.g., centrality).

*   **TCK-WEIGHTED-001: Weighted Procrustes Effect (T-9)**
    *   **Test:** Implement/call `orthogonal_procrustes_weighted`. Create scenario with noisy subgroup anchors. Compare residuals with/without weighting.
    *   **Assertions:** Weighted error drops for subgroup, unchanged for others.

*   **TCK-PROVENANCE-001: Provenance Storage in `hatsa_projector`**
    *   **Test:** Run `run_hatsa_core`.
    *   **Assertions:** Check `hatsa_obj$parameters` contains all input args (`k`, anchors, `k_conn`, `n_refine`). Check if QC flags (`σ_R`, etc.) are stored if computed.

---

This detailed list provides concrete, testable units covering the core logic, S3 interface, numerical stability, edge cases, and key conceptual components discussed in the blueprint and subsequent refinements.