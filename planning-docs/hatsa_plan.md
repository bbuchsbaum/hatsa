## Integration Plan: `hatsa` into `multivarious`

The goal is to represent the output of `run_hatsa_core` as an object that fits within the `multivarious` S3 class system, specifically deriving from `multiblock_biprojector`.

**Target Class:** `hatsa_projector` inheriting from `multiblock_biprojector`.

**Overall Strategy:**

1.  **Refactor `run_hatsa_core`:** Modify it to return a new `hatsa_projector` object instead of a plain list.
2.  **Define `hatsa_projector` Class & Constructor:** Create the S3 class and a constructor that takes the results from the HATSA computation and maps them to the slots required by `multiblock_biprojector`, also storing HATSA-specific results.
3.  **Implement Core S3 Methods:** Implement `print`, `summary`, `coef` (loadings), `scores` (aligned sketches), and `block_indices` for the new class.
4.  **Implement Key Functional Methods:** Implement `predict` (to align new subjects) and potentially `project_block`.
5.  **Address Preprocessing:** Acknowledge that HATSA's complex preprocessing (data -> graph -> L -> U_orig -> Rotate) doesn't fit the standard `multivarious` `preproc` model easily. The `predict` method will encapsulate this logic for new data. The `preproc` slot in the object might initially be simple (e.g., `prep(pass())`).

**Mapping HATSA Outputs to `multiblock_biprojector` Slots:**

Let `hatsa_results` be the list currently returned by `run_hatsa_core` (v0.3.0).
Let `N` be the number of subjects, `V_p` the number of parcels, and `k` the spectral rank.

*   **`v` (Loadings/Coefficients - `V_p x k`):** Represents the group-level structure in the aligned space.
    *   **Chosen Representation:** The *mean aligned sketch*.
    *   **Calculation:** `v = Reduce('+', hatsa_results$U_aligned_list) / N`
*   **`s` (Scores - `(N * V_p) x k`):** Represents the individual parcel embeddings within the aligned space.
    *   **Chosen Representation:** Stacked aligned sketches.
    *   **Calculation:** `s = do.call(rbind, hatsa_results$U_aligned_list)`
*   **`sdev` (Component Standard Deviations - length `k`):** Represents the scale or variance associated with each component.
    *   **Chosen Representation:** Default to `rep(1, k)`. HATSA doesn't directly output component variances in the same way PCA does. Could be potentially derived later from eigenvalues of an average aligned Laplacian or variance of scores `s`.
*   **`preproc`:** The preprocessing pipeline.
    *   **Chosen Representation:** Initially `prep(pass())`. This signifies that the object stores *already processed* results (aligned sketches). A custom `pre_processor` could be developed later *if* needed, but the `predict` method is likely sufficient.
*   **`block_indices`:** A list defining which rows in `s` belong to which subject block.
    *   **Calculation:** `block_indices = split(seq_len(N * V_p), rep(seq_len(N), each = V_p))`
*   **Additional Slots (Stored directly in the `hatsa_projector` list):**
    *   `R_final_list`: From `hatsa_results`.
    *   `U_original_list`: From `hatsa_results`.
    *   `T_anchor_final`: The final group anchor template.
    *   `parameters`: Input parameters used (`k`, anchors, sparsification, etc.).
    *   `method`: `"hatsa_core"`

**The `predict` Method Challenge:**

A standard `projector`'s `project(obj, newdata)` method typically implies `apply_transform(obj$preproc, newdata) %*% coef(obj)`. HATSA alignment for new data is more complex:

1.  `newdata` (time-series) -> Build Graph -> Compute Laplacian -> Get `U_orig_new`.
2.  Extract anchors `A_orig_new`.
3.  Find rotation `R_new` by aligning `A_orig_new` to the *stored* `T_anchor_final` from the fitted object.
4.  Compute aligned sketch `U_aligned_new = U_orig_new %*% R_new`.

This full process needs to be implemented in `predict.hatsa_projector`. It doesn't fit the standard `preproc %*% v` pattern.

## Consolidated Integration Plan

We'll integrate HATSA into `multivarious` by enhancing the `hatsa_projector` object and adding new specific methods for voxel projection. The core `hatsa_projector` will handle parcel-level data as initially planned, while new functionality will support voxel-level analysis. This approach leverages the `multivarious` S3 system for clean dispatch.

multivariois is an R package that has a system for S3 classes and methods for single and multiblock data.

https://github.com/bbuchsbaum/multivarious

**Key Goals:**

1.  **Establish `hatsa_projector`:** Create an S3 class `hatsa_projector` inheriting from `multiblock_biprojector`. This object will be the primary output of `run_hatsa_core`.
2.  **Standard S3 Methods:** Implement essential S3 methods for `hatsa_projector` (e.g., `print`, `coef`, `scores`, `summary`).
3.  **Parcel-Level Prediction:** Implement a `predict` method for aligning new parcel-level subject time-series data.
4.  **Voxel-Level Projection:** Implement a `project_voxels` method for projecting voxel-level time series into the common HATSA space. This involves:
    *   Storing original subject-specific spectral components (eigenvectors and eigenvalues from the parcel-level graph Laplacian) within the `hatsa_projector` object.
    *   Creating an internal helper function (`compute_voxel_basis_nystrom`) to compute the Nyström voxel basis, which maps voxel data to the parcel-based spectral sketch.
    *   Defining an S3 method `project_voxels.hatsa_projector` to perform the projection using this basis and then rotate the results into the common aligned space.

The `hatsa_projector` will remain the central fitted object. Voxel projection is treated as a transformation method specific to this object (`project_voxels`), distinct from the standard `project` method (which we reserve for parcel-level data). This approach keeps the core projector API clean while providing the necessary specialized functionality for HATSA's voxel-level interface. The internal helper `compute_voxel_basis_nystrom` encapsulates the Nyström mathematics.

### Consolidated Ticket List

**Phase 1: Core `hatsa_projector` Class & Functionality**

1.  **Ticket 1 (Define `hatsa_projector` S3 Class & Constructor):**
    *   Define the S3 class `hatsa_projector` inheriting from `multiblock_biprojector`.
    *   Create the constructor function `hatsa_projector(hatsa_core_results)`:
        *   Perform calculations to map HATSA outputs to `v` (mean aligned sketch), `s` (stacked aligned sketches), `sdev` (default to `rep(1,k)`), and `block_indices`.
        *   Set `preproc = prep(pass())`.
        *   Store essential HATSA-specific results: `R_final_list` (rotations), `U_original_list` (original unaligned sketches), `T_anchor_final` (group anchor template), and input `parameters`.
        *   **Crucially, store `Lambda_original_list` (list of original *k* eigenvalues per subject from parcel-level decomposition) required for Nyström voxel projection.**
        *   Set the class attribute for S3 dispatch.

2.  **Ticket 2 (Refactor `run_hatsa_core`):**
    *   Modify `run_hatsa_core` to:
        *   **Ensure it computes and makes available the original subject eigenvalues (`Lambda_orig_parcel_i`) alongside the eigenvectors (`U_orig_parcel_i`). The function performing the spectral sketch (e.g., `compute_spectral_sketch_sparse`) must return both values and vectors.**
        *   Gather all necessary results (including `U_aligned_list`, `R_final_list`, `U_original_list`, `Lambda_original_list`, `T_anchor_final`, `parameters`).
        *   Call the `hatsa_projector` constructor with these comprehensive results.
        *   Return the fully populated `hatsa_projector` object.

3.  **Ticket 3 (Implement Basic S3 Methods):**
    *   `print.hatsa_projector(x, ...)`: Display N subjects, Vp parcels, k components, method.
    *   `coef.hatsa_projector(object, ...)`: Return `object$v`.
    *   `scores.hatsa_projector(x, ...)`: Return `x$s`.
    *   `sdev.hatsa_projector(x, ...)`: Return `x$sdev`.
    *   `block_indices.hatsa_projector(x, ...)`: Return `x$block_indices`.
    *   Verify that `ncomp.projector(x)` and `shape.projector(x)` (from `multivarious`) work correctly via inheritance. (Verified by structural review: `hatsa_projector` stores `k` in `parameters$k` and `ncol(coef(object))`, and provides standard dimensional components `coef`, `scores`, `block_indices` that `shape.projector` would likely use. Definitive test requires runtime.)

**Phase 2: HATSA Generics for Parcel-Level Data**

4.  **Ticket 4 (Implement `predict.hatsa_projector` for Parcel Data):**
    *   Signature: `predict.hatsa_projector(object, newdata_list, ...)`
    *   Input `object` is the fitted `hatsa_projector`.
    *   Input `newdata_list` is a list of time-series matrices (for parcels) for new subjects.
    *   Implement the full HATSA alignment process for each new subject's parcel data to compute `U_aligned_new`.
    *   Return a list of these `U_aligned_new` matrices.

5.  **Ticket 5 (Implement `project_block.multiblock_projector`):**
    *   Signature: `project_block.multiblock_projector(x, newdata_ignored = NULL, block, ...)` where `x` is a `hatsa_projector`.
    *   If `newdata_ignored` is NULL, extract the stored aligned sketch for `block` from `x$s`.
    *   If `newdata_ignored` is provided (as parcel time-series for one subject), call `predict.hatsa_projector` to get the aligned sketch.
    *   Return the `V_p x k` aligned sketch matrix for the specified block.

6.  **Ticket 6 (Implement `summary.hatsa_projector`):**
    *   Calculate and return basic diagnostics, e.g., mean anchor alignment error using stored `U_original_list`, `R_final_list`, and `T_anchor_final`.

**Phase 3: Voxel Projection Implementation**

7.  **Ticket 7 (Implement Nyström Voxel Basis Helper: `compute_voxel_basis_nystrom`):**
    *   Create a new internal helper function (e.g., in `R/utils.R` or `R/voxel_projection.R`).
    *   Signature: `compute_voxel_basis_nystrom(voxel_coords, parcel_coords, U_orig_parcel, Lambda_orig_parcel, n_nearest_parcels, kernel_sigma, ...)`
    *   Inputs: `voxel_coords` (V_v x 3), `parcel_coords` (V_p x 3), `U_orig_parcel` (V_p x k subject original eigenvectors), `Lambda_orig_parcel` (length k subject original eigenvalues, ensure > 0).
    *   Steps:
        *   Find `n_nearest_parcels` for each voxel.
        *   Compute sparse affinity matrix `W_vox_parc` (V_v x V_p) using a Gaussian kernel. Consider row-normalization.
        *   Handle near-zero eigenvalues in `Lambda_orig_parcel` when inverting.
        *   Compute voxel basis `Phi_voxel = W_vox_parc %*% U_orig_parcel %*% diag(1/Lambda_orig_parcel_safe)`.
    *   Output: Dense matrix `Phi_voxel` (V_v x k).

8.  **Ticket 8 (Implement `project_voxels.hatsa_projector` S3 Method):**
    *   Define a new S3 generic `project_voxels`.
    *   Define method `project_voxels.hatsa_projector`.
    *   Signature: `project_voxels(object, voxel_timeseries_list, voxel_coords, parcel_coords, n_nearest_parcels = 10, kernel_sigma = 5.0, ...)`
    *   Inputs: fitted `hatsa_projector` `object`, list of voxel time-series matrices `voxel_timeseries_list` (T_i x V_v), `voxel_coords`, `parcel_coords`.
    *   Process for each subject `i`:
        *   Retrieve `U_orig_parcel_i`, `Lambda_orig_parcel_i`, and `R_i` from `object`.
        *   Call `compute_voxel_basis_nystrom` to get `Phi_voxel_i` (V_v x k).
        *   Get voxel time series `X_voxel_ts_i` (T_i x V_v).
        *   Compute subject-specific voxel spectral coefficients: `C_voxel_coeffs_i = (1/T_i) * Matrix::t(X_voxel_ts_i) %*% Phi_voxel_i`. Result: (V_v x k).
        *   Rotate to common space: `C_voxel_aligned_i = C_voxel_coeffs_i %*% object$R_final_list[[i]]`. Result: (V_v x k).
    *   Output: A list of aligned voxel coefficient matrices (`C_voxel_aligned_i`).

**Phase 4: Testing**

9.  **Ticket 9 (Unit Tests for Core Functionality - Phases 1 & 2):**
    *   A new test file `tests/testthat/test-hatsa_core_functionality.R` has been created.
    *   The file includes `test_that` blocks outlining tests for:
        *   `hatsa_projector` constructor integrity and `run_hatsa_core` output (S3 class, key components like `Lambda_original_list`, parameters).
        *   Output dimensions and types of S3 methods: `coef`, `scores`, `sdev`, `block_indices`.
        *   Output dimensions and types of `predict.hatsa_projector`, including basic error handling.
        *   Output dimensions and types of `project_block.hatsa_projector` (for stored and new data), including error handling.
        *   Edge cases: `spectral_rank_k = 1` and small `N_subjects` (e.g., N=2).
        *   Basic execution of `summary.hatsa_projector`.
    *   Helper functions for mock data generation are included.
    *   Further refinement of mock data and specific assertions will be needed.
    *   The `tests/testthat.R` helper file has also been created.
    *   **Status:** Outline complete. Requires user implementation with robust mock data and specific assertions.

10. **Ticket 10 (Unit Tests for Voxel Projection - Phase 3):**
    *   A new test file `tests/testthat/test-voxel_projection.R` has been created.
    *   The file includes `test_that` blocks outlining tests for `compute_voxel_basis_nystrom` and `project_voxels.hatsa_projector`, including basic functionality, dimensions, parameter options, edge cases, and error handling.
    *   Helper functions for generating mock data are included.
    *   Further refinement and consistency checks (e.g., effect of rotation matrices) will be needed.
    *   **Status:** Outline complete. Requires user implementation with robust mock data and specific assertions.

**Phase 5: Documentation**

11. **Ticket 11 (Documentation for Core Functionality - Phases 1 & 2):**
    *   The `@return` section of `run_hatsa_core` documentation has been updated to describe the `hatsa_projector` object.
    *   Class-level documentation for `hatsa_projector` has been added.
    *   Existing S3 method documentation (`print`, `summary`, `coef`, `scores`, `sdev`, `block_indices`, `project_block`) has been reviewed and deemed adequate for now (further examples could be added later).
    *   A comprehensive example has been added to `run_hatsa_core` documentation demonstrating basic usage.
    *   **Status: Complete.**

12. **Ticket 12 (Documentation for Voxel Projection - Phase 3):**
    *   Documentation for the `project_voxels` S3 generic and the `project_voxels.hatsa_projector` method (including the detailed `@section Methodological Details` added in V-D1) has been reviewed.
    *   The `@field` descriptions in the `hatsa_projector` class documentation cover the necessary components for voxel projection.
    *   An example demonstrating the setup (including mock `run_hatsa_core`) and usage of `project_voxels.hatsa_projector` has been added to its Roxygen documentation.
    *   **Status: Complete.**

## Phase 6: Validation Metrics for Core HATSA (using Toy Data)

This phase outlines key metrics for evaluating the performance of the core HATSA algorithm, particularly when using synthetic data with a known ground truth, such as the output from `make_toy_hatsa` (described in the `core-hatsa-toy-example` vignette).

These metrics help quantify how well HATSA recovers the underlying structure and parameters of the data.

1.  **Rotation Recovery (SO(k) Alignment):**
    *   **Description:** Assesses how accurately the subject-specific rotation matrices (`R_i`) are estimated relative to the true rotations used to generate the toy data.
    *   **Function:** `hatsa::misalign_deg(R_est, R_true)`.
    *   **Metric:** The geodesic distance on SO(k), reported in degrees. Lower values indicate better recovery.
    *   **Relevance:** Directly evaluates the efficacy of the Generalized Procrustes Analysis (GPA) step within HATSA.

2.  **Group Template (`v`) Recovery:**
    *   **Description:** Evaluates how well the estimated group-level spectral template (`object$v` in `hatsa_projector`) matches the true ground-truth spectral basis (`U_true` from toy data).
    *   **Method:**
        1.  Align the estimated `object$v` to `toy_data$U_true` using Procrustes analysis (e.g., via `vegan::procrustes`). The `U_true` is typically considered the target.
    *   **Metrics:**
        *   **Correlation:** Pearson correlation between the vectorized aligned `v` and vectorized `U_true`.
        *   **Frobenius Norm of Difference:** `norm(aligned_v - U_true, "F")`.
    *   **Relevance:** Indicates how well HATSA identifies the common underlying spectral structure shared across subjects.

3.  **Anchor Template (`T_anchor_final`) Recovery:**
    *   **Description:** Assesses the fidelity of the estimated group anchor template (`object$T_anchor_final`) compared to the true spectral basis at the anchor locations.
    *   **Method:**
        1.  Extract the true spectral basis restricted to the anchor parcels: `U_true_anchors = toy_data$U_true[anchor_indices_toy, ]`.
        2.  Align `object$T_anchor_final` to `U_true_anchors` using Procrustes analysis. (Note: `U_true_anchors` might itself need to be aligned to the same consensus frame as `object$v` if the overall orientation of `U_true` is arbitrary relative to HATSA's output consensus).
    *   **Metrics:**
        *   **Correlation:** Pearson correlation between vectorized aligned `T_anchor_final` and `U_true_anchors`.
        *   **Frobenius Norm of Difference:** `norm(aligned_T_anchor_final - U_true_anchors, "F")`.
    *   **Relevance:** Checks if the learned anchor representation accurately reflects the true underlying components at those specific locations, crucial for the alignment process.

4.  **Eigenvalue Spectrum Fidelity (for Core HATSA):**
    *   **Description:** Compares the subject-specific eigenvalues obtained from HATSA (`object$Lambda_original_list` for core HATSA) with the true eigenvalues of the subject-specific graph Laplacians from which the toy data's signals were generated.
    *   **Method:**
        1.  The `make_toy_hatsa` generator uses a common `L_true_ring` for the signal component `U0` for all subjects. The `object$Lambda_original_list` from a `core_hatsa` run on `toy_data$X_list` reflects eigenvalues from Laplacians derived from these noisy, rotated time series.
        2.  A more direct comparison would involve generating true per-subject Laplacians `L_true_subject_i` (e.g., if noise or structure varied slightly per subject before rotation in a more complex toy model) and comparing their first `k` eigenvalues to `object$Lambda_original_list[[i]]`.
        3.  For the current `make_toy_hatsa`, `toy_data$eigenvalues_true` are from the common `L_true_ring`. We can compare the distribution/values in `object$Lambda_original_list` to these, acknowledging that the former come from empirical graphs.
    *   **Metrics:**
        *   **Correlation:** For each subject, correlate their estimated `k` eigenvalues with the first `k` true eigenvalues (if applicable).
        *   **Mean Squared Error (MSE):** For each subject, MSE between estimated and true eigenvalues.
    *   **Relevance:** Indicates how well the spectral characteristics (e.g., energy distribution across components) of individual subjects' graph structures are captured. This is more nuanced as `Lambda_original_list` reflects the outcome of HATSA's graph construction from noisy data.

5.  **Overall Stability and Robustness Assessment:**
    *   **Description:** Evaluates how the recovery metrics (1-4) change as key parameters of the data generation or HATSA algorithm are varied.
    *   **Methodology:**
        *   Systematically vary parameters in `make_toy_hatsa`: `snr`, `k` (spectral rank in generator vs. `spectral_rank_k` in HATSA), `Nsubj`, `Tlen`.
        *   Systematically vary HATSA parameters: choice/number of `anchor_indices`, graph construction parameters (`k_conn_pos`, `alpha_laplacian`), `n_refine`.
    *   **Evaluation:** Observe trends and sensitivity of `misalign_deg`, group template recovery, anchor recovery, and eigenvalue fidelity.
    *   **Relevance:** Understands the operational envelope of HATSA, its sensitivity to noise, parameter choices, and potential failure modes.

Implementation of these metrics can be done within the `core-hatsa-toy-example.Rmd` vignette for demonstration, or as separate utility functions within the package if they are intended for broader use or more systematic benchmarking scripts.

## Unit Test Implementation Plan (User Task for Tickets 9 & 10)

This section outlines the next steps for implementing the unit tests based on the outlines created in `test-hatsa_core_functionality.R` and `test-voxel_projection.R`.

**General Approach:**
*   Use `set.seed()` for reproducibility of any random mock data.
*   Employ small, fixed datasets where possible for predictable results, especially for complex calculations.
*   Consider creating mock `hatsa_projector` objects directly to isolate tests for functions like `predict`, `project_voxels`, etc., avoiding repeated `run_hatsa_core` calls.
*   Use `testthat::expect_equal(tolerance=...)` for floating-point comparisons.
*   Use `testthat::skip_if_not_installed()` for optional dependencies (`RSpectra`, `RANN`).
*   Add specific `expect_error()` and `expect_warning()` checks for invalid inputs or conditions.

**Specific Implementation Steps:**

1.  **Ticket 9 (`test-hatsa_core_functionality.R`):**
    *   **`run_hatsa_core`/Constructor:** Refine mock data (potentially fixed, smaller), add dimension checks for all relevant outputs (`v`, `s`, `R_final_list`, `U_original_list`, `T_anchor_final`), check eigenvalue positivity.
    *   **S3 Methods:** Verify `sdev` values, check `block_indices` partitioning.
    *   **`predict`/`project_block`:** Use mocked `hatsa_projector`, test specific inputs/outputs (e.g., identical inputs -> identical outputs for predict, stored sketch retrieval vs. projection for project_block).
    *   **Edge Cases:** Verify behavior and dimensions for `k=1` and `N=2`.
    *   **`summary`:** Check calculation/presence of `mean_anchor_alignment_error`.

2.  **Ticket 10 (`test-voxel_projection.R`):**
    *   **`compute_voxel_basis_nystrom`:** Use small fixed inputs for numerical verification if possible. Test `row_normalize_W` effect. Check behavior near `eigenvalue_floor`.
    *   **`project_voxels`:** Use mocked `hatsa_projector`. Implement consistency checks (e.g., identical inputs -> outputs related by rotation matrices `R_i`). Add more error condition tests.

**Phase Z: Post-Review Enhancements (Addressing Voxel Projection Feedback)

This phase addresses actionable feedback from a second review of the voxel projection implementation, focusing on enhancing flexibility and robustness.

23. **Ticket Z-R1 (Flexible Voxel-Parcel Affinities):**
    *   **Goal:** Allow users to bypass internal k-NN/Gaussian kernel computation by providing a pre-computed affinity matrix.
    *   **Action:** Modified `compute_voxel_basis_nystrom` to accept an optional `W_vox_parc` argument (defaulting to `NULL`). If provided, use it directly (after optional row normalization). If `NULL`, use current k-NN/Gaussian logic.
    *   **Docs:** Updated `compute_voxel_basis_nystrom` documentation.
    *   **Status: Complete.**

24. **Ticket Z-R2 (Support Different Voxel Data Types):**
    *   **Goal:** Correctly handle projection of data other than time-series (e.g., beta maps) by controlling scaling.
    *   **Action:** Added a `data_type = c("timeseries", "coefficients")` argument to `project_voxels.hatsa_projector`. Only apply `(1/T_i)` scaling if `data_type == "timeseries"`. Corrected projection calculation to yield `T_i x k` coefficients.
    *   **Docs:** Updated `project_voxels.hatsa_projector` documentation (including `@return`).
    *   **Status: Complete.**

25. **Ticket Z-R3 (Rotation Matrix Validation):**
    *   **Goal:** Add safety checks for the rotation matrices used in projection.
    *   **Action:** In `project_voxels.hatsa_projector`, added checks after retrieving `R_i` to verify approximate orthonormality (`crossprod(R_i) == I`) and determinant (`det(R_i) approx +1`). Issues warnings if checks fail.
    *   **Docs:** Mention checks in `project_voxels.hatsa_projector` documentation.
    *   **Status: Complete.**

**Note on Other Feedback:**
*   **Sparse `Phi_voxel`:** The suggestion to return `Phi_voxel` as sparse is deferred. Analysis (Ticket V-R6) indicates `Phi_voxel` is generally dense, so sparse representation is unlikely to save memory. Memory usage is noted in docs.
*   **Threading Optimization:** Added as a future consideration (see Ticket 16 below).

**Phase 6: Advanced Integration (Future Considerations)**

13. **Ticket 13 (Future: Implement `reconstruct` methods):**
    *   Define `reconstruct.hatsa_projector` if a meaningful reconstruction from aligned sketches is desired (e.g., `s %*% t(v)` to approximate aligned parcel representations).

14. **Ticket 14 (Future: Implement `project_vars`):**
    *   Consider `project_vars.hatsa_biprojector` for projecting supplementary *parcel-level* variables into the aligned `k`-dimensional space.

15. **Ticket 15 (Future: Develop HATSA-specific `pre_processor`):**
    *   If a more canonical `projector` workflow is desired where `project(hatsa_obj, raw_timeseries)` handles everything, a complex HATSA-specific `pre_processor` object could be developed. This is challenging due to inter-subject dependencies.

16. **Ticket 16 (Future: Threading Optimization for Projection):**
    *   Investigate and potentially implement explicit multi-threading for large matrix multiplications in `project_voxels.hatsa_projector`, e.g., using `RhpcBLASctl` or confirming optimal use of threaded BLAS via `Matrix` package operations.

This plan cleanly integrates the voxel projection as a post-hoc method associated with the fitted HATSA model, using the `multivarious` S3 system for dispatch.

## Phase X: Nyström Implementation Refinements (Post Code Review)

This phase addresses feedback from a detailed code review of the Nyström voxel projection implementation, focusing on numerical robustness, statistical consistency, and algorithmic clarity for `compute_voxel_basis_nystrom` and related components.

16. **Ticket V-R1 (Coordinate Validation for Voxel Projection):**
    *   Implement checks or assertions in the voxel projection pipeline (e.g., within `project_voxels.hatsa_projector` or as a pre-flight check) to validate the consistency of input coordinate systems for `voxel_coords` and `parcel_coords` (e.g., units, approximate bounding box agreement if possible).
    *   Document assumptions about these coordinate systems clearly.

17. **Ticket V-R2 (Kernel Sigma Auto-Tuning):**
    *   Modify `compute_voxel_basis_nystrom` to accept an option like `kernel_sigma = "auto"`.
    *   If "auto", estimate `kernel_sigma` based on the distribution of distances to nearest parcels (e.g., `median_distance / sqrt(2)`, where `median_distance` is the median of the k_nn=1 nearest neighbor distances from `RANN::nn2`).
    *   Ensure the chosen heuristic is documented.

18. **Ticket V-R3 (Robust `W_vox_parc` Construction):**
    *   In `compute_voxel_basis_nystrom`, modify the construction of the sparse `W_vox_parc` matrix to robustly handle cases where multiple parcels might have identical coordinates. This could involve using `Matrix::sparseMatrix(..., repr = "T")` followed by an aggregation step (e.g., `Matrix::tgsAgg` or equivalent) to sum similarities for duplicate `(voxel_idx, parcel_idx)` pairs.
    *   Define and implement a clear strategy for voxels that have zero similarity to all `n_nearest_parcels` (resulting in all-zero rows in `W_vox_parc`). Options include: keeping them as zero-rows (documenting implications for `Phi_voxel`), or optionally removing such voxels from the `Phi_voxel` computation and output (requiring careful handling of voxel indices). Document the chosen approach.

19. **Ticket V-R4 (Laplacian Consistency for Nyström Formulation):**
    *   **Correction:** Based on external documentation, HATSA defaults to an **alpha-lazy random-walk normalized Laplacian (`L_rw_lazy`)**, not the unnormalized Laplacian (`L = D - W`) as previously inferred from the code.
    *   The function `compute_graph_laplacian_sparse` has been **modified** to compute the symmetric form of `L_rw_lazy = I - alpha * D^{-1} * W` (using `alpha=0.93` default, using absolute weights for D, ensuring symmetry).
    *   To ensure mathematical consistency for the Nyström extension with `L_rw_lazy`:
        *   The Nyström approach should involve row-normalization of the voxel-parcel affinities (`W_vox_parc`).
        *   This corresponds to using `row_normalize_W = TRUE` in `compute_voxel_basis_nystrom`.
        *   The default for `row_normalize_W` has been reverted to `TRUE`.
        *   Documentation for `compute_voxel_basis_nystrom` and `project_voxels.hatsa_projector` (@section) has been updated to reflect the `L_rw_lazy` Laplacian and the rationale for `row_normalize_W = TRUE`.

20. **Ticket V-R5 (DC Component Handling in Spectral Sketches):**
    *   Review of `compute_spectral_sketch_sparse` confirms the following (updated for `L_rw_lazy`):
        *   The function calculates eigenvectors for the symmetric alpha-lazy random-walk normalized Laplacian (`L_rw_lazy`).
        *   Eigenvectors are sorted by their eigenvalues (smallest first).
        *   An `eigenvalue_tol` argument (default `1e-8`) is used to identify and **explicitly remove** eigenvectors associated with the smallest eigenvalues (corresponding to the graph's trivial/DC component(s)).
        *   The final `k` components selected are the first `k` from this filtered set of "non-trivial" eigenvectors.
    *   Therefore, the `spectral_rank_k` parameter in HATSA effectively requests `k` **non-DC components** derived from `L_rw_lazy`.
    *   The `U_orig_parcel` and `Lambda_orig_parcel` passed to `compute_voxel_basis_nystrom` are indeed the intended non-DC components.
    *   No code changes required for this ticket beyond the Laplacian change in `compute_graph_laplacian_sparse` (covered by V-R4).

21. **Ticket V-R6 (Future - Memory Efficiency for `Phi_voxel`):**
    *   **Investigation:** `Phi_voxel` is computed as `W_vox_parc %*% (U_orig_parcel %*% Lambda_inv_diag)`. Since `U_orig_parcel` is dense, the product `U_orig_parcel %*% Lambda_inv_diag` is dense. Multiplying the sparse `W_vox_parc` (V_v x V_p) by this dense (V_p x k) matrix generally results in a **dense `Phi_voxel` matrix (V_v x k)**. Converting it to a sparse matrix format would not typically save memory.

## Phase Y: Voxel Projection Documentation Enhancements (Post Code Review)

This phase focuses on updating user-facing documentation to reflect the methodological choices and assumptions within the voxel projection pipeline, ensuring reproducibility and clarity for users and reviewers.

22. **Ticket V-D1 (Detailed Methodological Documentation for Voxel Projection):**
    *   A dedicated `@section Methodological Details for Voxel Projection:` has been added to the roxygen documentation of `project_voxels.hatsa_projector`.
    *   This section consolidates and clearly specifies:
        *   **Coordinate Systems:** Assumptions (e.g., MNI, millimeters) and validation checks.
        *   **Distance Metric:** Use of Euclidean distance by `RANN::nn2`.
        *   **Voxel Data Pre-processing:** Recommendations and assumptions for input voxel BOLD fMRI time-series.
        *   **Kernel Sigma (`kernel_sigma`):** Explanation of the parameter, guidance on selection, and the "auto" heuristic.
        *   **Nyström Formulation:** Statement reflects the use of the **alpha-lazy random-walk normalized Laplacian** for parcel components and explains the consistency with `row_normalize_W=TRUE` (default) in `compute_voxel_basis_nystrom`.
        *   **DC Component:** Clarification that parcel-level components are non-DC, as the DC component is filtered out during the initial spectral decomposition of `L_rw_lazy`.
    *   This completes the primary goal of creating comprehensive methodological documentation for the voxel projection feature within the R script's documentation.