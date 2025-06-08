# HATSA Enhancement Tickets: Task HATSA Robustness & Riemannian Geometry Integration

This document outlines two sets of tickets:
1.  **Part A (THFIX Series):** Focuses on refactoring, bug fixing, and robustifying the existing `task_hatsa` implementation, particularly `R/task_hatsa_helpers.R`, based on a detailed technical audit. These are prerequisites for reliable advanced feature development.
2.  **Part B (RGEOM Series):** Focuses on integrating Riemannian geometry metrics and functionalities into HATSA, building upon a stabilized `task_hatsa` foundation.

---

## **Part A: Refactoring and Bugfix Tickets for `task_hatsa` (Based on Audit)**

**Overall Goal for Part A:** Stabilize and robustify the `task_hatsa` implementation by addressing scope bugs, improving parallel safety, clarifying API contracts, and enhancing error handling, as identified in the audit.

---

### **Ticket THFIX-001: Refactor `validate_and_initialize_args()` for Completeness and Efficiency**

-   [ ] **Description:** Address issues related to `...` evaluation and ensure all necessary parameters used downstream are explicitly included in the returned `args` list.
-   [ ] **Tasks:**
    1.  [X] Evaluate `dots <- list(...)` only once at the beginning of the function.
    2.  [X] Explicitly add *all* parameters that are read from `args` in downstream functions (`process_single_subject`, `compute_task_matrices`, `shape_basis`, `perform_anchor_augmentation`, `prepare_and_run_gpa`, `perform_patch_alignment`, `construct_output`) to the `arg_out` list. This includes `k_conn_pos`, `alpha_laplacian`, `n_refine`, `task_data_list`, etc.
    3.  [X] Ensure `task_data_list` is correctly stored in `arg_out` if it's passed via `...` or as a named argument.
    4.  [ ] Consider splitting `validate_and_initialize_args` into smaller `validate_*` helpers for different parameter groups (e.g., `validate_core_params`, `validate_task_params`, `validate_gpa_params`) if it remains >100 lines after explicit field additions.
    5.  [ ] (Optional) Assign a class like `"hatsa_args"` to the returned list for easier type checking if beneficial.
-   [ ] **Audit Ref:** 2.1, 4 (Style/elegance for splitting).

---

### **Ticket THFIX-002: Improve Parallel Processing Safety and Error Handling in `process_subjects()`**

-   [X] **Description:** Address memory issues with captured globals in `future_lapply` and improve error reporting from workers.
-   [X] **Tasks:**
    1.  [X] Modify the `single_subject_processor_for_future` worker function:
        *   [X] Instead of relying on `common_args` to capture large lists like `subject_data_list` and `task_data_list`, pass only the current `idx` and have the worker function retrieve `subject_data_list[[idx]]` and `task_data_list[[idx]]` from these lists, which should be passed as *named arguments* to `future_lapply` (e.g., `future_lapply(X = subject_indices, FUN = ..., s_data_list = subject_data_list, t_data_list = task_data_list, common_args = args_subset_without_large_lists)`). This leverages `future.apply`'s more efficient handling of such arguments compared to global variable capture if `future.globals = FALSE` is not perfectly effective or too broad.
        *   [X] Alternatively, investigate using `future.globals::globalsByName()` to precisely specify smaller globals. (Covered by the chosen explicit passing strategy)
    2.  [X] Wrap the body of `single_subject_processor_for_future` (i.e., the call to `process_single_subject`) in a `tryCatch` block.
    3.  [X] If an error occurs in the worker, return a structured error object, e.g., `list(ok = FALSE, error_message = e$message, subject_idx = idx)`.
    4.  [X] In `process_subjects()` after `future_lapply` returns:
        *   [X] Iterate through `all_subject_results`.
        *   [X] Check for `!res$ok`. If an error is found, either:
            *   [-] Fail fast: Stop execution and report the first error encountered. (Chose collect errors)
            *   [X] Collect all errors: Continue processing, store `NULL` or a placeholder for failed subjects in the component lists (`W_conn_list`, etc.), and report all errors at the end or return them as part of the QC output. (The latter is often preferred for batch processing).
-   [ ] **Audit Ref:** 1 (Parallelism), 2.4.

---

### **Ticket THFIX-003: Ensure Correct Variable Scoping and Argument Passing to Helpers**

-   [X] **Description:** Fix identified bugs where helper functions attempt to access variables from `args` that are not present or rely on incorrect lexical scoping.
-   [X] **Tasks:**
    1.  [X] **`compute_task_matrices()`:**
        *   [X] Modify function signature to accept `L_conn_i` as an explicit argument: `compute_task_matrices(subject_idx, task_data_i, args, W_conn_i, L_conn_i)`.
        *   [X] Update the call site in `process_single_subject()` accordingly.
    2.  [X] **`shape_basis()`:**
        *   [X] Ensure all `args$...` values it reads (e.g., `args$lambda_blend_value`, `args$k_gev_dims`) are guaranteed to be in `args` due to THFIX-001.
        *   [X] When `L_task_i` is `NULL` (e.g., due to `W_task_i_raw` computation failure or `task_method="core_hatsa"` initially), and a fallback sketch is computed (e.g., from `L_conn_i`), ensure the `sketch$vectors` and `sketch$values` are correctly assigned to `result$U_original` and `result$Lambda_original`. The audit noted these might not be copied if `sketch` itself is `NULL`. The `if (is.null(sketch)) return(result)` handles this, but if `sketch` is valid, ensure assignment.
        *   [X] **Critical for Blending:** As per previous audit, if `task_method == "lambda_blend"`, this function must receive `W_conn_i` and `W_task_i` (both already z-scored) to compute `W_hybrid_i = (1-λ)W_conn_i + λW_task_i`, then `L_hybrid_i = compute_graph_laplacian_sparse(W_hybrid_i)`, and finally `compute_spectral_sketch_sparse(L_hybrid_i)`. The current implementation blends *Laplacians*, which is incorrect. This is a major correction related to THFIX-001 and the previous code audit. (Verified current code blends W matrices correctly).
    3.  [X] **`prepare_and_run_gpa()`:**
        *   [X] Ensure all `args$...` values are present due to THFIX-001.
        *   [X] If `omega_mode == "adaptive"` and `reliability_scores_for_gpa` (prepared by `prepare_reliability_scores`) ends up being `NULL` (because no valid scores were provided), either:
            *   [-] Default `omega_mode` to `"fixed"` within `perform_gpa_refinement` for that call if `reliability_scores_list` input is `NULL`. (Handled by next point)
            *   [X] Or, `prepare_and_run_gpa` should detect this and explicitly call `perform_gpa_refinement` with `omega_mode = "fixed"`. The current warning is good, but the downstream call needs to be robust.
-   [ ] **Audit Ref:** 2.2, 2.3, 2.5.

---

### **Ticket THFIX-004: API Clarity and Namespace Management**

-   [ ] **Description:** Improve API usability and R CMD CHECK compliance.
-   [ ] **Tasks:**
    1.  [ ] **Verbosity:** Implement a centralized way to handle the `verbose` flag. For instance, `message_stage` could be modified to check a package option (e.g., `getOption("hatsa.verbose", TRUE)`) or an environment variable, rather than passing `verbose` down through many functions. Alternatively, ensure `args$verbose` is consistently used. The current `interactive_only` flag in `message_stage` is good.
    2.  [ ] **Namespace Calls:** Systematically review all calls to functions from other packages (e.g., `Matrix::`, `RSpectra::`, `multivarious::`, `PRIMME::`, `RANN::`) and ensure they are explicit.
    3.  [ ] **Imports:** Add necessary `@importFrom` tags in the package Roxygen header (`hatsa-package.R`) for all imported functions to satisfy R CMD CHECK and clarify dependencies (e.g., `@importFrom Matrix Matrix Diagonal t crossprod sparseMatrix drop0 forceSymmetric is`). `multivarious::scores` is already handled by `all_generic.R`, which is fine.
    4.  [X] **Internal Function Export:** Ensure helper functions called within `future_lapply` workers (like `process_single_subject` and its own callees like `compute_subject_connectivity_graph_sparse`, etc.) are either exported from the package (if generally useful) or properly available in the worker's environment (usually handled if they are in the package namespace and the package is loaded). Using `future.globals = FALSE` means they *must* be in the namespace.
-   [ ] **Audit Ref:** 3.

---

### **Ticket THFIX-005: Code Style and Abstraction**

-   [ ] **Description:** Improve code readability and maintainability.
-   [ ] **Tasks:**
    1.  [ ] **`tryCatch` Abstraction:** Create a `safe_wrapper(expr, error_message_fmt, default_return = NULL)` function as suggested in the audit to reduce boilerplate `tryCatch` code.
    2.  [X] **Public API Object:** Confirm that `run_task_hatsa` returns only the S3 `task_hatsa_projector` object and that gigantic nested lists are confined to internal function returns. (Current `construct_output` already does this by calling the constructor).
-   [ ] **Audit Ref:** 4.

---

### **Ticket THFIX-006: Constructor Availability for `task_hatsa_projector`**

-   [X] **Description:** The audit notes `task_hatsa_projector()` constructor is called in `construct_output` (in `task_hatsa_helpers.R`) but might not be defined/imported in that file if it resides elsewhere (e.g., `task_hatsa_projector.R`).
-   [X] **Task:** Ensure `task_hatsa_projector()` is properly defined and accessible (e.g., exported from the package and re-imported via `@importFrom hatsa task_hatsa_projector` if necessary, or simply available in the package namespace if all files are part of the same package build). This is usually handled by R's package loading mechanism if `task_hatsa_projector.R` defines and exports the function.
-   [X] **Audit Ref:** 1 (Return object construction).

---

### **Ticket THFIX-007: Comprehensive End-to-End Unit Test for `task_hatsa`**

-   [ ] **Description:** Create a minimal, reproducible end-to-end unit test for `run_task_hatsa`.
-   [ ] **Tasks:**
    1.  [ ] Use synthetic data for 2-3 subjects with a very small `V_p` (e.g., 5-10 parcels) and small `k` (e.g., 2-3).
    2.  [ ] Test with a simple `task_method` (e.g., "core_hatsa" first, then "lambda_blend" with mock task data).
    3.  [ ] Run the test under `future::plan(future::multisession)` (or `sequential` for basic debugging) to catch parallel processing issues.
    4.  [ ] Assert the output object class, basic dimension consistency, and absence of errors.
-   [ ] **Audit Ref:** 6 (Minimal reproducible fix snippet and end-to-end unit test).

---
---

## **Part B: Revised Riemannian Geometry Integration Tickets (RGEOM Series)**

These tickets now assume a more robust `task_hatsa` foundation after Part A tickets are addressed.

**General Note on Existing Tools:** The R package `RiemBase` (Kisung You, CRAN) offers a suite of functions for manifold statistics and computations, including routines relevant to SPD matrices and Grassmann manifolds. While `hatsa` is implementing specific geometric tools tailored to its needs, `RiemBase` serves as a notable existing resource and could be considered for comparative validation or for accessing more optimized/alternative algorithms if future requirements arise.

---

### **RGEOM-001: Implement Core Riemannian Metric Functions for SPD Matrices**
-   [X] **(Largely Complete)** (As detailed in previous discussions: `regularize_spd_matrix` (as `.regularize_spd`), `matrix_sqrt_spd`, `matrix_logm_spd`, `matrix_expm_spd`, `airm_distance`, `logeuclidean_distance` (as `riemannian_distance_spd`), `frechet_mean_spd`, `logmap_spd` (specific versions for LogEuclidean and AIRM: `logmap_spd_logeuclidean`, `logmap_spd_airm`), `expmap_spd` (specific versions: `expmap_spd_logeuclidean`, `expmap_spd_airm`) implemented in `R/riemannian_geometry.R`)

---

### **RGEOM-002: Implement Distance Functions for Orthonormal Bases**
-   [X] **(Largely Complete)** (As detailed: `grassmann_distance`, potentially `stiefel_distance` in `R/riemannian_geometry.R`)
    -   [X] `grassmann_distance` implemented in `R/riemannian_geometry.R`.
    -   [ ] `stiefel_distance`: Consider if a Procrustes-based distance (`sqrt(2k - 2 * sum(svd(V^T U)$d))`) or simpler chordal distance (`||U-V||_F`) is needed. Geodesic Stiefel distance is more complex.

---

### **RGEOM-003: S3 Methods for Extracting SPD Representations from Projector Objects (Revised for Robustness)**
-   [X] **(Largely Complete)** **Description:** Create `get_spd_representations.hatsa_projector`.
    -   [X] S3 generic `get_spd_representations` defined in `R/spd_representations.R`.
    -   [X] Method `get_spd_representations.hatsa_projector` implemented in `R/spd_representations.R` for `type = "cov_coeffs"` and `"fc_conn"`.
    -   [X] Method `get_spd_representations.task_hatsa_projector` implemented in `R/spd_representations.R` for `type = "fc_task"`, `"fc_hybrid"` (retrieval-based) and calls parent method for others.
    -   [X] Helpers `.get_subject_aligned_sketch` and `.compute_cov_spectral_coeffs` implemented.
    -   [ ] TODO: Consider re-computation logic for `fc_task`/`fc_hybrid` in `task_hatsa_projector` method if source matrices are not stored.
-   [ ] **Function Signature & Logic (Revised):**
    *   `get_spd_representations(object, type = c("cov_coeffs", "fc_conn", "fc_task", "fc_hybrid"), subject_idx = NULL, regularize_epsilon = 1e-6, subject_data_list_for_fc = NULL, k_conn_params_for_fc = list(pos=10,neg=10), lambda_blend_for_hybrid = NULL, ...)`
    *   `type = "cov_coeffs"`: Computes `cov(.get_subject_aligned_sketch(object, i))` where `.get_subject_aligned_sketch` robustly extracts the `V_p x k` matrix for subject `i` from `object$s` and `object$block_indices`, handling potential NAs in `object$s`.
    *   `type = "fc_conn"`: Requires `subject_data_list_for_fc`. Computes `W_conn_i = compute_subject_connectivity_graph_sparse(subject_data_list_for_fc[[i]], ...)`. *Crucially, regularize `W_conn_i` to SPD*. The `W` matrices are not inherently SPD. The proposal notes "If the W matrices can be negative, add δI with δ ≳ |λmin|+ε to make them SPD". This must be implemented carefully.
    *   `type = "fc_task"` / `"fc_hybrid"` (for `task_hatsa_projector`): Retrieve *stored* `W_task_list[[i]]` or `W_hybrid_list[[i]]` (requires THFIX tickets to ensure these are computed and potentially stored or re-computable). Then regularize to SPD.
-   [ ] **Return:** List of SPD matrices.

---

### **RGEOM-004: S3 Methods for Riemannian Distance Matrices and Dispersion**
-   [X] **(Largely Complete)** (As detailed, but relies on robust RGEOM-003).
    -   [X] S3 generic `riemannian_distance_matrix_spd` and method `riemannian_distance_matrix_spd.hatsa_projector` implemented in `R/riemannian_methods_hatsa.R`.
    -   [X] S3 generic `riemannian_dispersion_spd` and method `riemannian_dispersion_spd.hatsa_projector` implemented in `R/riemannian_methods_hatsa.R`.
    -   [ ] TODO: Implement/test `task_hatsa_projector` specific methods or ensure `NextMethod()` behavior is sufficient.

---

### **RGEOM-005: S3 Method for Tangent Space Embedding**
-   [X] **(Largely Complete)** **Description:** Implement an S3 method `get_tangent_space_coords` to project SPD representations to a common tangent space.
    -   [X] Define S3 generic `get_tangent_space_coords` in `R/riemannian_methods_hatsa.R`.
    -   [X] Implement `get_tangent_space_coords.hatsa_projector` in `R/riemannian_methods_hatsa.R`:
        -   [X] Retrieves SPD matrices using `get_spd_representations`.
        -   [X] Computes their Fréchet mean (`frechet_mean_spd`) as the tangent point.
        -   [X] Uses the appropriate `logmap_spd_*` function to project each SPD matrix to the tangent space at the mean.
        -   [X] Returns a list including the tangent vectors (list of symmetric matrices), the mean, and the metric used.
    -   [ ] Plan for vectorization of tangent vectors for use in standard multivariate tools (e.g., as a helper or option).
    -   [ ] Consider `task_hatsa_projector` specific method or `NextMethod()`.
-   [X] **Function Signature & Logic:** (As detailed in `Riemanian_plan.md`, Section 7, and implemented in `R/riemannian_methods_hatsa.R`).

---

### **RGEOM-006: Integration of Riemannian Metrics into Summary, QC, and Parameter Guidance (Revised for `k` selection)**
-   [X] **(Largely Complete)** **Description:** Enhance `summary` methods and create dedicated QC/guidance tools.
-   [ ] **Tasks:**
    1.  [X] **`summary` methods:** Optionally include mean/median `riemannian_dispersion(object, representation_type = "cov_coeffs")`.
        -   [X] `summary.hatsa_projector` modified in `R/hatsa_projector.R` to accept `compute_riemannian_dispersion`, `riemannian_dispersion_type`, and `riemannian_dispersion_options` arguments.
        -   [X] If `compute_riemannian_dispersion = TRUE`, calls `riemannian_dispersion_spd` and stores results.
        -   [X] `print.summary.hatsa_projector` modified in `R/hatsa_projector.R` to display these.
        -   [ ] TODO: Test this functionality once `task_hatsa_projector` inherits or has its own `summary` updated. (Currently relies on `hatsa_projector`'s summary via `NextMethod()` if not overridden, but explicit check needed for task-specific data pass-through).
    2.  [X] **`k`-Stability Plot Function:** Create `plot_k_stability_hatsa(projector_list_over_k, metrics_to_plot = c("Hk", "CV_eigen_cov"))`.
        -   [X] Implemented as `plot_k_stability_hatsa` in `R/hatsa_qc_plots.R`.
        -   [X] Takes a list of `hatsa_projector` objects (for different `k`).
        -   [X] Calculates `H(k)` using `riemannian_dispersion_spd(type = "cov_coeffs")`.
        -   [X] Calculates `CV(σ_eig(Cov_coeff_i(k)))` (mean CV of eigenvalues of `cov_coeffs` matrices).
        -   [X] Generates `ggplot2` plots of specified metrics vs. `k`.
        -   [ ] TODO: Develop test data (a list of `hatsa_projector` objects for different `k`) to thoroughly test this plotting function.
    3.  [X] **MDS for QC/Visualization (SPD matrices):** (As per M&B, Sec. 3.4)
        -   [X] Implemented `plot_mds_spd_subjects` function in `R/hatsa_qc_plots.R`.
        -   [X] Function takes a `projector_object` and options for distance calculation (`spd_representation_type`, `dist_mat_options`).
        -   [X] Calls `hatsa::riemannian_distance_matrix_spd`.
        -   [X] Performs MDS using `stats::cmdscale`.
        -   [X] Generates a 2D scatter plot of subjects using `ggplot2`.
        -   [X] Allows integration of `subject_info` for labels, color, and shape aesthetics.
        -   [X] Returns a list with the plot object and MDS results.
        -   [ ] TODO: Develop test data (a `hatsa_projector` object and optionally `subject_info`) to test this plotting function.

---

### **RGEOM-007: Experimental Geo-HATSA Implementation**
-   [X] **(Largely Complete)** Relies on RGEOM-001 (SPD/SO(k) metrics) and optimizers for SO(k).
    -   [X] Generalized `frechet_mean_so_k` moved from `summary.hatsa_projector` to `R/riemannian_geometry.R` and `summary.hatsa_projector` updated to use it.
    -   [X] New GPA function `perform_geometric_gpa_refinement` implemented in `R/procrustes_alignment.R`. This function uses `frechet_mean_so_k` to find `R_bar` (mean rotation) and updates the consensus anchor template `T_anchor_geo` by averaging subject anchors after aligning them to `R_bar`'s orientation.
    -   [X] New core algorithm function `run_geo_hatsa_core` implemented in `R/geo_hatsa_core_algorithm.R`. This function is analogous to `run_hatsa_core` but calls `perform_geometric_gpa_refinement`.
    -   [X] `run_geo_hatsa_core` returns `R_bar_final` in its results list, in addition to standard HATSA outputs.
    -   [ ] TODO: Consider if/how `R_bar_final` and the "geo_hatsa_core" method type should be integrated into the `hatsa_projector` object (e.g., new class `geo_hatsa_projector` or extend existing one). For now, it's an extra item in the list returned by `run_geo_hatsa_core`.
    -   [ ] TODO: Extensive testing and validation of Geo-HATSA against standard HATSA would be needed to assess its properties and potential benefits.

---

### **RGEOM-008: Implement MRA-Select (Manifold-Regularized Anchor Selection - Revised for Inputs)**
-   [X] **(Largely Complete)** **Description:** Create `select_anchors_mra(...)`.
    -   [X] Function `select_anchors_mra` created in `R/anchor_selection_mra.R`.
-   [X] **Inputs:** 
    -   [X] `U_original_list_pilot`, `k_spectral_rank`, `m_target`, `total_parcels`.
    -   [X] `max_kappa`, `weight_inv_kappa`, `weight_dispersion`.
    -   [X] Optional: `initial_selection`, `candidate_pool`, `riemannian_dispersion_options`, `min_anchors_for_metrics`.
    -   [ ] `parcel_quality_info` is an argument but currently unused in core logic (reserved).
-   [X] **Logic:**
    1.  [X] Greedy forward selection implemented.
    2.  [X] At each step, evaluates candidate anchors by impact on:
        *   [X] `κ(EuclideanMean_pilots({U_pilot[A_sel_tentative, :]}))`: Condition number of the Euclidean mean of extracted anchor rows from pilot subjects' original sketches. Helper `.calculate_kappa` implemented.
        *   [X] `Dispersion_MB_cov_coeffs(A_sel_tentative)`: For each pilot, computes `S_i = cov(U_pilot[A_sel_tentative,:])`. Then calculates mean squared Riemannian distance of these `S_i` to their Fréchet mean. This is done via direct calls to `hatsa::frechet_mean_spd` and `hatsa::riemannian_distance_spd` within a helper (`.evaluate_anchor_set`).
            -   [ ] TODO: Refine dispersion calculation to more directly use `riemannian_dispersion_spd` if that function is adapted to handle lists of SPD matrices or if a dummy object construction is preferred for MRA-Select.
    3.  [X] Selects anchor maximizing `weight_inv_kappa * (1/κ) - weight_dispersion * dispersion_val`, subject to `kappa <= max_kappa`.
-   [X] **Output:** Sorted vector of selected anchor indices.
-   [ ] TODO: Develop comprehensive test cases with synthetic pilot data to validate selection behavior under different weighting schemes and input conditions.

---

### **RGEOM-009: Implement Bures-Wasserstein Barycenter for SPD Matrices and Integrate into HATSA**
-   [ ] **Description:** Implement a function to compute the Bures-Wasserstein (BW) barycenter of a list of SPD matrices. This barycenter is the geodesic mean under the BW metric and offers a robust way to average SPD matrices like graph Laplacians. Explore its application in HATSA for:
    *   Seeding the initial group spectral template: Using the eigenvectors of the BW-mean of subject Laplacians (after parcel registration) to form an initial `T_anchor` basis or to initialize the overall HATSA fit more broadly, replacing or augmenting current arithmetic/convex blend approaches.
    *   Providing a robust global prior for task-patch GEV (using mean Laplacian as the `A` matrix in GEV).
    *   Incrementally updating a global template for drift correction in a "data warehouse" scenario.
-   [ ] **Tasks:**
    1.  [ ] **Core BW Barycenter Function:**
        *   [ ] Implement `bures_wasserstein_barycenter(S_list, weights = NULL, initial_mean = NULL, max_iter = 50, tol = 1e-7, regularize_epsilon = 1e-6, verbose = FALSE)`.
        *   [ ] Base implementation on iterative algorithms (e.g., fixed-point iteration involving matrix square roots and `solve()`).
        *   [ ] Ensure input matrices `S_list` are regularized to SPD.
        *   [ ] Support optional `weights` for a weighted barycenter.
        *   [ ] Allow an `initial_mean` estimate.
    2.  [ ] **Incremental BW Barycenter Update:**
        *   [ ] Implement `bures_wasserstein_barycenter_update(current_barycenter, S_new, n_old_samples, weight_new = 1, ...)` based on a streaming update rule (e.g., Haasler & Frossard, Alg. 1).
    3.  [ ] **Laplacian Type Investigation:**
        *   [ ] Investigate and document suitable Laplacian types for BW averaging in HATSA (e.g., combinatorial `D-W`, normalized `I - D^-1/2 W D^-1/2`, random-walk `I - D^-1 W`).
        *   [ ] Ensure chosen Laplacians are made SPD (e.g., via `+ J/N` or other regularization) before BW computation, if necessary.
    4.  [ ] **Integration Point A (Initial Group Template):**
        *   [ ] Modify `run_hatsa_core` (or create `run_bw_hatsa_core` or add an option) to optionally compute subject Laplacians, find their BW barycenter, and use the eigenvectors of this barycenter as the initial `T_anchor` basis for GPA.
    5.  [ ] **Integration Point B (GEV Prior):**
        *   [ ] Explore using the BW-mean of subject Laplacians (or relevant task-derived Laplacians) as the `A` matrix (global prior) when solving GEV for task patches (e.g., in `shape_basis` or related logic in `task_hatsa_helpers.R`).
    6.  [ ] **Integration Point C (Drift Correction / Warehouse Sync):**
        *   [ ] Design a workflow or standalone function that uses `bures_wasserstein_barycenter_update` to refresh a global template (e.g., mean Laplacian or its eigenvectors) as new subject data becomes available.
    7.  [ ] **Unit Testing:**
        *   [ ] Create unit tests for `bures_wasserstein_barycenter` with known small examples.
        *   [ ] Create unit tests for `bures_wasserstein_barycenter_update`.
    8.  [ ] **Benchmarking:**
        *   [ ] Benchmark performance for typical HATSA parcel matrix sizes (e.g., `V_p` ~ 400).
-   [ ] **References/Notes:**
    *   Considers the natural geometry of SPD matrices.
    *   Node ordering of input matrices (e.g., Laplacians) must be homologous.
    *   Iterative steps typically involve `O(V_p^3)` operations (matrix square roots, inversions).
    *   Eigenvectors of the BW-mean Laplacian (obtained via `RSpectra::eigs_sym` or similar) can serve as a basis.
    *   The BW mean can be robust and less prone to being skewed by outliers compared to arithmetic mean.
    *   Potentially better conditioning and stability for downstream tasks.
    *   Consider impact on capturing local activations versus spatial bandwidth of eigen-kernels.

---

### **Future RGEOM Considerations / Potential Tickets**

Beyond the currently defined RGEOM tickets, the following areas represent potential future enhancements and research directions for the Riemannian geometry aspects of the HATSA package, drawing from the Geo-HATSA audit and ongoing development:

*   **RGEOM-FUTURE-001: Covariance Normalization in Geometric GPA**
    *   **Description:** Investigate and optionally implement row-wise normalization (e.g., `sqrt(degree)` weighting, `A_i_norm = D_i^(1/2) A_i`) for anchor matrices `A_i` prior to the Procrustes steps within `perform_geometric_gpa_refinement`.
    *   **Rationale:** As suggested in Geo-HATSA audit (Gap 2e), this could improve stability and performance if anchor contributions are heterogeneous.
    *   **Tasks:**
        1.  Research appropriate normalization schemes (e.g., based on row norms of `A_i`).
        2.  Implement as an optional step in `perform_geometric_gpa_refinement`.
        3.  Benchmark its impact on alignment quality and convergence.

*   **RGEOM-FUTURE-002: Full Riemannian Optimization for SO(k) GPA**
    *   **Description:** Explore and potentially implement a full Riemannian optimization approach (e.g., Riemannian gradient descent, trust-region methods) for the GPA step on SO(k) within `perform_geometric_gpa_refinement`.
    *   **Rationale:** To directly minimize a geodesic distance-based cost function on the SO(k) manifold, offering a more geometrically precise (though potentially more computationally intensive) alternative to the current SVD-based Euclidean approximation.
    *   **Tasks:**
        1.  Research suitable algorithms and their implementations (e.g., leveraging `geomstats` concepts or `manopt`).
        2.  Implement an alternative update rule for rotations `R_i` in `perform_geometric_gpa_refinement`.
        3.  Compare performance, convergence, and solution quality against the SVD-based approach.

*   **RGEOM-FUTURE-003: Comprehensive Benchmarking of Geo-HATSA**
    *   **Description:** Conduct thorough benchmarking of `run_geo_hatsa_core` (and the full Geo-HATSA pipeline if it evolves) against the standard Euclidean `run_hatsa_core`.
    *   **Rationale:** To systematically evaluate the benefits, trade-offs, accuracy, robustness, and computational cost of the geometric approach.
    *   **Tasks:**
        1.  Develop synthetic datasets with known ground truth for rotations and templates.
        2.  Apply to standard real-world neuroimaging datasets (e.g., HCP-15).
        3.  Define key metrics for comparison (e.g., ICC, template recovery error, rotation recovery error, computational time).
        4.  Publish or document findings.

*   **RGEOM-FUTURE-004: Enhanced Vectorization and Utility for Tangent Space Coordinates**
    *   **Description:** Improve the utility of tangent space embeddings obtained from `get_tangent_space_coords`.
    *   **Rationale:** Facilitate easier integration with standard R multivariate analysis tools.
    *   **Tasks:**
        1.  Implement a robust helper function (e.g., `.vectorize_tangent_vectors`) or add an option to `get_tangent_space_coords` to return a subject-by-feature data matrix from the list of tangent vector matrices.
        2.  Provide examples or vignette sections on using these vectorized coordinates for downstream analyses like PCA, clustering, or regression.
