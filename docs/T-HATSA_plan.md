The chosen name based on the discussion will be **task_hatsa (Task-Informed HATSA)** for its clarity and broad applicability, aligning with the `hatsa` package context.

---

## Proposal: task_hatsa (Task-Informed or "super-charged" HATSA) Implementation

**Objective:** Implement the task_hatsa variant within the R `hatsa` package. task_hatsa extends Core HATSA by incorporating task-derived information (e.g., from stimulus labels or experimental conditions) to produce a task-aware alignment, improving comparability and utility for downstream analyses involving evoked activity or stimulus encoding.

*Downstream "Rosetta Stone" Effect:* Every basis element admitted through validation—whether from λ-blend or the labGEV patch—can become a public axis in a shared template (like MA-HATSA). This allows future datasets to project onto these axes for directly comparable "task coordinates", effectively growing the shared "brain-function dictionary".

**Target Audience:** R Package Engineer / Developer

**Core Concept:** task_hatsa modifies Core HATSA primarily via two synergistic stages:
1.  **Stage A: Basis Shaping (Defining `U_original_list[[i]]`)**: Creates a `spectral_rank_k`-dimensional spectral basis reflecting intrinsic connectivity (`W_conn`) and potentially task-related similarity (`W_task`). Methods include:
    *   *(`task_method = "lambda_blend"`)*: **Linear `λ`-Blend:** Combine `W_conn` and `W_task` (potentially residualized `W_task_res`) into `W_hybrid` before computing the Laplacian and its eigenvectors.
    *   *(`task_method = "gev_patch"`)*: **GEV Patch:** Compute `U_core` from `L_conn`. Separately, use Generalized Eigenvalue Decomposition (`L_task v = λ L_conn v`) to find task-specific eigenvectors `U_taskGEV` orthogonal to connectivity, stored as a separate patch.
    *   *(`task_method = "core_hatsa"`)*: Uses only `W_conn`, replicating Core HATSA basis shaping.
2.  **Stage B: Alignment Fine-Tuning (Using `A_aug_i` in Procrustes)**: Uses explicit task feature representations to guide the subject-specific rotation `R_i`. Method:
    *   *(`row_augmentation = TRUE`)*: **Row Augmentation:** Add condition-specific representations (e.g., projected activation vectors) as extra rows to the anchor matrix (`A_aug_i`) used in the Procrustes alignment, often with differential weighting (`omega_weights`).

**Algorithmic Overview (task_hatsa Workflow):**

**(Assumes Core HATSA functions like `compute_subject_connectivity_graph_sparse`, `compute_graph_laplacian_sparse`, `compute_spectral_sketch_sparse`, `solve_procrustes_rotation`, `perform_gpa_refinement` from the `hatsa` package are available and robust).**

1.  **Inputs for `run_task_hatsa`:**
    *   `subject_data_list`: List of `T x V_p` matrices (parcel time series). (Same as `run_hatsa_core`)
    *   `anchor_indices`: Vector of `m` canonical anchor indices. (Same as `run_hatsa_core`)
    *   `spectral_rank_k`: Integer, dimensionality of the primary spectral sketch. (Primary dimension for output `task_hatsa_projector`).
    *   `task_data_list`: List of pre-processed task-related data per subject (e.g., a `C x V_p` matrix of activations/betas per subject). Required if `task_method != "core_hatsa"`. Data is assumed cross-validated *before* input if needed. If `omega_mode == "adaptive"`, this list may need to include or allow computation of reliability scores (e.g., `R^2_p`) for each beta map.
    *   `task_method`: Character string: `"core_hatsa"`, `"lambda_blend"`, `"gev_patch"`. Default: `"lambda_blend"`.
    *   `lambda_blend_value`: Numeric `λ ∈ [0,1]`. Used if `task_method == "lambda_blend"`. Default 0.15.
    *   `k_gev_dims`: Integer, requested dimension for GEV patches. Used if `task_method == "gev_patch"`. Default 10-15.
    *   `row_augmentation`: Logical, whether to add projected task features to anchors. Default `TRUE` if `task_data_list` contains suitable data. Can be controlled by `--no_row_aug` equivalent.
    *   `residualize_condition_anchors`: Logical, whether to residualize projected condition anchors (`Z_i`) against parcel anchors (`A_parc_i`). Default `FALSE`. Experimental feature.
    *   `omega_weights`: List specifying weights for Procrustes alignment if `omega_mode == "fixed"`. E.g., `list(parcel = 1.0, condition = 0.5)`. If `NULL` and `omega_mode == "fixed"`, defaults to `list(parcel = 1.0, condition = 0.5)`. The diagonal weighting matrix `Ω` derived from these will be trace-rescaled.
    *   `omega_mode`: Character string: `"fixed"` or `"adaptive"`. Default `"fixed"`. If `"adaptive"`, weights for condition rows are adapted based on reliability (e.g., `R^2_p`) of input beta maps.
    *   `alpha_laplacian`: Numeric, laziness parameter for `compute_graph_laplacian_sparse`. Default 0.93.
    *   `k_conn_pos`, `k_conn_neg`: Sparsification for `W_conn`. (Same as `run_hatsa_core`).
    *   `k_conn_task_pos`, `k_conn_task_neg`: Sparsification for `W_task` (used by `compute_W_task_*` helpers).
    *   `similarity_method_task`: Method ('pearson', 'spearman', function) for `compute_W_task_*` helpers. Default 'pearson'.
    *   `n_refine`: Number of GPA iterations. (Same as `run_hatsa_core`).
    *   `check_redundancy`: Logical, whether to perform the `W_conn`/`W_task` redundancy check. Default `TRUE` if `task_method` is not `core_hatsa`.
    *   `redundancy_threshold`: Numeric, Spearman `rho` threshold for triggering residualization. Default **0.45**.
    *   `residualize_k_conn_proj`: Integer, number of `L_conn` eigenvectors to project `W_task` out of during residualization. Default **64**.
    *   `residualize_k_conn_labels`: Integer, k-NN value for re-sparsifying `W_task_res` after residualization.
    *   `gev_lambda_max`: Numeric, maximum GEV eigenvalue `λ` to retain for patches. Default **0.8**.
    *   `gev_stability_r_min`: Numeric, minimum split-half stability `r` for GEV eigenvectors. Default **0.6**.
    *   `gev_epsilon_reg`: Small numeric for regularizing `L_conn` in GEV. Default 1e-6.

2.  **Preprocessing & Redundancy Handling per Subject `i`:**
    *   `X_res_i = subject_data_list[[i]]`.
    *   Compute `W_conn_i` using `compute_subject_connectivity_graph_sparse(X_res_i, ..., k_conn_pos, k_conn_neg)`. (Output is sparse, z-scored).
    *   If `task_method != "core_hatsa"`: Compute initial `W_task_i_raw` using appropriate helper (e.g., `compute_W_task_from_activations(task_data_list[[i]], ..., k_conn_task_pos = k_conn_task_pos, k_conn_task_neg = k_conn_task_neg, similarity_method = similarity_method_task)`). (Output is sparse, z-scored).
    *   **Redundancy Check (if `check_redundancy == TRUE` and `W_task_i_raw` exists):**
        *   Compute `rho_i = compute_graph_correlation(W_conn_i, W_task_i_raw)` (Spearman, union+zero-fill).
        *   Store `rho_i`. Set `residualize_flag_i = FALSE`.
        *   If `rho_i >= redundancy_threshold`:
            *   Compute `L_conn_i = compute_graph_laplacian_sparse(W_conn_i, alpha = alpha_laplacian)`.
            *   **Residualize `W_task` (Plan B):** `W_task_i = residualize_graph_on_subspace(W_graph_to_residualize = W_task_i_raw, L_graph_for_projection = L_conn_i, k_eigenvectors_to_remove = residualize_k_conn_proj, k_nn_resparsify = residualize_k_conn_labels)`. (Output is sparse, re-z-scored).
            *   Set `residualize_flag_i = TRUE`.
        *   Else: `W_task_i = W_task_i_raw`.
    *   Else (`task_method == "core_hatsa"` or `check_redundancy == FALSE`): `W_task_i = W_task_i_raw` (if exists), `rho_i = NA`, `residualize_flag_i = NA`.

3.  **Basis Shaping per Subject `i` (Based on `task_method`):**
    *   **If `task_method == "core_hatsa"`:**
        *   `L_i = compute_graph_laplacian_sparse(W_conn_i, alpha = alpha_laplacian)`.
        *   `sketch = compute_spectral_sketch_sparse(L_i, spectral_rank_k, eigenvalue_tol = 1e-8)`.
        *   `U_original_list[[i]] = sketch$vectors`. `Lambda_original_list[[i]] = sketch$values`. `U_patch_list[[i]] = NULL`.
    *   **If `task_method == "lambda_blend"`:**
        *   `L_conn_i = compute_graph_laplacian_sparse(W_conn_i, alpha = alpha_laplacian)`.
        *   `L_task_i = compute_graph_laplacian_sparse(W_task_i, alpha = alpha_laplacian)` (using potentially residualized `W_task_i`).
        *   `L_hybrid_i = (1-lambda_blend_value)*L_conn_i + lambda_blend_value*L_task_i`.
        *   Check `eigengap_ratio_k` of `L_hybrid_i` (diagnostic, cf. threshold **~1.3-1.4**, possibly SNR dependent).
        *   `sketch = compute_spectral_sketch_sparse(L_hybrid_i, spectral_rank_k, eigenvalue_tol = 1e-8)`.
        *   `U_original_list[[i]] = sketch$vectors`. `Lambda_original_list[[i]] = sketch$values`. `U_patch_list[[i]] = NULL`.
    *   **If `task_method == "gev_patch"`:**
        *   `L_conn_i = compute_graph_laplacian_sparse(W_conn_i, alpha = alpha_laplacian)`.
        *   `L_task_i = compute_graph_laplacian_sparse(W_task_i, alpha = alpha_laplacian)` (using potentially residualized `W_task_i`).
        *   `gev_results = solve_gev_laplacian_primme(L_task_i, L_conn_i, k_request = k_gev_dims*2, lambda_max_thresh = gev_lambda_max, stability_r_thresh = gev_stability_r_min, epsilon_reg_L_conn = gev_epsilon_reg)`. Filter eigenvectors based on `lambda` and stability `r`. Store `U_patch_taskGEV_i` (actual `k_gev_dims` passing filters) and `Lambda_GEV_i`.
        *   Compute core sketch: `core_sketch = compute_spectral_sketch_sparse(L_conn_i, spectral_rank_k, eigenvalue_tol = 1e-8)`.
        *   `U_original_list[[i]] = core_sketch$vectors`. `Lambda_original_list[[i]] = core_sketch$values`.
        *   Store `U_patch_taskGEV_i`, `Lambda_GEV_i`, and GEV diagnostics separately.

4.  **Anchor Augmentation (If `row_augmentation == TRUE`):**
    *   For each subject `i`:
        *   `U_basis_i = U_original_list[[i]]` (this is the `U_i` from Appendix §, chosen by Basis Shaping).
        *   Extract parcel anchor rows `A_parc_i = U_basis_i[anchor_indices, , drop=FALSE]` (`m × k` matrix).
        *   If `task_data_list[[i]]` provides suitable features (e.g., a `C x V_p` activation matrix `Act_i`):
            *   Project features: `Z_i_projected = project_features_to_spectral_space(feature_matrix = t(Act_i), U_basis = U_basis_i)`. Note: `project_features_to_spectral_space` returns features in `k_dims_basis x C` format if input `feature_matrix` was `V_p x C`. We need `Z_i` as `C x k` for `rbind`. So, `Z_i = t(Z_i_projected)`.
            *   **Optional Residualization of Condition Anchors:** If `residualize_condition_anchors == TRUE`:
                *   `Z_i = residualize_matrix_on_subspace(matrix_to_residualize = Z_i, subspace_basis_matrix = A_parc_i)`. (This uses the new helper TCK-AUGALI-001.5).
            *   `A_aug_orig_i = rbind(A_parc_i, Z_i)`. Store associated condition labels.
        *   Else: `A_aug_orig_i = A_parc_i`.
    *   The list of `A_aug_orig_i` (or `A_parc_i`) matrices becomes `A_originals_list_for_gpa`.

5.  **Iterative Refinement (GPA):**
    *   Initialize `T_aug_anchor` from mean of `A_originals_list_for_gpa`.
    *   Call modified `perform_gpa_refinement(A_originals_list_for_gpa, n_refine, k=spectral_rank_k, omega_weights = omega_weights, omega_mode = omega_mode, reliability_scores_for_adaptive_omega = ...)` The modified function uses `solve_procrustes_rotation_weighted` internally.
    *   Final core rotation `R_final_list[[i]]`. `T_anchor_final` is the final `T_aug_anchor`.

6.  **Patch Alignment (If `task_method == "gev_patch"` and patches `U_patch_taskGEV_i` exist):**
    *   Align `U_patch_taskGEV_i` across subjects using their anchor rows and `solve_procrustes_rotation` to get `R_patch_lists[[taskGEV]][[i]]`. (See **Appendix G.4** for details on patch alignment architecture).

7.  **Construct Output (`task_hatsa_projector` object):**
    *   Create `task_hatsa_projector` object, inheriting from `hatsa_projector`.
    *   Populate standard slots using core components (`U_original_list`, `Lambda_original_list`, `R_final_list`, `T_anchor_final`, computed `U_aligned_list`, `v`, `s`, `sdev`, `block_indices`).
    *   Store parameters: include `method = "task_hatsa"`, and all task_hatsa specific inputs (`task_method`, `lambda_blend_value`, `k_gev_dims`, `alpha_laplacian`, `similarity_method_task`, `check_redundancy`, `redundancy_threshold`, `residualize_k_conn_proj`, `residualize_k_conn_labels`, `gev_lambda_max`, `gev_stability_r_min`).
    *   Add new slots for:
        *   `qc_metrics`: list per subject containing `rho_redundancy`, `was_residualized`.
        *   If GEV patch: `gev_patch_data` (list containing `U_patch_lists`, `R_patch_lists`, `Lambda_GEV_lists`, `gev_diagnostics`).
        *   If anchor augmentation: `anchor_augmentation_info` (list containing e.g., `condition_labels_for_anchors`, `T_aug_anchor_final`, `omega_weights_used` or `omega_mode` and related parameters, `condition_anchors_residualized_flag`).

---

**Granular Tickets for R Implementation (task_hatsa Focus)**

**(Assumes Core HATSA functions from `hatsa` package are mostly implemented and tested per previous tickets. New functions should follow existing naming conventions, e.g., `verb_noun_details`.)**

**Module: Graph Construction & Blending (Task-Specific)**
*   `[X]` **TCK-TSKGR-001:** Implement `compute_W_task_from_activations` function.
    *   Inputs: `activation_matrix` (a `C x V_p` matrix of pre-processed activation coefficients for subject `i`), `parcel_names` (`V_p` length), `k_conn_task_pos` (integer), `k_conn_task_neg` (integer), `similarity_method` (character string: "pearson" (default), "spearman", or a function that takes `activation_matrix` and returns a `V_p x V_p` similarity matrix).
    *   Internals: Compute `V_p x V_p` similarity matrix from `activation_matrix` using the specified `similarity_method`. Sparsify this matrix using `k_conn_task_pos` and `k_conn_task_neg`. Symmetrize and z-score the non-zero edge weights.
    *   Output: Sparse `W_task_i` (`dgCMatrix`, `V_p x V_p`).
*   `[X]` **TCK-TSKGR-002:** Implement `compute_W_task_from_encoding` function.
    *   Inputs: `encoding_weights_matrix` (a `V_p x N_features` matrix of pre-trained encoding model weights for subject `i`), `parcel_names` (`V_p` length), `k_conn_task_pos` (integer), `k_conn_task_neg` (integer), `similarity_method` (character string: "pearson" (default), "spearman", or a function that takes `encoding_weights_matrix` as input and returns a `V_p x V_p` similarity matrix).
    *   Internals: Compute `V_p x V_p` similarity matrix from `encoding_weights_matrix` using the specified `similarity_method`. For "pearson" or "spearman", this involves correlating `t(encoding_weights_matrix)` to compare parcel profiles. Sparsify this matrix using `k_conn_task_pos` and `k_conn_task_neg`. Symmetrize and z-score the non-zero edge weights.
    *   Output: Sparse `W_task_i` (`dgCMatrix`, `V_p x V_p`).
*   `[X]` **TCK-TSKGR-003:** Implement `compute_graph_correlation` function.
    *   Inputs: `W_graph1` (`dgCMatrix`), `W_graph2` (`dgCMatrix`).
    *   Output: Spearman correlation `rho` between vectorized non-zero elements (union, zero-filled) of the two graphs' upper triangles.
*   `[X]` **TCK-TSKGR-004:** Implement `residualize_graph_on_subspace` function.
    *   Inputs: `W_graph_to_residualize` (`dgCMatrix`), `L_graph_for_projection` (`dgCMatrix`), `k_eigenvectors_to_remove` (integer, default **64**), `k_nn_resparsify` (integer, k for NN sparsification after residualization).
    *   Internals: Compute first `k_eigenvectors_to_remove` smallest eigenvectors `U` of `L_graph_for_projection`. Compute residual `W_res = W_graph - U(U^T W_graph) - (W_graph U)U^T + U(U^T W_graph U)U^T`. Symmetrize `W_res = (W_res + t(W_res))/2`. Re-sparsify `W_res` using `k_nn_resparsify`-NN method. Re-z-score the final sparse `W_res`.
    *   Output: Residualized, sparse, symmetric, z-scored `W_graph_res` (`dgCMatrix`).
*   `[X]` **TCK-TSKGR-005:** (Potentially skip if direct arithmetic used) Implement `blend_laplacians` function.
    *   Inputs: `L_conn`, `L_task`, `lambda_blend_value`, `method` (e.g., "linear", "geo").
    *   Output: `L_hybrid` (`dgCMatrix`).
*   `[X]` **TCK-TSKGR-006:** Confirm `compute_graph_laplacian_sparse` uses `alpha_laplacian` correctly (likely already does).

**Module: Graph Construction & Blending (Advanced Task-Specific - Intra-Parcel Inputs)**

*   `[ ]` **TCK-TSKADVGR-001:** Implement `compute_multivar_encoding_weights_ridge` function.
    *   Inputs: `X_stimulus` (`T x M`), `Y_p_voxels` (`T x Vp'` - voxels for one parcel), `lambda_ridge`.
    *   Output: Weight matrix `W_p_encoding` (`M x Vp'`).
    *   Note: Likely external pre-processing.
*   `[ ]` **TCK-TSKADVGR-002:** Implement `summarize_encoding_signature` function.
    *   Inputs: `W_p_encoding` (`M x Vp'`), `method` ("mean", "pca", "svd", etc.), `k_reduce`.
    *   Output: Summarized signature `w̃_p` (`M x k_reduce` or `M x 1`).
*   `[ ]` **TCK-TSKADVGR-003:** Implement `compute_W_task_from_intra_parcel_patterns` function.
    *   Inputs: `task_data_list_advanced` (providing `Y_p` or `W_p_encoding`), `method_name`, `k_conn_task_pos`, `k_conn_task_neg`, etc.
    *   Internals: Orchestrates calls based on `method_name` (e.g., Method 1: calls helpers, computes similarity, sparsifies, z-scores).
    *   Output: Sparse `W_task_i` (`dgCMatrix`).
*   `[ ]` **TCK-TSKADVGR-004:** (Optional) Implement `compute_rv_coefficient(matrix_A, matrix_B)` function.
*   `[ ]` **TCK-TSKADVGR-005:** (Optional) Implement `compute_cca_similarity_parcel_patterns(Y_p_voxels, Y_q_voxels, n_cca_dims=1)` function.

**Module: Generalized Eigenvalue Decomposition (GEV)**
*   `[X]` **TCK-GEV-001:** Implement `solve_gev_laplacian_primme` function.
    *   Inputs: `L_task`, `L_conn` (both sparse, symmetric `dgCMatrix`), `k_request` (number of eigenpairs), `lambda_max_thresh` (**synced value, e.g., 0.8**), `stability_r_thresh` (e.g., 0.6, for external filtering), `epsilon_reg_L_conn`.
    *   Internals: Use `PRIMME::eigs_sym(A=L_task, B=L_conn_regularized, NEig=k_solve, which="SM", ...)` Filter output based on `abs(lambda) < lambda_max_thresh`.
    *   Output: List containing filtered `vectors`, `values`, `n_converged`, `n_filtered`.
*   `[X]` **TCK-GEV-002:** Implement `compute_gev_spectrum_diagnostics` function.
    *   Input: `Lambda_GEV` vector, `lambda_max_thresh`.
    *   Output: List/JSON with stats.

**Module: Anchor Augmentation & Weighted Alignment**
*   `[X]` **TCK-AUGALI-001:** Implement `project_features_to_spectral_space` function.
    *   Inputs: `feature_matrix` (`C x V_p` or `V_p x C`), `U_basis` (`V_p x k_dims_basis`).
    *   Internals: Handle potential transposition of `feature_matrix`. Perform **orthogonal projection**: check if `crossprod(U_basis)` is identity; if yes, use `crossprod(U_basis, features)`; if no, use `solve(crossprod(U_basis), crossprod(U_basis, features))`. 
    *   Output: `projected_features` (e.g., `k_dims_basis x C` if input was `V_p x C`).
*   `[X]` **TCK-AUGALI-001.5 (NEW):** Implement `residualize_matrix_on_subspace` helper function.
    *   Inputs: `matrix_to_residualize` (e.g., `Z_i`, `C x k`), `subspace_basis_matrix` (e.g., `A_parc_i`, `m x k`).
    *   Internals: Compute `Y_res = Y - X_subspace %*% solve(crossprod(X_subspace), crossprod(X_subspace, Y))`, where `Y` is `matrix_to_residualize` and `X_subspace` is `subspace_basis_matrix`. Ensure dimensions are handled correctly for matrix multiplication.
    *   Output: `matrix_residualized` (same dimensions as `matrix_to_residualize`).
*   `[X]` **TCK-AUGALI-002:** Implement `build_augmented_anchor_matrix` function.
    *   Inputs: `A_parcel_anchors` (`m_parcels x k_dims`), `Z_task_features_projected` (`m_task_features x k_dims`).
    *   Output: `A_augmented` (`(m_parcels + m_task_features) x k_dims`) by `rbind`.
*   `[X]` **TCK-AUGALI-003:** Implement `solve_procrustes_rotation_weighted` function.
    *   Inputs: `A_source`, `T_target`, `omega_weights` (list, e.g. `list(parcel=w_p, condition=w_c)`), `omega_mode` (char: "fixed", "adaptive"), `reliability_scores` (numeric vector/list, needed if `omega_mode == "adaptive"`).
    *   Internals:
        *   If `omega_mode == "fixed"`: Construct diagonal weight matrix `Omega`. Default to parcel=1.0, condition=0.5 if `omega_weights` is `NULL` or unspecified for these. Rescale `Omega` so `trace(Omega)` equals number of rows in `A_source`.
        *   If `omega_mode == "adaptive"`: Calculate condition weights based on `reliability_scores`. Parcel weights typically remain fixed (e.g., 1.0). Construct `Omega` and rescale trace. (Detailed logic for adaptive calculation needs to be specified, possibly based on `R^2_p` values).
        *   Apply weights: `A_w = Omega %*% A_source`, `T_w = Omega %*% T_target`.
        *   Call `solve_procrustes_rotation(A_w, T_w)`.
    *   Output: Rotation matrix `R`.
*   `[X]` **TCK-AUGALI-004:** Modify `perform_gpa_refinement` function (in `R/procrustes_alignment.R`).
    *   Add optional parameters `omega_weights`, `omega_mode`, `reliability_scores_for_adaptive_omega`.
    *   If `omega_weights` or `omega_mode` indicate weighting, call `solve_procrustes_rotation_weighted` internally, passing through necessary parameters.
    *   Ensure template update step remains unweighted or appropriately handled if weights influence the target.

**Module: task_hatsa Main Workflow & Object**
*   `[X]` **TCK-TSKWKFL-001:** Implement main `run_task_hatsa()` function.
    *   Inputs: Include all relevant parameters from Step 1 of Overview.
    *   Logic: Orchestrate calls reflecting the detailed Algorithmic Overview (Steps 2-7), including redundancy check/residualization.
*   `[X]` **TCK-TSKWKFL-001.P (NEW):** Parallelize per-subject processing in `run_task_hatsa()`.
    *   Refactor the per-subject loop (graph construction, Laplacian computation, basis shaping) to use `future.apply::future_lapply`.
    *   Define an internal helper function to encapsulate single-subject processing, returning only necessary outputs to reduce memory transfer.
    *   Ensure necessary packages (`Matrix`, `RSpectra`, `PRIMME`, etc.) are available to parallel workers.
*   `[X]` **TCK-TSKWKFL-002:** Define `task_hatsa_projector` S3 class and constructor.
    *   Class definition: `c("task_hatsa_projector", "hatsa_projector", ...)`
    *   New slots/parameters storage: Ensure all relevant parameters and QC metrics (rho, residualization status) are stored.
*   `[X]` **TCK-TSKWKFL-003:** Implement/Update S3 methods for `task_hatsa_projector` (`print`, `summary`).
*   `[ ]` **TCK-TSKWKFL-004 (Advanced/Optional):** Implement Adaptive Selection logic.
*   `[ ]` **TCK-TSKWKFL-005 (Advanced/Optional):** Implement `compute_alignment_transferability_isc`. Note mandatory check in validation (`median >= -1%`).

**Module: Testing (task_hatsa Specific)**
*   `[ ]` **TCK-TSKTST-001:** Test `compute_W_task_from_activations` and `compute_W_task_from_encoding` (incl. similarity methods).
*   `[ ]` **TCK-TSKTST-002:** Test `compute_graph_correlation` (union/zero-fill Spearman).
*   `[ ]` **TCK-TSKTST-004:** Test `residualize_graph_on_subspace` (low-rank formula, re-sparsify, re-z-score).
*   `[ ]` **TCK-TSKTST-005:** Test `project_features_to_spectral_space` (orthogonal projection logic).
*   `[ ]` **TCK-TSKTST-005.5 (NEW):** Test `residualize_matrix_on_subspace`.
*   `[ ]` **TCK-TSKTST-006:** Test `solve_procrustes_rotation_weighted` and `perform_gpa_refinement` modifications.
    *   Cover fixed default weights with trace rescaling.
    *   Cover user-specified fixed weights.
    *   Cover adaptive omega mode (once reliability metric and adaptation logic are defined).
*   `[ ]` **TCK-TSKTST-007:** End-to-end test `run_task_hatsa` (`lambda_blend`, including redundancy check path).
    *   Include variations: no anchor aug, anchor aug with fixed default omega, anchor aug with residualized condition anchors.
*   `[ ]` **TCK-TSKTST-008:** End-to-end test `run_task_hatsa` (`gev_patch`, including redundancy check path).
    *   Include variations similar to TCK-TSKTST-007.
*   `[ ]` **TCK-TSKTST-009 (Advanced/Optional):** Test Adaptive Selection.
*   `[ ]` **TCK-TSKTST-010 (Advanced/Optional):** Test `compute_alignment_transferability_isc`.

**Module: Testing (Advanced Task-Specific - Intra-Parcel Inputs)**
*   `[ ]` **TCK-TSKADVTST-001:** Test `compute_multivar_encoding_weights_ridge`.
*   `[ ]` **TCK-TSKADVTST-002:** Test `summarize_encoding_signature`.
*   `[ ]` **TCK-TSKADVTST-003:** Test end-to-end `compute_W_task_from_intra_parcel_patterns`.
*   `[ ]` **TCK-TSKADVTST-004:** (Optional) Test `compute_rv_coefficient`.
*   `[ ]` **TCK-TSKADVTST-005:** (Optional) Test `compute_cca_similarity_parcel_patterns`.

This set of tickets provides a granular plan for implementing task_hatsa, including its advanced GEV options and considerations for weighted alignment, within the R package structure, emphasizing alignment with existing code.