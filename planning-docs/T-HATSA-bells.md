# HATSA Core Metrics & Anchor Optimization Support Tickets

**Overall Goal:** To provide a rich set of generally useful diagnostic metrics and component extractions from fitted HATSA objects, which can secondarily be used to inform and evaluate anchor selection strategies, without embedding anchor optimization logic directly into the core fitting routines.

---

## **Ticket HMET-001: Enhance `summary.hatsa_projector` with Detailed Alignment & Variance Metrics**

- [ ] **Description:** The current `summary.hatsa_projector` computes mean anchor alignment error. This ticket proposes expanding it to include more comprehensive metrics about the alignment quality and the spectral sketches.
- [ ] **Specific Metrics to Add:**
    - [ ] **Rotation Dispersion (`σ_R`):**
        - Calculate the Fréchet mean `R̄` of all subject rotation matrices `{R_i}` on SO(k).
        - Report the mean and standard deviation of geodesic distances `d(R_i, R̄)` across subjects. This quantifies how much individual rotations deviate from the group mean rotation.
        - *Vetting Commentary (Original):* Implementation detail: do not iterate a full Karcher mean each call. Cache `R̄` inside the projector object once computed (e.g., `object$._cache$R_bar`). For Large N (N ≫ 500), Fréchet mean and `σ_R` will be the slowest part. Use the chordal‐mean approximation (`R̃ = UΣVᵀ`, where `Σ = diag(sign(det(UVᵀ))))` with one refinement-step; error <0.01 rad for k≤40.
        - *Developer Notes & Vetted Implementation for Fréchet Mean on SO(k) (HMET-001 Update):*
            - **Context & Recommendations Summary:**
                - The iterative Karcher mean (matrix logarithm -> average in tangent space -> matrix exponential) is mathematically sound for SO(k).
                - **Performance is critical**: `R_bar` (the Fréchet mean) **must not** be recomputed on every `summary()` call. It should be **cached** (e.g., computed once during object creation or on first access via a getter).
                - The **default computation** should be a fast chordal mean (SVD-based) followed by one Karcher refinement step. This approach offers an excellent balance of speed and accuracy (RMS error typically < 10⁻³ radians for k≤50 with one refinement).
                - A full iterative Karcher mean should be an **optional alternative** for users requiring maximum precision.
                - Relying on `expm` (an existing dependency) for a self-contained helper function is preferable to introducing a new package dependency like `frechet` for this specific SO(k) mean, as the latter might not offer significant performance benefits here due to its internal data conversions.
            - **Recommended Implementation Strategy:**
                1.  **Helper Functions**: Implement `frechet_mean_so_fast` (chordal mean + optional single refinement, defaults to full Karcher if one step is insufficient) and `frechet_mean_so_karcher` (full iterative Karcher mean) as provided below. These should be placed in a suitable utility file (e.g., `R/math_so_k_helpers.R` or `R/utils.R`).
                2.  **Caching `R_bar`**:
                    - Initialize a `._cache` list within the `hatsa_projector` object (e.g., `proj$._cache <- list()`).
                    - Compute and store `R_bar` using `frechet_mean_so_fast(R_final_list, refine = TRUE)` during object construction (e.g., in `run_hatsa_core` before calling the `hatsa_projector` constructor, or within the constructor itself) and store it, e.g., `proj$._cache$R_frechet_mean`.
                    - Alternatively, compute it on first access via an internal getter function called by `summary.hatsa_projector`.
                3.  **`summary.hatsa_projector` Usage**:
                    - The `summary` method should retrieve the cached `R_bar`.
                    - It could optionally expose an argument like `use_full_karcher_mean = FALSE` to allow users to force recalculation using `frechet_mean_so_karcher`.
            - **Vetted R Helper Functions for SO(k) Fréchet Mean:**
              ```R
              # (To be placed in a utility file like R/math_so_k_helpers.R)

              # -------- Chordal (SVD) Mean with Optional One-Step Refinement ---------
              # Computes the Fréchet mean of a list of SO(k) rotation matrices.
              # Defaults to a fast SVD-based chordal mean plus one Karcher refinement step.
              # If refine=TRUE and one step is insufficient, it falls back to full Karcher.
              #
              # Args:
              #   Rlist: A list of k x k SO(k) rotation matrices.
              #   refine: Logical, whether to perform one Karcher refinement step.
              #   tol: Tolerance for convergence checks.
              #
              # Returns: A k x k matrix representing the Fréchet mean (R_bar).
              frechet_mean_so_fast <- function(Rlist, refine = TRUE, tol = 1e-8) {
                if (!all(sapply(Rlist, function(R) is.matrix(R) && inherits(R, "matrix")))) {
                  stop("All elements in Rlist must be standard R matrices.")
                }
                valid_Rlist <- Filter(Negate(is.null), Rlist)
                if (length(valid_Rlist) == 0) {
                  warning("Rlist is empty or contains only NULLs. Cannot compute Fréchet mean.")
                  # Return an identity matrix of undetermined dimension or handle error appropriately.
                  # k_fallback <- if(length(Rlist)>0 && is.matrix(Rlist[[1]])) nrow(Rlist[[1]]) else 0
                  # return(if(k_fallback > 0) diag(k_fallback) else matrix(numeric(0),0,0))
                  stop("Cannot determine dimension for Fréchet mean with empty or all-NULL Rlist.")
                }
                k <- nrow(valid_Rlist[[1]])
                if (k == 0) return(matrix(numeric(0),0,0)) # Handle k=0 case

                N <- length(valid_Rlist)
                M <- matrix(0, nrow = k, ncol = k) # Sum of all R_i

                for (R_i in valid_Rlist) {
                  if (is.matrix(R_i) && all(dim(R_i) == c(k, k))) {
                    M <- M + R_i
                  } else {
                    warning(sprintf("Skipping non-conformant matrix in Rlist: expected %d x %d.", k, k))
                  }
                }
                M <- M / N # Arithmetic mean

                sv <- svd(M)
                S_diag_vals <- rep(1, k)
                # Ensure determinant is +1 for SO(k) projection
                if (det(sv$u %*% t(sv$v)) < 0) S_diag_vals[k] <- -1
                Rbar_chordal <- sv$u %*% diag(S_diag_vals, nrow = k, ncol = k) %*% t(sv$v)

                if (!refine) return(Rbar_chordal)

                # --- One Karcher refinement step ---
                log_sum_tangent <- matrix(0, nrow = k, ncol = k)
                num_valid_for_refine <- 0
                for (R_i in valid_Rlist) {
                   if (is.matrix(R_i) && all(dim(R_i) == c(k, k))) {
                    A <- Rbar_chordal %*% t(R_i) # Error rotation from R_i to Rbar_chordal (in R_i's frame)
                    # log_A <- expm::logm(A)      # Matrix log in Lie algebra so(k)
                    log_A <- tryCatch(expm::logm(A), error = function(e) {
                        warning("logm failed in frechet_mean_so_fast refinement. Using zero matrix for this tangent vector."); matrix(0,k,k)
                    })
                    log_sum_tangent <- log_sum_tangent + log_A
                    num_valid_for_refine <- num_valid_for_refine + 1
                  }
                }
                
                if(num_valid_for_refine == 0) {
                    warning("No valid matrices for refinement step in frechet_mean_so_fast. Returning chordal mean.")
                    return(Rbar_chordal)
                }
                mean_update_log <- log_sum_tangent / num_valid_for_refine
                
                Rbar_refined <- expm::expm(mean_update_log) %*% Rbar_chordal
                
                if (max(abs(mean_update_log)) < tol) return(Rbar_refined)
                
                # Fall back to full Karcher if one step didn't converge and refinement was requested
                # message("Fréchet mean (fast): One step not converged to tol, falling back to full Karcher.")
                return(frechet_mean_so_karcher(valid_Rlist, R_init = Rbar_refined, tol = tol))
              }

              # -------- Full Karcher Mean (Slow but Accurate) -----------------------
              # Iteratively computes the Karcher mean of a list of SO(k) matrices.
              #
              # Args:
              #   Rlist: A list of k x k SO(k) rotation matrices.
              #   R_init: Optional initial estimate for R_bar. If NULL, uses chordal mean.
              #   tol: Tolerance for convergence.
              #   maxit: Maximum number of iterations.
              #
              # Returns: A k x k matrix for the Karcher mean.
              frechet_mean_so_karcher <- function(Rlist, R_init = NULL, tol = 1e-10, maxit = 100L) {
                valid_Rlist <- Filter(Negate(is.null), Rlist)
                 if (length(valid_Rlist) == 0) {
                  k_dim_fallback <- if(!is.null(R_init) && is.matrix(R_init)) nrow(R_init) else 0
                  warning("Rlist is empty or contains only NULLs for Karcher mean.")
                  if (k_dim_fallback > 0) return(diag(k_dim_fallback)) else stop("Cannot determine dimension for Karcher mean.")
                }
                k <- nrow(valid_Rlist[[1]])
                if (k == 0) return(matrix(numeric(0),0,0))

                Rbar_current <- if (is.null(R_init)) {
                  frechet_mean_so_fast(valid_Rlist, refine = FALSE) # Initialize with chordal mean
                } else {
                  R_init
                }
                N <- length(valid_Rlist)

                for (it in seq_len(maxit)) {
                  log_sum_tangent <- matrix(0, k, k)
                  num_valid_R <- 0
                  for (R_i in valid_Rlist) {
                    if (is.matrix(R_i) && all(dim(R_i) == c(k, k))) {
                      # A <- Rbar_current %*% t(R_i) # For update R_new = expm(mean(logm(A))) %*% R_current
                      # log_A_i <- expm::logm(A)
                      # For update R_new = R_current %*% expm(mean(logm(t(R_current)%*%R_i))):
                      err_rot_in_tangent_at_I <- t(Rbar_current) %*% R_i
                      # log_A_i <- expm::logm(err_rot_in_tangent_at_I)
                      log_A_i <- tryCatch(expm::logm(err_rot_in_tangent_at_I), error = function(e) {
                          warning(sprintf("logm failed in Karcher iteration %d. Using zero matrix for tangent vector.", it)); matrix(0,k,k)
                      })
                      log_sum_tangent <- log_sum_tangent + log_A_i
                      num_valid_R <- num_valid_R + 1
                    }
                  }
                  if(num_valid_R == 0) { # Should have been caught by initial check if N > 0
                      warning("No valid matrices left in Rlist during Karcher iteration. Returning current Rbar.")
                      return(Rbar_current)
                  }
                  
                  mean_tangent_update <- log_sum_tangent / num_valid_R
                  
                  # Rbar_new <- expm::expm(mean_tangent_update) %*% Rbar_current # If A = Rbar %*% t(R_i)
                  Rbar_new <- Rbar_current %*% expm::expm(mean_tangent_update)   # If A = t(Rbar) %*% R_i

                  if (max(abs(mean_tangent_update)) < tol) {
                    Rbar_current <- Rbar_new
                    break
                  }
                  Rbar_current <- Rbar_new
                  if (it == maxit) warning("Karcher mean did not converge in ", maxit, " iterations.")
                }
                return(Rbar_current)
              }
              ```
    - [ ] **Variance Explained by `k` Components (per subject):**
        - If the original graph Laplacian `L_i` and its full spectrum were available (or estimated via trace), calculate the proportion of total graph variance (sum of all eigenvalues of `L_i`) captured by the `k` selected eigenvalues `Λ_original_list[[i]]`.
        - *Alternative (if full spectrum is too costly):* Report the sum/mean of the `k` selected eigenvalues `Λ_original_list[[i]]` as an indicator of captured "energy."
        - *Vetting Commentary:* Clarify that you'll compute `∑_{j≤k} λ_j / trace(L_i)` only if the trace was persisted during fitting; otherwise fall back to "energy in first k eigenvalues".
    - [ ] **Eigengap Ratios (per subject):**
        - For each subject\'s `Lambda_original_list[[i]]`, calculate `gap_ratio_j = (λ_{j+1} - λ_j) / λ_j` for `j=1..k-1`.
        - Report summary statistics (e.g., median, min) of `gap_ratio_{k-1}` (the gap before the last component) or `gap_ratio_k` (if `k+1` eigenvalues were temporarily computed).
        - *Vetting Commentary:* Good. Compute once during fit (cheap) and store in `hatsa_obj$Lambda_original_list_gaps` (or similar) to avoid recomputation in `summary()`.
- [ ] **Implementation:**
    - Modify `summary.hatsa_projector` to compute and include these new fields.
    - Helper functions might be needed for Fréchet mean on SO(k) and geodesic distance (consider Rcpp for performance).
- [ ] **General Usefulness:** Provides deeper insight into alignment consistency and the quality of the spectral representation for any HATSA run.
- [ ] **Anchor Optimization Use:** `σ_R` is a direct input. Eigengap ratios and variance explained help assess if the chosen `k` and `V_p` (which influence anchors) are stable.
- [ ] **Backward Compatibility Note:** `summary.hatsa_projector` gets new fields. Ensure downstream code that consumes its output (e.g., by specific named fields or column order if `as.data.frame` is used) is robust to these additions or updated.

---

## **Ticket HMET-002: S3 Method `reconstruction_error()` for `hatsa_projector`**

- [ ] **Description:** Implement an S3 method `reconstruction_error(object, type = "anchors" | "all_parcels" | "non_anchors", ...)` to quantify how well the aligned spectral sketches (or a subset of them, like anchors) can reconstruct a target, or how well non-anchors can be predicted from anchors.
- [ ] **Specific Metrics & Types:**
    - [ ] **`type = "anchors"` (Default):** This is essentially what `summary.hatsa_projector` already calculates (Frobenius error of `A_orig_i @ R_i` vs `T_anchor_final`). This method would provide it directly, perhaps with more detail (e.g., per-anchor-row error).
    - [ ] **`type = "all_parcels"`:**
        - Reconstruct each subject\'s original full sketch `U_orig_i` from their aligned sketch `U_aligned_i` using the inverse rotation: `U_reconstructed_i = U_aligned_i %*% t(R_i)`.
        - Report `||U_orig_i - U_reconstructed_i||_F`. (Should be near zero, a sanity check on `R_i` orthogonality).
        - *Vetting Commentary:* Strictly returns ~0 by construction (`U_aligned %*% Rᵀ = U_orig`). Still useful as a smoke-test, but consider wrapping it in `if (!object$parameters$debug_mode) return(NULL)` or similar to skip by default unless explicitly requested.
    - [ ] **`type = "non_anchors"` (Anchor Residual `ε̄`):**
        - For each subject `i`:
            - Let `U_A_aligned_i = U_aligned_i[A, :]` (aligned anchor rows).
            - Let `U_NA_aligned_i = U_aligned_i[-A, :]` (aligned non-anchor rows).
            - Predict `U_NA_aligned_i_pred` from `U_A_aligned_i` using a chosen method. The proposal mentioned "barycentric weights." This implies finding weights `W_bary` such that `U_NA_orig_i ≈ W_bary %*% U_A_orig_i` (learned from `U_orig_i` on pilot data or from group mean `v`), then applying these weights in aligned space: `U_NA_aligned_i_pred = W_bary %*% U_A_aligned_i`.
            - Report `||U_NA_aligned_i - U_NA_aligned_i_pred||_F^2` averaged across non-anchor parcels and subjects.
            - *Vetting Commentary:* Solid. Provide two weighting modes: (a) barycentric (default) and (b) ridge-regression (e.g., `lambda=1e-4`) from `U_A_aligned_i` to `U_NA_aligned_i` (or `U_A_orig_i` to `U_NA_orig_i` if pilot data is used) to keep things stable when anchors are nearly collinear.
- [ ] **Implementation:**
    - Create a new S3 generic `reconstruction_error` and its method `reconstruction_error.hatsa_projector`.
    - Helper for barycentric/regression weight calculation might be needed.
- [ ] **Return Value:**
    - *Vetting Commentary:* Return a named list with elements like `per_subject_error`, `mean_error`, and optionally `per_parcel_error` (if requested for `type = "non_anchors"`). This structure facilitates easy plotting (e.g., violin plots) in downstream analysis notebooks.
- [ ] **General Usefulness:** `type="anchors"` is a direct alignment quality check. `type="all_parcels"` is a sanity check. `type="non_anchors"` assesses how well the chosen `k`-dimensional space and anchor set capture information in the rest of the brain/parcels.
- [ ] **Anchor Optimization Use:** `ε̄` (from `type="non_anchors"`) is a direct input.

---

## **Ticket HMET-003: Expose Subject-Specific Anchor Sketch Matrices**

- [ ] **Description:** While `T_anchor_final` (the group template for anchors) is stored, the subject-specific *original* anchor sketches `A_orig_i = U_original_list[[i]][anchor_indices, ]` and *aligned* anchor sketches `A_aligned_i = U_aligned_list[[i]][anchor_indices, ]` are not directly top-level outputs, though they are computed internally during GPA. Make these easily accessible.
- [ ] **Implementation Options:**
    1.  Add to `hatsa_projector` object: Store `A_original_list_anchors` and `A_aligned_list_anchors` as new top-level lists in the projector object.
    2.  S3 Extractor Method: Create `get_anchor_sketches(object, type = "original" | "aligned", subject_idx = NULL, cache_in_object = FALSE)`.
- [ ] **Recommendation & Vetting Commentary:** Option 2 (S3 extractor) is cleaner. Add `cache_in_object = TRUE` argument. If `TRUE`, the first computation caches the result (e.g., `object$._cache$A_original_list_anchors`) to avoid re-computation. Storing `A_original_list_anchors` (m × k × N numeric) directly in the object could significantly increase memory footprint. Lazy computation with optional caching is preferred.
- [ ] **General Usefulness:** Allows detailed inspection of how individual subjects\' anchor regions are represented and how well they align to the group template.
- [ ] **Anchor Optimization Use:** Needed for computing condition numbers of individual `A_orig_i` or for strategies that analyze anchor properties per subject.

---

## **Ticket HMET-004: S3 Method `get_rotations()` for `hatsa_projector`**

- [ ] **Description:** The `R_final_list` is already stored. This ticket is about providing a formal S3 accessor for consistency.
- [ ] **Implementation:**
    - Create S3 generic `get_rotations` and `get_rotations.hatsa_projector(object, subject_idx = NULL, as_array = FALSE, ...)`.
    - If `subject_idx` is `NULL`, return the full list. If an index is given, return that specific rotation matrix.
    - *Vetting Commentary:* Add `as_array = TRUE` argument. If `TRUE`, return a 3D array `[k, k, N_subjects]`.
- [ ] **General Usefulness:** Standardized access to the core transformation outputs.
- [ ] **Anchor Optimization Use:** `R_i` matrices are needed to compute `σ_R`.

---

## **Ticket HMET-005: Store and Access Graph Construction QC in `task_hatsa_projector`**

- [ ] **Description:** `task_hatsa` (specifically `task_hatsa_helpers.R/compute_task_matrices`) already computes `rho_redundancy` between `W_conn` and `W_task_raw`, and a flag `was_residualized`. These are stored in `object$qc_metrics`. This ticket ensures they are easily accessible and potentially expands them.
- [ ] **Specific Metrics:**
    - [ ] `rho_redundancy` (already there).
    - [ ] `was_residualized` (already there).
    - [ ] **(New)** Properties of `W_conn_i` and `W_task_i` (if computed): e.g., sparsity (`nnz / V_p^2`), mean/std of non-zero edge weights *before* z-scoring and *before* `drop0()`.
        - *Vetting Commentary:* These metrics can help diagnose issues if `lambda_blend` behaves unexpectedly (e.g., "all-zeros because of k-conn mis-specification" bugs).
    - [ ] **(New for Hybrid Graph)** If `task_method == "lambda_blend"`: Store the eigengap ratios of the *hybrid Laplacian* `L_hybrid_i` for each subject (similar to HMET-001 for `L_conn`).
        - *Vetting Commentary:* Also store the condition number of `(I - αD⁻¹W_hybrid)`. This is critical for monitoring the stability of the blended basis, as it can indicate problems before the eigengap if `lambda_blend` (or `alpha_task` from proposal) is too high.
- [ ] **Implementation:**
    - Modify `task_hatsa_helpers.R/process_single_subject` (or its callees like `compute_task_matrices` and `shape_basis`) to calculate and return these additional graph/Laplacian properties.
    - Store them within the `qc_metrics` list for each subject in the `task_hatsa_projector` object.
    - Enhance `summary.task_hatsa_projector` to report summaries of these new metrics (e.g., median hybrid eigengap ratio, median cond(I - αD⁻¹W_hybrid)).
- [ ] **Integration Note:**
    - *Vetting Commentary:* For mixed `task_hatsa` / `hatsa_core` analysis pipelines, ensure these `qc_metrics` fields are gracefully handled (e.g., present as `NA` or an empty list) in `hatsa_projector` objects that didn\'t go through `task_hatsa` fitting. Avoid structures that might break `as.data.frame()` or similar operations on combined results.
- [ ] **General Usefulness:** Provides critical diagnostics for the graph construction and basis shaping stages, especially when using label information in `task_hatsa`.
- [ ] **Anchor Optimization Use:** Hybrid eigengap and condition number are key metrics for choosing `lambda_blend`, which might interact with anchor choices.

---

## **Ticket HMET-006: Expose Condition Number of Final Group Anchor Template**

- [ ] **Description:** The condition number `κ(T_anchor_final)` is a good indicator of the stability of the final reference space for anchors.
- [ ] **Implementation:**
    - Add `condition_number_T_anchor_final = pracma::cond(object$T_anchor_final)` to `summary.hatsa_projector`.
    - (The κ-aware greedy selection in Appendix E proposal focuses on `κ` of the *pilot mean anchor sketch* during selection; this ticket is about `κ` of the *final template from the full cohort fit*).
    - *Vetting Commentary:* Good. Add a simple color-coded warning message in the `summary()` print output if `κ > 1e4` (or a configurable threshold).
- [ ] **General Usefulness:** Overall diagnostic of the Procrustes problem\'s stability.
- [ ] **Anchor Optimization Use:** Can be a final check; if `κ` of the template derived from a candidate anchor set is too high on the full normative sample, that anchor set might be problematic despite passing pilot checks.

---

**Summary of How These Tickets Support Anchor Optimization (as per Appendix E Roadmap):**

*   **`ε̄` (Anchor Reconstruction Residual):** HMET-002 (`reconstruction_error(type="non_anchors")`).
*   **`σ_R` (Rotation Dispersion):** HMET-001 (via `summary()`, uses `R_i` from HMET-004 `get_rotations()`).
*   **Downstream Score (Decoding/ISC):** This requires external code that uses the aligned sketches (`object$s` or `project_block(object)` outputs). The HATSA object provides the necessary aligned data.
*   **Laplacian Eigengap:**
    *   For `L_conn` (subject-specific): HMET-001 (via `summary()`, using `Lambda_original_list` and derived `Lambda_original_list_gaps`).
    *   For `L_hybrid` (subject-specific, in `task_hatsa`): HMET-005.
*   **Time-to-Solution:** Can be wrapped externally around `run_hatsa_core` or `run_task_hatsa`.
*   **Condition Number `κ` for anchor selection (during pilot phase):** The κ-aware greedy selection described in Appendix E operates *before* the main HATSA fit, using pilot data. HMET-003 (`get_anchor_sketches(type="original")`) would be useful if operating on individual subject sketches from pilot data.
*   **Condition Number `κ` for final template evaluation:** HMET-006.

These tickets focus on making the necessary *components and summary statistics* available from a standard HATSA fit, allowing external scripts/notebooks to implement anchor optimization experiments without cluttering the core HATSA functions.
