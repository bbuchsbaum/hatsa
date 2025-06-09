#' HATSA Workflow Examples
#'
#' @description
#' This documentation provides comprehensive examples of common HATSA workflows,
#' from basic alignment to advanced task-informed analyses.
#'
#' @section Basic Workflow:
#' \preformatted{
#' # Load your data (list of subject matrices)
#' data <- load_fmri_data()  # Returns list of time × voxel matrices
#' 
#' # Get parameter suggestions
#' params <- hatsa_suggest(data)
#' 
#' # Run basic HATSA
#' result <- hatsa(data, components = params$components)
#' 
#' # Extract aligned data
#' aligned_data <- get_aligned_data(result)
#' template <- get_template(result)
#' 
#' # Check quality
#' hatsa_summary(result)
#' plot_hatsa(result, type = "eigenvalues")
#' }
#'
#' @section Task-Informed Workflow:
#' \preformatted{
#' # Load task design matrices
#' task_data <- load_task_designs()  # List of time × condition matrices
#' 
#' # Using lambda blend method
#' result <- task_hatsa(subject_data_list = data,
#'                      task_data_list = task_data,
#'                      task_method = "lambda_blend")
#' 
#' # Or use GEV method with custom options
#' opts <- task_hatsa_opts(lambda_blend_value = 0.2)
#' result <- task_hatsa(subject_data_list = data,
#'                      task_data_list = task_data,
#'                      task_method = "gev_patch",
#'                      opts = opts)
#' 
#' # Analyze task-specific alignment
#' task_metrics <- get_task_alignment_metrics(result)
#' }
#'
#' @section Advanced Anchor Selection:
#' \preformatted{
#' # Manual anchor selection based on ROI
#' roi_indices <- get_roi_voxels("visual_cortex")
#' result <- hatsa(data, anchors = roi_indices)
#' 
#' # Multi-resolution anchor selection
#' anchors <- select_anchors_mra(
#'   U_list = preliminary_decomposition,
#'   n_anchors = 100,
#'   n_resolutions = 5
#' )
#' result <- hatsa(data, anchors = anchors)
#' }
#'
#' @section Preprocessing Integration:
#' \preformatted{
#' # HATSA works best with preprocessed data
#' data_clean <- lapply(data, function(X) {
#'   X <- scale(X)  # Z-score time series
#'   X[is.na(X)] <- 0  # Handle missing data
#'   X
#' })
#' 
#' # Run with custom preprocessing
#' result <- hatsa(data_clean, preset = "accurate")
#' }
#'
#' @section Group Analysis:
#' \preformatted{
#' # Align multiple groups separately
#' result_controls <- hatsa(data[control_idx])
#' result_patients <- hatsa(data[patient_idx])
#' 
#' # Compare alignment quality
#' compare_alignments(result_controls, result_patients)
#' 
#' # Project new subjects to existing space
#' new_aligned <- predict(result_controls, newdata_list = new_subjects)
#' }
#'
#' @section Performance Optimization:
#' \preformatted{
#' # For large datasets, use fast preset
#' result <- hatsa(big_data, preset = "fast")
#' 
#' # For parallel processing (if available)
#' options(hatsa.parallel = TRUE)
#' options(hatsa.cores = 4)
#' result <- hatsa(data)
#' 
#' # For very high-dimensional data
#' # First reduce dimensions
#' data_reduced <- lapply(data, function(X) {
#'   svd_X <- svd(X, nu = 100, nv = 100)
#'   svd_X$u %*% diag(svd_X$d[1:100])
#' })
#' result <- hatsa(data_reduced)
#' }
#'
#' @name hatsa-workflows
#' @aliases HATSA-workflows hatsa-examples
NULL

#' Common HATSA Pitfalls and Solutions
#'
#' @description
#' Guide to avoiding common issues when using HATSA.
#'
#' @section Pitfall 1 - Too Many Components:
#' Using too many components relative to your data size can lead to overfitting.
#' 
#' Solution: Use hatsa_suggest() or keep components < 10% of voxels.
#'
#' @section Pitfall 2 - Mismatched Dimensions:
#' All subjects must have the same number of voxels (spatial alignment).
#' 
#' Solution: Ensure all subjects are in the same space/parcellation.
#'
#' @section Pitfall 3 - Including Bad Subjects:
#' Subjects with excessive motion or artifacts can degrade alignment.
#' 
#' Solution: Pre-screen subjects and exclude outliers.
#'
#' @section Pitfall 4 - Wrong Preprocessing:
#' HATSA expects centered data (zero mean per voxel).
#' 
#' Solution: Always scale/center your data before HATSA.
#'
#' @section Pitfall 5 - Ignoring Convergence Warnings:
#' Warnings about convergence indicate potential issues.
#' 
#' Solution: Check your data quality and try different parameters.
#'
#' @name hatsa-pitfalls
NULL