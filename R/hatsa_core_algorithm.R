#' Run the Core HATSA Algorithm
#'
#' Provides a backwards-compatible wrapper for running the standard
#' HATSA alignment without any task information. Internally this
#' calls \code{task_hatsa()} with \code{task_method = "core_hatsa"} and
#' returns a \code{hatsa_projector} object.
#'
#' @param subject_data_list List of subject time-series matrices.
#' @param anchor_indices Numeric vector of anchor parcel indices.
#' @param spectral_rank_k Integer spectral rank to retain.
#' @param k_conn_pos Integer number of positive connections for k-NN sparsification.
#' @param k_conn_neg Integer number of negative connections for k-NN sparsification.
#' @param n_refine Integer number of GPA refinement iterations.
#' @param use_dtw Logical, currently ignored. Included for API compatibility.
#' @param n_cores Integer number of cores (unused; parallelism is not implemented here).
#'
#' @return An object of class \code{hatsa_projector}.
#' @export
run_hatsa_core <- function(subject_data_list,
                           anchor_indices,
                           spectral_rank_k,
                           k_conn_pos,
                           k_conn_neg,
                           n_refine,
                           use_dtw = FALSE,
                           n_cores = 1L) {
  opts <- task_hatsa_opts(
    lambda_blend_value = 0,
    row_augmentation = FALSE,
    k_conn_pos = k_conn_pos,
    k_conn_neg = k_conn_neg,
    n_refine = n_refine
  )

  proj <- task_hatsa(
    subject_data_list = subject_data_list,
    anchor_indices = anchor_indices,
    spectral_rank_k = spectral_rank_k,
    task_data_list = NULL,
    task_method = "core_hatsa",
    opts = opts,
    verbose = TRUE
  )

  proj$method <- "hatsa_core"
  if (is.list(proj$parameters)) proj$parameters$method <- "hatsa_core"
  class(proj) <- c("hatsa_projector",
                   setdiff(class(proj), "task_hatsa_projector"))
  proj
}
