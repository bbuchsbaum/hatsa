% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hatsa_projector.R
\name{hatsa_projector}
\alias{hatsa_projector}
\title{HATSA Projector Object}
\usage{
hatsa_projector(hatsa_core_results, parameters)
}
\arguments{
\item{hatsa_core_results}{A list containing the core outputs from the HATSA
algorithm. Expected elements include:
\itemize{
  \item{\code{U_aligned_list}: List of subject-specific aligned sketch matrices (V_p x k).}
  \item{\code{R_final_list}: List of subject-specific rotation matrices (k x k).}
  \item{\code{U_original_list}: List of subject-specific original sketch matrices (V_p x k).}
  \item{\code{Lambda_original_list}: List of subject-specific original eigenvalues (length k).}
  \item{\code{Lambda_original_gaps_list}: List of subject-specific eigengap ratios (length k-1).}
  \item{\code{T_anchor_final}: The final group anchor template matrix (N_anchors x k).}
}}

\item{parameters}{A list of parameters used to run HATSA. Expected elements include:
\itemize{
  \item{\code{k}: The number of spectral components (rank).}
  \item{\code{N_subjects}: The number of subjects.}
  \item{\code{V_p}: The number of parcels/vertices per subject.}
  \item{\code{method}: A string, typically \code{hatsa_core}.}
}}
}
\value{
A `hatsa_projector` object.

An object of class \code{c("hatsa_projector", "multiblock_biprojector", "projector", "list")}.
  This object contains:
  \itemize{
    \item{\code{v}: The group-level loading matrix (mean aligned sketch, V_p x k).}
    \item{\code{s}: The stacked scores matrix (concatenated aligned sketches, (N*V_p) x k).}
    \item{\code{sdev}: Component standard deviations (defaulted to 1s, length k).}
    \item{\code{preproc}: A preprocessing object, set to \code{multivarious::prep(multivarious::pass())}.}
    \item{\code{block_indices}: A list defining rows in \code{s} corresponding to each subject block.}
    \item{\code{R_final_list}: Stored from input.}
    \item{\code{U_original_list}: Stored from input.}
    \item{\code{Lambda_original_list}: Stored from input (crucial for voxel projection).}
    \item{\code{Lambda_original_gaps_list}: Stored from input.}
    \item{\code{T_anchor_final}: Stored from input.}
    \item{\code{parameters}: Stored from input.}
    \item{\code{method}: Stored from input parameters, typically \code{hatsa_core}.}
  }
}
\description{
An S3 object of class \code{hatsa_projector} that stores the results of a
Harmonized Tensors SVD Alignment (HATSA) analysis. This object inherits from
\code{multiblock_biprojector} (from the `multivarious` package) and is designed
to integrate HATSA outputs into a common framework for multiblock data analysis.

Creates a \code{hatsa_projector} object, which stores the results of the
HATSA algorithm and is designed to integrate with the \code{multivarious}
package, inheriting from \code{multiblock_biprojector}.
}
\section{Fields}{

\describe{
\item{\code{v}}{A numeric matrix (V_p x k) representing the mean aligned sketch,
serving as the group-level template or common loadings.}

\item{\code{s}}{A numeric matrix ((N*V_p) x k) of stacked aligned sketches for all
subjects. These are the subject-specific parcel scores in the common space.}

\item{\code{sdev}}{A numeric vector of length k, representing component-wise standard
deviations (or scales). Currently defaults to a vector of 1s.}

\item{\code{preproc}}{A \code{pre_processor} object (from `multivarious`). For HATSA,
this is typically `prep(pass())` as the input data to `run_hatsa_core` is
already processed up to the point of raw time-series per subject.}

\item{\code{block_indices}}{A list defining which rows in the scores matrix `s` belong
to which subject (block).}

\item{\code{R_final_list}}{A list of subject-specific rotation matrices (k x k) used to
align each subject's original sketch to the common space.}

\item{\code{U_original_list}}{A list of subject-specific original (unaligned) sketch
matrices (V_p x k) derived from their parcel-level graph Laplacians.}

\item{\code{Lambda_original_list}}{A list of numeric vectors, where each vector contains
the k original eigenvalues corresponding to the eigenvectors in `U_original_list`
for that subject. These are crucial for Nyström voxel projection.}

\item{\code{Lambda_original_gaps_list}}{A list of numeric vectors. Each vector contains
the k-1 eigengap ratios `(λ_{j+1} - λ_j) / λ_j` for the corresponding subject's
original eigenvalues. Useful for assessing spectral stability.}

\item{\code{T_anchor_final}}{A numeric matrix (V_a x k, where V_a is the number of
anchors) representing the final group anchor template after Procrustes alignment.}

\item{\code{parameters}}{A list containing the input parameters used for the HATSA run
(e.g., `k`, `V_p`, `N_subjects`, `anchor_indices`, `k_conn_pos`, `k_conn_neg`, `n_refine`).}

\item{\code{method}}{A character string, typically "hatsa_core", indicating the method
used to generate the projector.}

\item{\code{U_aligned_list}}{(Internal) A list of subject-specific aligned sketch matrices (V_p x k).
While `s` provides the stacked version, this list might be retained internally from the
`run_hatsa_core` output passed to the constructor. For user access to aligned sketches per subject,
one would typically use `project_block(object, block = i)` or segment `scores(object)`
using `block_indices(object)`.}
}}

\examples{
# This is a conceptual example, as real data structures are complex.
# Assuming hatsa_results and params are populated from run_hatsa_core:
# projector_obj <- hatsa_projector(
#   hatsa_core_results = list(
#     U_aligned_list = replicate(5, matrix(rnorm(100*10), 100, 10), simplify=FALSE),
#     R_final_list = replicate(5, diag(10), simplify=FALSE),
#     U_original_list = replicate(5, matrix(rnorm(100*10), 100, 10), simplify=FALSE),
#     Lambda_original_list = replicate(5, runif(10, 0.1, 1), simplify=FALSE), # example
#     Lambda_original_gaps_list = replicate(5, runif(9, 0.05, 0.5), simplify=FALSE), # example
#     T_anchor_final = matrix(rnorm(5*10), 5, 10)
#   ),
#   parameters = list(
#     k=10,
#     N_subjects=5,
#     V_p=100,
#     method="hatsa_core"
#   )
# )
# class(projector_obj)
# names(projector_obj)
}
\seealso{
\code{\link{run_hatsa_core}}, \code{\link{predict.hatsa_projector}}, \code{\link{project_voxels.hatsa_projector}}
}
