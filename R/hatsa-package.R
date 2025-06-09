#' @keywords internal
"_PACKAGE"

#' hatsa: Core HATSA for Functional Connectivity Alignment
#'
#' The `hatsa` package implements the Core HATSA (Hybrid Anchor-based
#' Time-Series Alignment) algorithm. This algorithm is designed to align
#' functional connectivity patterns, represented as spectral sketches derived
#' from graph Laplacians, across multiple subjects. It leverages anchor parcels
#' and Generalized Procrustes Analysis for robust alignment, employing sparse
#' matrix operations for computational efficiency.
#'
#' @section Core HATSA Algorithm:
#' The primary function is \code{\link{hatsa}}. It processes
#' subject-specific time-series data and alignment parameters to produce:
#' \itemize{
#'   \item Original spectral sketches for each subject.
#'   \item Aligned spectral sketches for each subject.
#'   \item The corresponding rotation matrices (in SO(k)) used for alignment.
#' }
#'
#' The algorithm follows three main stages:
#' 1.  **Initial Spectral Sketching (per subject):**
#'     \itemize{
#'       \item A correlation matrix is computed from the time-series.
#'       \item This graph is sparsified based on strongest positive/negative
#'             connections per parcel.
#'       \item The sparsified graph is symmetrized using a weighted averaging scheme.
#'       \item Non-zero edge weights of the symmetrized graph are z-scored,
#'             ensuring symmetry is preserved.
#'       \item A sparse graph Laplacian (`L = D - W`) is constructed.
#'       \item The `k` eigenvectors corresponding to the smallest, non-zero
#'             eigenvalues of `L` form the subject's original spectral sketch.
#'             Eigenvectors for eigenvalues numerically close to zero are discarded.
#'             This step uses efficient methods for sparse matrices.
#'     }
#' 2.  **Iterative Refinement (Generalized Procrustes Analysis - GPA):**
#'     \itemize{
#'       \item Rows corresponding to pre-defined anchor parcels are extracted from
#'             each subject's original spectral sketch.
#'       \item A group anchor template is initialized (typically as the mean of
#'             subjects' anchor sketches).
#'       \item Iteratively:
#'             \itemize{
#'               \item Subject-specific orthogonal rotation matrices (`R_i \in SO(k)`)
#'                     are computed to best align each subject's anchor sketch to the
#'                     current group template (Orthogonal Procrustes Problem).
#'               \item The group anchor template is updated as the mean of the
#'                     subjects' rotated anchor sketches.
#'             }
#'     }
#' 3.  **Apply Final Rotations:**
#'     \itemize{
#'       \item The final rotation matrices from the GPA are applied to each
#'             subject's full original spectral sketch to obtain the aligned
#'             spectral sketches.
#'     }
#'
#' @docType _PACKAGE
#' @name hatsa-package
#' @aliases hatsa hatsaR
NULL