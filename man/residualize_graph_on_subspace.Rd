% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/task_graph_construction.R
\name{residualize_graph_on_subspace}
\alias{residualize_graph_on_subspace}
\title{Residualize Graph B based on Subspace from Graph A's Laplacian}
\usage{
residualize_graph_on_subspace(
  W_graph_to_residualize,
  L_graph_for_projection,
  k_eigenvectors_to_remove = 64,
  k_nn_resparsify,
  eigenvalue_tol = 1e-09
)
}
\arguments{
\item{W_graph_to_residualize}{The sparse graph matrix (`dgCMatrix`) to be
residualized (e.g., `W_task`).}

\item{L_graph_for_projection}{The sparse Laplacian matrix (`dgCMatrix`) from
which the projection subspace is derived (e.g., `L_conn`). Must be symmetric.}

\item{k_eigenvectors_to_remove}{Integer, the number of smallest (by magnitude)
eigenvectors of `L_graph_for_projection` to define the subspace for projection.
Defaults to 64.}

\item{k_nn_resparsify}{Integer, the k value for k-NN sparsification applied
to the residual graph (same k used for positive and negative edges).}

\item{eigenvalue_tol}{Numeric, tolerance for eigenvalue decomposition convergence
and for identifying near-zero eigenvalues if needed (though projection uses the space).
Default 1e-9.}
}
\value{
A sparse, symmetric, z-scored `dgCMatrix` representing the residualized graph.
}
\description{
Projects `W_graph_to_residualize` onto the subspace spanned by the first
`k_eigenvectors_to_remove` smallest eigenvectors of
`L_graph_for_projection` and subtracts this projection.
The residual graph is then re-sparsified and re-z-scored.
}
\keyword{internal}
