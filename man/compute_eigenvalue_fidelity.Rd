% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hatsa_validation_metrics.R
\name{compute_eigenvalue_fidelity}
\alias{compute_eigenvalue_fidelity}
\title{Compute Eigenvalue Spectrum Fidelity Metrics}
\usage{
compute_eigenvalue_fidelity(
  hatsa_object,
  true_eigenvalues_list,
  k_to_compare = NULL
)
}
\arguments{
\item{hatsa_object}{A fitted `hatsa_projector` object, which contains
`Lambda_original_list` (list of numeric vectors of eigenvalues).}

\item{true_eigenvalues_list}{A list of numeric vectors, where each vector contains
the true eigenvalues for the corresponding subject. If a single vector is provided,
it is assumed to be the common true eigenvalues for all subjects.}

\item{k_to_compare}{Integer or NULL. The number of top eigenvalues to compare.
If NULL (default), all available eigenvalues are compared (up to the length of
the shorter of the estimated or true eigenvalue vectors for that subject).}
}
\value{
A list of lists, one for each subject. Each inner list contains:
  \item{correlation}{Pearson correlation between estimated and true eigenvalues.}
  \item{mse}{Mean Squared Error between estimated and true eigenvalues.}
  \item{num_compared}{The number of eigenvalue pairs actually compared.}
  Returns NULL for a subject if inputs are invalid for that subject.
}
\description{
Compares subject-specific eigenvalues from HATSA (`Lambda_original_list`)
with true eigenvalues (e.g., from a toy data generator).
}
\examples{
# See core-hatsa-toy-example.Rmd for usage.
}
