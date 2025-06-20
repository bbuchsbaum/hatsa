% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projection_helpers.R
\name{build_augmented_anchor_matrix}
\alias{build_augmented_anchor_matrix}
\title{Build Augmented Anchor Matrix}
\usage{
build_augmented_anchor_matrix(A_parcel_anchors, Z_task_features_projected)
}
\arguments{
\item{A_parcel_anchors}{A numeric matrix representing parcel anchors,
typically with dimensions `m_parcels x k_dims`.}

\item{Z_task_features_projected}{A numeric matrix representing projected task
features (e.g., condition means in the spectral space), typically with
dimensions `m_task_features x k_dims`. Can be `NULL` or have 0 rows if
no task features are to be added.}
}
\value{
An augmented numeric matrix with dimensions
  `(m_parcels + m_task_features) x k_dims`. 
  Row names are combined from inputs (made unique with `make.unique`). 
  Column names are taken from `A_parcel_anchors` if it has them and rows; 
  otherwise from `Z_task_features_projected`. An error is thrown if both have
  differing non-NULL column names.
  If `Z_task_features_projected` is `NULL` or has 0 rows, `A_parcel_anchors` is returned.
  If `A_parcel_anchors` has 0 rows and `Z_task_features_projected` is valid,
  `Z_task_features_projected` is returned.
}
\description{
Combines parcel anchors and projected task features into a single
  augmented anchor matrix for Procrustes alignment.
  If inputs have different numeric types (e.g. integer and double), the result
  is coerced per R's default rules (usually to the more general type).
}
\examples{
m_parcels <- 5; m_task <- 3; k_dims <- 4
A_p <- matrix(rnorm(m_parcels*k_dims), m_parcels, k_dims, dimnames=list(paste0("p",1:m_parcels), paste0("k",1:k_dims)))
Z_t <- matrix(rnorm(m_task*k_dims), m_task, k_dims, dimnames=list(paste0("t",1:m_task), paste0("k",1:k_dims)))
A_aug <- build_augmented_anchor_matrix(A_p, Z_t)
print(dimnames(A_aug))

# Differing colnames should error
Z_t_bad_colnames <- Z_t; colnames(Z_t_bad_colnames) <- paste0("dim",1:k_dims)
try(build_augmented_anchor_matrix(A_p, Z_t_bad_colnames))

}
