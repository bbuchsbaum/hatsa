% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/riemannian_geometry.R
\name{airm_distance}
\alias{airm_distance}
\title{Affine-Invariant Riemannian Metric (AIRM) Distance}
\usage{
airm_distance(S1, S2, regularize_epsilon = 1e-06)
}
\arguments{
\item{S1}{A numeric, symmetric positive-definite matrix.}

\item{S2}{A numeric, symmetric positive-definite matrix.}

\item{regularize_epsilon}{Epsilon for regularization. Default 1e-6.}
}
\value{
The AIRM distance (a non-negative scalar).
}
\description{
Computes the Affine-Invariant Riemannian Metric (AIRM) distance between two
symmetric positive-definite (SPD) matrices S1 and S2.
The distance is defined as: `||logm(S1^(-1/2) %*% S2 %*% S1^(-1/2))||_F`,
where `||.||_F` is the Frobenius norm.
}
\examples{
S1 <- matrix(c(2.3, -0.3, -0.3, 3.6), 2, 2)
S2 <- matrix(c(3.7, 1.9, 1.9, 2.8), 2, 2)
# dist_airm <- airm_distance(S1, S2)
# print(dist_airm)
# dist_logeuclidean <- riemannian_distance_spd(S1, S2)
# print(dist_logeuclidean) # Note: AIRM and LogEuclidean are different metrics.
}
