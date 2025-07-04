% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mds_helpers.R
\name{prepare_mds_plot_df}
\alias{prepare_mds_plot_df}
\title{Prepare Data for MDS Plot}
\usage{
prepare_mds_plot_df(
  mds_fit,
  valid_subject_indices,
  valid_subjects_mask,
  n_total_subjects,
  subject_info = NULL,
  plot_labels = FALSE
)
}
\arguments{
\item{mds_fit}{Result from `run_cmdscale_safe`.}

\item{valid_subject_indices}{Indices of subjects included in the analysis.}

\item{valid_subjects_mask}{Logical vector identifying valid subjects in the
original distance matrix.}

\item{n_total_subjects}{Total number of subjects in the projector object.}

\item{subject_info}{Optional data frame with one row per subject.}

\item{plot_labels}{Logical indicating if labels should be prepared.}
}
\value{
A data frame suitable for `ggplot2`.
}
\description{
Creates a data frame of MDS coordinates merged with optional
subject information for plotting.
}
\keyword{internal}
