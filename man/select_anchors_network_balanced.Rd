% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anchor_selection_mra_enhanced.R
\name{select_anchors_network_balanced}
\alias{select_anchors_network_balanced}
\title{Network-Informed Anchor Selection}
\usage{
select_anchors_network_balanced(
  U_original_list_pilot,
  spectral_rank_k,
  m_target,
  network_labels,
  anchors_per_network = 2,
  ...
)
}
\arguments{
\item{U_original_list_pilot}{Pilot subject sketch matrices}

\item{spectral_rank_k}{Number of spectral components}

\item{m_target}{Target number of anchors}

\item{network_labels}{Vector of network labels for each parcel}

\item{anchors_per_network}{Minimum anchors per network}

\item{...}{Additional arguments for select_anchors_mra}
}
\value{
Selected anchor indices
}
\description{
Select anchors that represent different functional networks or
brain regions.
}
