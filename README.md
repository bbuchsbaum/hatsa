# HATSA: Hyperalignment via Task-informed Shared Analysis

[![R-CMD-check](https://github.com/bbuchsbaum/hatsa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/hatsa/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/bbuchsbaum/hatsa/branch/main/graph/badge.svg)](https://app.codecov.io/gh/bbuchsbaum/hatsa?branch=main)

## Overview

HATSA (Hyperalignment via Task-informed Shared Analysis) is an R package that implements advanced methods for aligning functional neuroimaging data across subjects. The package addresses the fundamental challenge in functional neuroimaging of establishing correspondence between brain regions across individuals, enabling group-level analyses while preserving subject-specific functional organization.

The core innovation of HATSA lies in its integration of task-informed features with graph-based spectral alignment, providing a principled approach to functional hyperalignment that leverages both resting-state connectivity and task-evoked activation patterns.

## Key Features

### Core HATSA Algorithm
- **Sparse graph-based alignment** using spectral decomposition of connectivity Laplacians
- **Generalized Procrustes Analysis (GPA)** for optimal orthogonal transformations
- **Anchor-based alignment** ensuring robust correspondence across subjects
- **SO(k) manifold operations** with geodesic distance computations

### Task-informed Extensions
- **Lambda blending**: Weighted combination of connectivity and task-activation graphs
- **GEV patches**: Generalized eigenvalue-based incorporation of task information
- **Residualization methods** for handling redundancy between connectivity and task features

### Advanced Capabilities
- **Nyström extension** for projecting high-dimensional voxel data
- **Riemannian geometry operations** on symmetric positive definite (SPD) manifolds
- **Quality control metrics** and alignment validation
- **Multivariate projection framework** integration

## Installation

```r
# Install from GitHub
devtools::install_github("bbuchsbaum/hatsa")
```

## Quick Start

### Basic HATSA Alignment

```r
library(hatsa)

# Core HATSA with automatic parameter selection
result <- hatsa(subject_data_list)

# Extract aligned data and examine results
aligned_data <- get_aligned_data(result)
summary(result)
```

### Task-informed HATSA

```r
# Task-informed alignment with lambda blending
task_result <- hatsa_task(
  subject_data_list = connectivity_data,
  task_data_list = activation_maps,
  method = "lambda_blend",
  lambda = 0.3
)

# GEV-based task integration
gev_result <- hatsa_task(
  subject_data_list = connectivity_data,
  task_data_list = activation_maps,
  method = "gev_patch",
  k_gev_dims = 5
)
```

### Advanced Usage

```r
# Detailed control over hyperalignment parameters
result <- run_hatsa_core(
  subject_data_list = timeseries_data,
  anchor_indices = c(1, 15, 23, 45, 67),
  spectral_rank_k = 20,
  k_conn_pos = 7,
  k_conn_neg = 7,
  n_refine = 3
)

# Project high-dimensional voxel data
voxel_projections <- project_voxels(
  result, 
  voxel_timeseries_list = voxel_data,
  voxel_coords = voxel_coordinates,
  parcel_coords = parcel_centroids
)
```

## Mathematical Foundation

HATSA implements a spectral approach to functional alignment based on the eigendecomposition of graph Laplacians derived from functional connectivity. The core algorithm:

1. **Graph Construction**: Builds sparse k-NN graphs from subject-specific connectivity matrices
2. **Spectral Decomposition**: Computes low-rank spectral sketches via Laplacian eigendecomposition
3. **Anchor Alignment**: Establishes correspondence using anchor parcels across subjects
4. **Procrustes Refinement**: Iteratively refines alignment via generalized Procrustes analysis
5. **Projection**: Maps data to the common representational space

The task-informed extensions incorporate activation patterns through:
- **Weighted graph combination** (lambda blending)
- **Generalized eigenvalue problems** (GEV patches)
- **Manifold-based integration** on SPD matrices

## Core Functions

### Primary Interfaces
- `hatsa()`: High-level interface for core HATSA alignment
- `hatsa_task()`: Task-informed hyperalignment
- `run_hatsa_core()`: Low-level core algorithm control
- `run_task_hatsa()`: Detailed task-informed alignment

### Spectral Graph Methods
- `compute_subject_connectivity_graph_sparse()`: Sparse connectivity graphs
- `compute_graph_laplacian_sparse()`: Laplacian computation
- `compute_spectral_sketch_sparse()`: Eigendecomposition
- `solve_gev_laplacian_primme()`: Generalized eigenvalue solver

### Alignment and Projection
- `solve_procrustes_rotation()`: Orthogonal alignment
- `project_voxels()`: Nyström extension for voxel projection
- `misalign_deg()`: SO(k) geodesic distances

### Quality Control
- `hatsa_summary()`: Alignment quality metrics
- `plot_hatsa()`: Visualization methods
- `reconstruction_error()`: Validation metrics

## Dependencies

### Required
- `Matrix`: Sparse matrix operations
- `PRIMME`: Generalized eigenvalue problems
- `expm`: Matrix exponential/logarithm operations
- `multivarious`: S3 projector framework
- `RANN`: Fast nearest neighbor search

### Suggested
- `vegan`: Procrustes analysis
- `ggplot2`: Visualization
- `testthat`: Testing framework

## Mathematical Details

The package implements computations on several mathematical structures:

- **Special Orthogonal Group SO(k)**: Rotation matrices with geodesic distance metrics
- **SPD Manifolds**: Riemannian operations on positive definite matrices  
- **Graph Laplacians**: Spectral graph theory for connectivity analysis
- **Stiefel Manifolds**: Orthogonal matrix optimization via Procrustes methods

Numerical stability is ensured through:
- Eigenvalue tolerance checking
- SVD-based orthogonalization
- Condition number monitoring
- Robust matrix logarithms

## Citations

When using HATSA in your research, please cite:

```
Buchsbaum, B.R. (2024). HATSA: Hyperalignment via Task-informed Shared Analysis. 
R package version 0.1.0. https://github.com/bbuchsbaum/hatsa
```

## License

MIT License - See LICENSE file for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines. 