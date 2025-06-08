# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

HATSA (Hyperalignment via Task-informed Shared Analysis) is an R package for advanced functional connectivity alignment in neuroimaging data. The package is ~75% complete with 5 critical core functions missing.

## Development Commands

```bash
# Install dependencies (including missing ones from DESCRIPTION)
R -e "install.packages(c('Matrix', 'RANN', 'stats', 'multivarious', 'expm', 'RSpectra', 'PRIMME', 'testthat', 'vegan'))"

# Run all tests
R -e "devtools::test()"

# Run specific test file
R -e "devtools::test_file('tests/testthat/test-spectral_graph_construction.R')"

# Check package (build, tests, examples, documentation)
R -e "devtools::check()"

# Build and install package locally
R -e "devtools::install()"

# Document package (rebuild .Rd files from roxygen2 comments)
R -e "devtools::document()"

# Load package for interactive development
R -e "devtools::load_all()"
```

## Architecture & Key Components

### S3 Class Hierarchy
- `hatsa_projector`: Main projector class, inherits from `multivarious::multiblock_biprojector`
- `task_hatsa_projector`: Task-informed variant with additional methods

### Core Algorithms
1. **HATSA Core** (`run_hatsa_core`): Base hyperalignment using sparse graph representations
2. **Task-HATSA** (`run_task_hatsa`): Extensions with lambda blending, GEV patches, anchor augmentation
3. **Voxel Projection**: Nyström extension for mapping full voxel data to spectral space

### Critical Missing Functions (Priority)
Located in `R/spectral_graph_construction.R`:
1. `compute_subject_connectivity_graph_sparse()` - Build sparse correlation graphs
2. `compute_graph_laplacian_sparse()` - Compute normalized Laplacian
3. `compute_spectral_sketch_sparse()` - Eigendecomposition with DC filtering
4. `misalign_deg()` - SO(k) geodesic distance metric
5. `solve_gev_laplacian_primme()` - Generalized eigenvalue solver

### Important Note on Dependencies
- `PRIMME`: Used for both sparse eigendecomposition and generalized eigenvalue problems (preferred over RSpectra for speed)
- All dependencies are now properly declared in DESCRIPTION

### Testing Structure
- Framework: testthat v3
- Test files mirror source structure
- 400+ tests already written (including for missing functions)
- Snapshot tests use `_snaps/` directory

### Important Design Patterns
1. All sparse matrices use Matrix::dgCMatrix format
2. Subject data format: list of T×V matrices (time × voxels)
3. Anchor indices are 1-based (R convention)
4. Methods follow S3 dispatch pattern from multivarious package
5. Extensive input validation in all public functions

## Mathematical & Implementation Details

### Task-HATSA Methods
1. **Lambda blending**: Combines connectivity and task graphs with W = (1-λ)W_conn + λW_task (default λ=0.15)
2. **GEV patches**: Uses generalized eigenvalue decomposition for task-orthogonal eigenvectors
3. **Anchor augmentation**: Appends task features to anchor matrix with differential weighting
4. **Auto-residualization**: When task/connectivity correlation > 0.45, task graph is residualized

### Graph Construction
- **Laplacian formula**: L_rw_lazy = I - α*D^(-1)*W (α=0.93 default)
- **k-NN sparsification**: Separate positive/negative edge selection
- **Z-scoring**: Applied to non-zero edge weights
- **DC filtering**: First eigenvector removed in spectral decomposition

### Riemannian Geometry Components
- **SPD manifold operations**: Log-Euclidean and AIRM metrics implemented
- **Mitteroecker & Bookstein distance**: Primary metric for HATSA output analysis
- **Geo-HATSA variant**: Uses Fréchet mean on SO(k) for geometric rotation updates
- **Tangent space projection**: Maps SPD matrices to common tangent space for analysis

### Validation & Diagnostics
- **Anchor selection metrics**: Reconstruction error, rotation dispersion, condition number
- **MRA-Select algorithm**: Balances condition number and Riemannian dispersion
- **Quality metrics**: ISC transferability, eigenvalue fidelity, graph correlation
- **Diagnostic plots**: k-stability, MDS visualization, eigengap analysis