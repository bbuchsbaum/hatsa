# HATSA: Hyperalignment via Task-informed Shared Analysis

**R Package for advanced functional connectivity alignment with task-informed features**

[![Status](https://img.shields.io/badge/Status-75%25%20Complete-orange)](./HATSA_PROJECT_STATUS.md)
[![Critical Path](https://img.shields.io/badge/Critical%20Path-5%20Functions-red)](./HATSA_PROJECT_STATUS.md#critical-path-to-completion)

---

## 🎯 **Current Status**

The HATSA package implements advanced hyperalignment methods for functional neuroimaging data, with **major architectural components complete** but **5 critical core functions missing**.

**📊 Progress**: ~75% complete  
**📋 Current Focus**: [5 Critical Path tickets](./HATSA_PROJECT_STATUS.md#immediate-priorities-week-1)  
**⏱️ Estimated Completion**: 1-2 weeks  

---

## 🚀 **Quick Start (Once Complete)**

```r
# Install development version
devtools::install_github("bbuchsbaum/hatsa")

# Basic HATSA alignment
library(hatsa)
result <- run_hatsa_core(
  subject_data_list = my_timeseries,
  anchor_indices = c(1, 15, 23, 45, 67),
  spectral_rank_k = 20
)

# Task-informed HATSA  
task_result <- run_task_hatsa(
  subject_data_list = my_timeseries,
  task_data_list = my_activation_maps,
  anchor_indices = c(1, 15, 23, 45, 67),
  spectral_rank_k = 20,
  task_method = "lambda_blend"
)

# Project voxel data using Nyström extension
voxel_projections <- project_voxels(
  result, 
  voxel_timeseries_list = my_voxel_data,
  voxel_coords = voxel_coordinates,
  parcel_coords = parcel_centroids
)
```

---

## ✅ **What's Implemented**

### **Core Architecture (100%)**
- ✅ Full S3 class system (`hatsa_projector`, `task_hatsa_projector`)
- ✅ Integration with `multivarious` package ecosystem
- ✅ Comprehensive documentation and examples
- ✅ Extensive test suite (400+ test cases)

### **Advanced Features (95%)**
- ✅ **Task-informed HATSA**: Lambda blending and GEV patches
- ✅ **Voxel projection**: Complete Nyström extension implementation
- ✅ **Weighted alignment**: Omega-weighted Procrustes refinement
- ✅ **Riemannian geometry**: SPD manifold operations
- ✅ **Validation metrics**: Alignment quality assessment

### **Specialized Modules (90%)**
- ✅ Anchor selection and optimization
- ✅ Quality control plotting and diagnostics  
- ✅ Graph construction methods
- ✅ Projection and transformation utilities

---

## ❌ **What's Missing (Critical Path)**

**5 core functions** that everything depends on:

1. `compute_subject_connectivity_graph_sparse()` - Sparse correlation graphs
2. `compute_graph_laplacian_sparse()` - Laplacian computation  
3. `compute_spectral_sketch_sparse()` - Eigendecomposition
4. `misalign_deg()` - SO(k) geodesic distance
5. `solve_gev_laplacian_primme()` - Generalized eigenvalue solver

**📋 Detailed tickets**: See [HATSA_PROJECT_STATUS.md](./HATSA_PROJECT_STATUS.md)

---

## 📖 **Documentation**

- **[Project Status & Roadmap](./HATSA_PROJECT_STATUS.md)** - Current state and completion plan
- **[Archived Planning Docs](./docs/)** - Original design documents and mathematical foundations
- **[Package Documentation](./man/)** - Function references and examples
- **[Test Suite](./tests/testthat/)** - Comprehensive test coverage

---

## 🔧 **Development Setup**

```bash
# Clone and setup
git clone https://github.com/bbuchsbaum/hatsa.git
cd hatsa

# Install dependencies
R -e "install.packages(c('Matrix', 'RSpectra', 'PRIMME', 'expm', 'multivarious'))"

# Test current state (will show missing functions)
R -e "devtools::test()"

# Check package structure
R -e "devtools::check()"
```

### **Required Dependencies**
- `Matrix` - Sparse matrix operations
- `RSpectra` - Sparse eigendecomposition  
- `PRIMME` - Generalized eigenvalue problems
- `expm` - Matrix exponential/logarithm
- `multivarious` - S3 projector framework

---

## 🎯 **Next Steps for Contributors**

### **Immediate Priority (Week 1)**
Implement the 5 critical functions following the detailed specifications in [HATSA_PROJECT_STATUS.md](./HATSA_PROJECT_STATUS.md):

1. **Day 1-2**: Connectivity graphs and Laplacian computation
2. **Day 3-4**: Spectral sketches and rotation metrics  
3. **Day 5**: GEV solver

### **Integration Testing (Week 2)**
1. Validate full HATSA pipeline
2. Test task-informed extensions
3. Verify voxel projection functionality
4. Performance optimization

### **Ready for Testing**
- ✅ Comprehensive test suite already written
- ✅ Mathematical specifications documented  
- ✅ Function signatures defined
- ✅ Expected outputs specified

---

## 📊 **Package Statistics**

| Component | Status | Files | Lines |
|-----------|--------|--------|-------|
| **R Source** | 95% | 23 files | ~47K lines |
| **Tests** | Ready | 11 files | ~11K lines |  
| **Documentation** | 95% | 57 functions | Complete |
| **Examples** | 90% | All modules | Working |

---

## 🏗️ **Architecture Overview**

```
hatsa/
├── R/
│   ├── hatsa_core_algorithm.R      # Main HATSA workflow
│   ├── task_hatsa_main.R           # Task-informed extensions  
│   ├── voxel_projection.R          # Nyström voxel mapping
│   ├── spectral_graph_construction.R  # ❌ Missing core functions
│   ├── riemannian_geometry.R       # SPD manifold operations
│   └── [18 other modules]          # Specialized functionality
├── tests/testthat/
│   ├── test-hatsa_core_functionality.R
│   ├── test-spectral_graph_construction.R  # Ready to test
│   └── [9 other test files]
└── docs/                           # Archived planning documents
```

---

## 🤝 **Contributing**

1. **Focus on Critical Path**: Implement the 5 missing core functions first
2. **Follow Test-Driven Development**: Comprehensive tests already exist
3. **Reference Documentation**: Mathematical specs in `docs/` folder
4. **Check Integration**: Ensure functions work with existing S3 methods

---

## 📜 **License**

MIT License - See LICENSE file for details.

---

**🎯 Ready to complete this project? Start with [HATSA_PROJECT_STATUS.md](./HATSA_PROJECT_STATUS.md) for detailed implementation tickets.** 