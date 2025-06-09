# Contributing to HATSA

Thank you for your interest in contributing to HATSA! This document provides guidelines for contributing to the project.

## 🚀 Getting Started

### Prerequisites

- R (≥ 4.0.0)
- Git
- GitHub account

### Development Setup

1. **Fork and clone the repository**
   ```bash
   git clone https://github.com/your-username/hatsa.git
   cd hatsa
   ```

2. **Install development dependencies**
   ```r
   # Install package development tools
   install.packages(c("devtools", "roxygen2", "testthat", "styler", "covr"))
   
   # Install package dependencies
   devtools::install_deps(dependencies = TRUE)
   ```

3. **Load the package for development**
   ```r
   devtools::load_all()
   ```

## 🛠️ Development Workflow

### Code Standards

- **Style**: We use the [tidyverse style guide](https://style.tidyverse.org/)
- **Formatting**: Run `styler::style_pkg()` before committing
- **Documentation**: All functions must have complete roxygen2 documentation
- **Testing**: New code should include comprehensive tests

### Making Changes

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Write your code**
   - Follow existing code patterns
   - Add comprehensive documentation
   - Include appropriate tests

3. **Test your changes**
   ```r
   # Run tests
   devtools::test()
   
   # Check package
   devtools::check()
   
   # Check test coverage
   covr::package_coverage()
   ```

4. **Style your code**
   ```r
   styler::style_pkg()
   ```

### Documentation

- **Functions**: All exported functions need `@param`, `@return`, `@examples`, and `@export`
- **Classes**: S3 classes need clear documentation of structure and methods
- **Vignettes**: Complex features should include vignettes with examples

#### Example Documentation Template
```r
#' Brief function description
#'
#' Longer description with details about what the function does.
#'
#' @param param1 Description of parameter 1
#' @param param2 Description of parameter 2
#' @return Description of what the function returns
#' @export
#' @examples
#' # Example usage
#' result <- my_function(param1 = "value", param2 = 123)
```

### Testing

- **Coverage**: Aim for >90% test coverage
- **Test types**: Unit tests, integration tests, and edge cases
- **Test location**: Place tests in `tests/testthat/test-*.R`

#### Test Structure
```r
test_that("function does what it should", {
  # Setup
  input <- create_test_data()
  
  # Execute
  result <- my_function(input)
  
  # Verify
  expect_equal(result$expected_field, expected_value)
  expect_true(is.matrix(result$matrix_field))
})
```

## 📝 Pull Request Process

### Before Submitting

1. **Ensure all checks pass**
   ```r
   devtools::check()
   ```

2. **Run the full test suite**
   ```r
   devtools::test()
   ```

3. **Check code coverage**
   ```r
   covr::package_coverage()
   ```

4. **Update documentation**
   ```r
   devtools::document()
   ```

### PR Guidelines

- **Title**: Use a clear, descriptive title
- **Description**: Explain what changes you made and why
- **Testing**: Describe how you tested your changes
- **Breaking changes**: Clearly mark any breaking changes

### PR Template
```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

## Testing
- [ ] Tests pass locally
- [ ] Added tests for new functionality
- [ ] Updated documentation

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] Tests added/updated
```

## 🐛 Reporting Issues

### Bug Reports

When reporting bugs, include:

- **R version**: Output of `sessionInfo()`
- **HATSA version**: `packageVersion("hatsa")`
- **Minimal reproducible example**
- **Expected vs actual behavior**
- **Error messages** (full traceback)

### Feature Requests

For feature requests, provide:

- **Use case**: Why is this feature needed?
- **Proposed solution**: How should it work?
- **Alternatives**: What other approaches have you considered?

## 🏗️ Project Structure

```
hatsa/
├── R/                          # Source code
│   ├── hatsa_core_algorithm.R  # Main HATSA workflow
│   ├── task_hatsa_main.R       # Task-informed extensions
│   ├── voxel_projection.R      # Nyström voxel mapping
│   └── ...                     # Other modules
├── tests/testthat/             # Test files
├── man/                        # Generated documentation
├── vignettes/                  # Package vignettes
├── .github/workflows/          # CI/CD workflows
├── DESCRIPTION                 # Package metadata
└── NAMESPACE                   # Exported functions
```

## 📚 Resources

### Learning Resources

- [R Packages (2nd ed)](https://r-pkgs.org/) - Comprehensive guide to R package development
- [Advanced R](https://adv-r.hadley.nz/) - Deep dive into R programming
- [testthat documentation](https://testthat.r-lib.org/) - Testing framework

### HATSA-Specific Resources

- **Mathematical background**: See `docs/` folder for algorithm specifications
- **API documentation**: https://bbuchsbaum.github.io/hatsa/
- **Project status**: [HATSA_PROJECT_STATUS.md](./HATSA_PROJECT_STATUS.md)

## 🤝 Code of Conduct

### Our Standards

- **Respectful**: Be respectful and constructive in discussions
- **Inclusive**: Welcome contributors from all backgrounds
- **Collaborative**: Help others learn and improve
- **Professional**: Maintain professional conduct in all interactions

### Reporting Issues

If you experience or witness inappropriate behavior, please report it to [brad.buchsbaum@gmail.com](mailto:brad.buchsbaum@gmail.com).

## 🙏 Recognition

Contributors will be acknowledged in:

- Release notes
- Package documentation
- GitHub contributors list

Thank you for helping make HATSA better! 🎉 