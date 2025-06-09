# GitHub Actions CI/CD for HATSA

This directory contains GitHub Actions workflows for continuous integration and deployment of the HATSA R package.

## Workflows

### 1. R-CMD-check.yaml
**Purpose**: Comprehensive R package checking across multiple platforms
- **Triggers**: Push/PR to main/master branches
- **Platforms**: Ubuntu (devel, release, oldrel-1), Windows (release), macOS (release)
- **Actions**: 
  - Install R and dependencies
  - Run `R CMD check`
  - Upload test snapshots

### 2. test-coverage.yaml
**Purpose**: Code coverage measurement and reporting
- **Triggers**: Push/PR to main/master branches
- **Platform**: Ubuntu latest
- **Actions**:
  - Run `covr::package_coverage()`
  - Generate Cobertura XML report
  - Upload to Codecov

### 3. pkgdown.yaml
**Purpose**: Build and deploy package documentation website
- **Triggers**: Push to main/master, PRs, releases, manual dispatch
- **Platform**: Ubuntu latest
- **Actions**:
  - Build pkgdown site
  - Deploy to GitHub Pages (on push to main only)

### 4. style.yaml
**Purpose**: Automated code styling with styler
- **Triggers**: Push/PR to main/master branches
- **Platform**: Ubuntu latest
- **Actions**:
  - Run `styler::style_pkg()`
  - Commit and push style changes (on push only)

## Setup Requirements

### For Repository Owner
1. **Enable GitHub Pages**: Go to Settings > Pages > Source: GitHub Actions
2. **Codecov Token** (optional): Add `CODECOV_TOKEN` secret for private repos
3. **Branch Protection**: Consider requiring status checks to pass before merging

### For Contributors
- No additional setup required
- All workflows run automatically on PRs
- Style changes only applied to push events (not PRs)

## Configuration Files

- **codecov.yml**: Codecov configuration with coverage targets and ignore patterns
- **DESCRIPTION**: Updated with `covr` and `styler` in Suggests field

## Badges

The following badges are included in the main README:
- R-CMD-check status
- Codecov test coverage
- Links to GitHub Actions and Codecov dashboard

## Status

✅ All workflows are ready to run
✅ Cross-platform testing (Windows, macOS, Linux)  
✅ Multiple R versions (devel, release, oldrel-1)
✅ Code coverage reporting
✅ Automated documentation deployment
✅ Code style enforcement 