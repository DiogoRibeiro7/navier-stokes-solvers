# Changelog

All notable changes to this project will be documented in this file, following [Semantic Versioning](https://semver.org/).

## v1.0.0 - 2025-01-16

### Added
- Newton–Raphson finite difference solver with configurable boundary conditions and multigrid pressure correction.
- Fourier pseudospectral solver featuring 4th-order Runge–Kutta integration, hyperviscosity, and stochastic forcing.
- Publication pack: SIAM-style manuscript, figures, captions, blog post, conference slides, and metadata (CITATION.cff, Zenodo JSON).
- Automated configuration system supporting layered JSON/INI overrides with schema validation.
- Benchmark generation scripts and performance plots for cavity and Taylor–Green validation.

### Changed
- Standardised repository layout under `paper/`, `scripts/`, `config/`, and `docker/`.
- Enhanced tests to cover configuration loading, multigrid Poisson solves, and hyperviscosity damping.

### Fixed
- Improved portability for figure-generation scripts (ASCII-only, mathtext-safe titles).

### Contributors
- Diogo Ribeiro (@diogoribeiro7)
- Open-source community contributors

> **Breaking changes:** Initial stable release (`v1.0.0`). Future changes will note incompatibilities here.
