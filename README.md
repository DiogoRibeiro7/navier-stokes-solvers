# Navier–Stokes Solvers

[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-latest-brightgreen.svg)](TESTING_DOCUMENTATION.md)
[![CI](https://img.shields.io/github/actions/workflow/status/diogoribeiro7/navier-stokes-solvers/test-and-document.yml?branch=main&label=CI)](https://github.com/diogoribeiro7/navier-stokes-solvers/actions/workflows/test-and-document.yml)
[![Coverage](https://img.shields.io/badge/coverage-pending-lightgrey.svg)](docs/test_badges.md)
[![Tests](https://img.shields.io/badge/tests-automated-blue.svg)](TESTING_DOCUMENTATION.md)

High-order implementations of the two-dimensional incompressible Navier–Stokes equations featuring:

- An implicit Newton–Raphson finite difference solver with multigrid pressure correction.
- A Fourier pseudospectral solver with hyperviscosity, stochastic forcing, and adaptive time stepping.
- Reproducible configuration files, automated test reporting, and publication-ready artefacts (manuscript, slides, figures).

---

## What’s New

- **Continuous integration** via GitHub Actions (`test-and-document.yml`) covering Ubuntu/macos and GCC/Clang versions, Valgrind checks, and coverage.
- **Automated testing report**: the workflow publishes `TESTING_DOCUMENTATION.md` with pass/fail counts, trends, and platform metadata.
- **Configuration parser** (`ns_config`) supporting layered JSON/INI overrides with schema validation.
- **Benchmark toolkit** (`scripts/run_benchmarks.sh`) and figure scripts under `paper/figures/scripts/` to regenerate publication graphics.
- **Docker development environment** (`docker/Dockerfile.dev`, `docker/docker-compose-dev.yml`) with compilers, profilers, and Jupyter Lab.
- **Release automation** (`scripts/make_release.sh`, `CHANGELOG.md`, `VERSION`) aligned with semantic versioning.

---

## Repository Layout

```
.
├── config/                  # Base solver parameters (INI/JSON overrides)
├── docker/                  # Development container definitions
├── paper/                   # Manuscript, figures, captions, references
├── presentation/            # Slide deck sources and generator script
├── scripts/                 # Automation (install, validate, benchmarks, release)
├── src/                     # Finite difference & spectral solver sources
├── tests/                   # C regression tests + config fixtures
├── TESTING_DOCUMENTATION.md # Generated CI test report
└── docs/                    # Test history, badges, supplementary notes
```

---

## Getting Started

### 1. Install Dependencies

```bash
git clone https://github.com/diogoribeiro7/navier-stokes-solvers.git
cd navier-stokes-solvers
./scripts/install_dependencies.sh
```

The script detects Ubuntu/macOS, installs compilers, FFTW, Python tools, and verifies the toolchain.

### 2. Build Solvers

```bash
make spectral   # Fourier pseudospectral solver
make fd         # Finite difference solver
```

### 3. Quick Validation

```bash
./scripts/validate_installation.sh
```

Runs smoke builds and validation cases, reporting compiler, FFTW, and Python status.

### 4. Optional: Docker Development

```bash
docker compose -f docker/docker-compose-dev.yml up --build
```

Provides compiler, profiler, and Jupyter services with source volumes mounted at `/workspace`.

---

## Configuration

Runtime settings are driven by human-readable configuration files loaded by `ns_config`:

- `config/parameters.ini`: baseline solver parameters.
- Additional `.ini`/`.json` files can be passed via `--config` to override values.
- Validation schemas enforce consistency (positive Reynolds numbers, admissible forcing types, etc.).

Example:

```bash
bin/ns_spectral_solver \
  --config config/parameters.ini \
  --config overrides/high_re.json
```

---

## Testing & Continuous Integration

- Comprehensive test suites live in `tests/` (`test_solvers`, `test_conservation`, `test_taylor_green`).
- GitHub Actions matrix covers GCC/Clang (9–15) on Ubuntu 20.04/22.04 and macOS.
- Valgrind leak checks and optional coverage (lcov) run on Linux.
- Results, logs, and coverage reports are archived and summarised in **[TESTING_DOCUMENTATION.md](TESTING_DOCUMENTATION.md)**.
- Historical trends and badge snippets are stored under `docs/`.

To run locally:

```bash
make test          # C test suite
pytest tests/python  # If Python tests are added later
./scripts/run_benchmarks.sh  # Multi-resolution benchmark sweep
```

---

## Publication Artefacts

- **Manuscript**: `paper/manuscript.tex` (SIAM style) + `paper/references.bib`.
- **Figures**: reproducible via `paper/figures/scripts/*.py` (exports to `paper/figures/generated/` at 300 DPI).
- **Captions**: stored in `paper/figures/captions/` for editorial workflows.
- **Slides**: `presentation/slides.tex` and `presentation/generate_slides.py` produce conference-ready decks.
- **Blog draft**: `blog_post.md` summarises motivations, features, and visuals.

---

## Citation

```bibtex
@software{ribeiro2025navierstokes,
  author       = {Ribeiro, Diogo},
  title        = {High-Order Navier-Stokes Solvers: Newton-Raphson Finite
                  Difference and Fourier Spectral Methods},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://github.com/diogoribeiro7/navier-stokes-solvers}
}
```

Refer to `CITATION.cff` for machine-readable metadata and `docs/ZENODO_GUIDE.md` for deposit instructions.

---

## Contributing

We welcome issues, feature proposals, and pull requests! Please review:

- [`CONTRIBUTING.md`](CONTRIBUTING.md) for workflow, testing, and style.
- `.github/ISSUE_TEMPLATE/*` for bug/feature report templates.
- `.github/PULL_REQUEST_TEMPLATE.md` for PR checklist.

Key expectations:

- Include regression tests for new functionality.
- Update documentation/figures when behaviour changes.
- Regenerate `TESTING_DOCUMENTATION.md` via CI (auto-handled on merge).

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Acknowledgements

- **Diogo Ribeiro** (ESMAD – Instituto Politécnico do Porto) and community contributors.
- Reference implementations build upon classical benchmarks by Ghia et al. (1982), Taylor & Green (1937), and spectral methodology by Canuto et al. (2007).

—

_If these tools accelerate your research, please star the repository, cite the software, or share your results with the community!_
