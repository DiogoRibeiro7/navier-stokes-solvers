# Contributing Guidelines

Thank you for your interest in improving the Navier–Stokes Solvers project!

## Code of Conduct

We follow the [Contributor Covenant](https://www.contributor-covenant.org/). Be respectful to all contributors. Report unacceptable behaviour to dfr@esmad.ipp.pt.

## Getting Started

1. **Fork** the repository and create a feature branch (`git checkout -b feature/my-improvement`).
2. **Install dependencies** via `./scripts/install_dependencies.sh`.
3. **Run validation** with `./scripts/validate_installation.sh`.
4. **Develop and test**:
   - Add or update regression tests in `tests/`.
   - Regenerate figures or benchmarks if outputs change.
5. **Lint and format**: use `clang-format`, `black`, and `flake8` where applicable.
6. **Submit a pull request** using the template provided; reference related issues and attach relevant logs/figures.

## Testing Checklist

- `make test`
- `./scripts/run_benchmarks.sh` (if performance changes)
- `./scripts/validate_installation.sh`
- Ensure `TESTING_DOCUMENTATION.md` is updated by CI.

## Reporting Issues

Open an issue using the templates in `.github/ISSUE_TEMPLATE/`. Provide reproduction steps, configurations, logs, and system details.

## Style Guidelines

- Prefer ANSI C11, limit dependencies in core solvers.
- Use descriptive commit messages (Imperative mood).
- Document public interfaces in headers and update relevant docs.

We appreciate contributions of all forms—code, documentation, testing, and triage. Thank you!
