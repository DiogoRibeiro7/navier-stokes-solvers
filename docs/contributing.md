---
layout: single
title: "Contributing"
permalink: /contributing/
toc: true
---

# Contributing Guide

We welcome pull requests that improve the solvers, add physics models, or expand documentation. Please read this guide before submitting changes.

## Getting Started

1. Fork the repository on GitHub.
2. Clone your fork:
   ```bash
   git clone https://github.com/<your-username>/navier-stokes-solvers.git
   cd navier-stokes-solvers
   ```
3. Configure Git remotes:
   ```bash
   git remote add upstream https://github.com/original/navier-stokes-solvers.git
   ```

## Branching and Workflow

- Create feature branches from `main`:
  ```bash
  git checkout -b feature/my-improvement
  ```
- Keep branches focused; one logical change per pull request.
- Rebase on top of `main` before opening a PR:
  ```bash
  git fetch upstream
  git rebase upstream/main
  ```
- Use descriptive commit messages (imperative mood).

## Coding Standards

- **C code:** follow `clang-format` defaults with 4-space indentation, limit lines to 100 characters, and keep header guards consistent with existing style.
- **Python scripts:** conform to PEP 8. Run `black` (see `scripts/requirements.txt`) before committing.
- **Documentation:** Markdown with GitHub-flavoured extensions; cross-link related pages (see [docs/index](/usage-guide/)).

Refer to `include/common_types.h` for naming conventions (snake_case for functions, PascalCase for structs).

## Adding New Features

1. **Plan:** open a GitHub issue describing the feature and design approach.
2. **Implement:** extend the solver modules or add new ones under `src/`.
3. **Test:** add regression tests in `tests/` and update `tests/Makefile` if new binaries are required.
4. **Document:** update relevant Markdown files (e.g., [Usage Guide](/usage-guide/), [API Reference](/api-reference/)).
5. **Submit PR:** link to the original issue, summarise changes, and note any follow-up work.

## Testing Requirements

- Run `make test` (or the relevant subset) before submission; attach output in the PR.
- For numerical features, include verification scripts or benchmark comparisons (see [Benchmark Results](/benchmark-results/)).
- Ensure new scripts are covered by basic smoke tests (e.g., `python -m compileall`).

## Documentation Standards

- Every new feature must be described in the usage or theory docs.
- Add figures to `docs/images/` and reference them with `![caption](images/<file>.png)`.
- Keep language concise and accessible to graduate-level readers.

## Code Review Checklist

Before requesting review ensure:

- [ ] All tests pass locally.
- [ ] Linting/formatting applied.
- [ ] Documentation updated.
- [ ] No unrelated changes in the diff.

## Contact

For major design discussions open a GitHub discussion or email the maintainers listed in `README.md`. We aim to respond to PRs within five working days.
