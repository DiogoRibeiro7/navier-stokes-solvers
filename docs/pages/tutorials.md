---
layout: single
title: "Tutorials"
permalink: /tutorials/
toc: true
---

Hands-on tutorials guide you from first run to rigorous validation. Supplementary notebooks and datasets live in the `tutorials/` folder (planned). Video walk-throughs reference placeholder clips under `docs/assets/videos/` for future uploads.

## Tutorial 1: Your First Simulation

**Objective:** Run the finite-difference lid-driven cavity and visualise velocity profiles.

1. Build the solver: `make fd`
2. Execute a short run:
   ```bash
   ./bin/ns_fd_solver --steps 500 --output results/fd_quickstart.dat
   ```
3. Plot velocity magnitude and streamlines:
   ```bash
   python scripts/plot_solution.py results/fd_quickstart.dat \
       --output-prefix figures/fd_quickstart
   ```
4. Compare centreline profiles with literature using `scripts/compare_methods.py`.

<div class="video-callout">
  <span>Video:</span> <em>Coming soon – quickstart overview (placeholder).</em>
</div>

## Tutorial 2: Understanding Results

**Objective:** Interpret solver diagnostics, energy budgets, and convergence.

- Inspect `SolverStats` printed every 10 steps (`energy`, `residual`).
- Enable verbose logging by setting `VERBOSE=1` in the Makefile and rebuilding.
- Plot vorticity isolines and surface visualisations:
  ```bash
  python scripts/plot_solution.py results/fd_final_solution.dat \
      --vorticity-colormap coolwarm --surface-view 35 -60
  ```
- Read [Usage Guide](/usage-guide/) for file format details.

## Tutorial 3: Parameter Studies

**Objective:** Automate Reynolds number sweeps and analyse convergence.

1. Use the spectral driver for Taylor–Green vortex:
   ```bash
   ./bin/ns_spectral_solver --config configs/taylor_green_128.yml
   ```
2. Automate sweeps via `scripts/convergence_study.py`:
   ```bash
   python scripts/convergence_study.py --resolutions 64 96 128 192 \
       --fd-results-pattern "results/fd_re100_{nx:03d}.dat" \
       --output-dir studies/re100
   ```
3. Generate convergence plots and tables under `studies/re100`.
4. Capture summary figures for documentation (add to [Gallery](/gallery/)).

## Tutorial 4: Validation Workflow

**Objective:** Reproduce benchmark tables against Ghia et al. and spectral references.

- Run the high-resolution cavity example (`examples/lid_driven_cavity_re1000.c`).
- Export centreline comparisons:
  ```bash
  python scripts/compare_methods.py results/fd_final_solution.dat \
      results/spectral_final_solution.dat -o figures/fd_vs_spec
  ```
- Populate entries in [Benchmark Results](/benchmark-results/).
- Submit your findings via a pull request following the [Contributing Guide](/contributing/).

<div class="video-callout">
  <span>Video:</span> <em>Validation deep dive (placeholder link).</em>
</div>

## Additional Resources

- **Notebook templates:** (coming soon) interactive parameter sweeps.
- **Community sessions:** Monthly office hours announced on the project wiki.
