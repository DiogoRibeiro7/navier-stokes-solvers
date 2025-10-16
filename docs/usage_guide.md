---
layout: single
title: "Usage Guide"
permalink: /usage-guide/
toc: true
---

# Usage Guide

> New to the codebase? Skim [Theory Background](/theory-background/) and [Numerical Methods](/numerical-methods/) for context before diving in.

## 5-Minute Quick Start

```bash
# 1. Clone and enter the project directory
git clone https://github.com/your-org/navier-stokes-solvers.git
cd navier-stokes-solvers

# 2. Build both solvers (requires gcc/clang, FFTW3)
make all

# 3. Run the finite-difference lid-driven cavity
./bin/ns_fd_solver

# 4. Visualise the result
python scripts/plot_solution.py results/fd_final_solution.dat -o figures/fd_run
```

Expected outputs are written under `results/`. For a richer demo, see `examples/` once populated.

## First Simulation (Step-by-Step)

1. **Configure parameters** in `src/finite_difference/main_fd.c` (grid size, Reynolds number, final time).
2. **Build the solver:** `make fd`.
3. **Execute:** `./bin/ns_fd_solver > logs/fd_run.txt 2>&1`.
4. **Check diagnostics:** inspect residual history printed every 10 steps.
5. **Plot key fields:**
   ```bash
   python scripts/plot_solution.py results/fd_final_solution.dat \
       --output-prefix figures/fd_final --contour-levels 40
   ```
6. **Compare methods:**
   ```bash
   python scripts/compare_methods.py results/fd_final_solution.dat \
       results/spectral_final_solution.dat -o figures/fd_vs_spec
   ```

Cross-link: deeper details in [Benchmark Results](/benchmark-results/).

## Parameter Selection Guidelines

- **Grid resolution:** start with $64 \times 64$ for laminar runs; increase to $192 \times 192$ or more for $\mathrm{Re} \ge 1000$.
- **Time step limits:** rely on adaptive controllers; tighten `MAX_DT` for highly unsteady flows.
- **Spectral runs:** maintain square grids with powers of two for FFT efficiency.
- **Physical scaling:** ensure nondimensional variables map to realistic physical units if comparing experiments (see [Mathematical Formulation](mathematical_formulation.md#non-dimensionalisation-and-reynolds-number)).

## Understanding Output Files

| File | Description |
|------|-------------|
| `results/fd_final_solution.dat` | Space-delimited table `x y u v [omega] [p]` depending on solver |
| `results/spectral_solution_step_XXXX.dat` | Snapshot at RK4 substep `XXXX` |
| `results/taylor_green_errors.dat` | Time series of L2/Linf errors |
| `figures/*` | PNG figures produced by plotting scripts |

Metadata is embedded via header comments (`# key: value`). Use `head -n 5 file.dat` to preview.

## Visualisation Workflow

1. **Static plots:** `scripts/plot_solution.py` (contours, surface plots, vorticity isolines).
2. **Animations:**
   ```bash
   python scripts/animate_flow.py "results/spectral_step_*.dat" \
       -o figures/flow.mp4 --fps 30 --particle-count 200
   ```
3. **Comparisons:** leverage `scripts/compare_methods.py` for FD vs spectral metrics.
4. **Spectral energy:** `python scripts/energy_spectrum.py results/spectral_final_solution.dat`.

Suggested image placeholders live under `docs/images/` so you can drop high-resolution renders for documentation.

## Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| Solver diverges (`SOLVER_DIVERGED`) | Reynolds number too high for given grid | Refine grid, reduce `MAX_DT`, or switch to spectral solver |
| `fftw_plan` failure | FFTW missing or misconfigured | Install FFTW3, re-run `make clean && make spectral` |
| Segmentation fault on start | Grid sizes not initialized | Check `ns_fd_allocate` return value; ensure memory succeeds |
| Flat velocity field | Forgot to set lid speed or initial condition | Revisit initialization routine |

More FAQs in [docs/faq.md](/faq/).

## Performance Optimisation

- **Compiler flags:** enable `-O3 -march=native` (edit `Makefile`).
- **Threaded FFTs:** configure FFTW with `--enable-openmp` and set `FFTW_NUM_THREADS`.
- **Batch snapshots:** lower `output_interval` to reduce I/O.
- **Matrix reuse:** for steady FD problems, reuse Jacobian structure to warm-start iterative solvers.

## Memory Estimates

For structured grids memory scales as $\mathcal{O}(N_x N_y)$. Approximate usage:

| Grid | FD (MiB) | Spectral (MiB) |
|------|----------|----------------|
| $64^2$ | ~5 | ~12 |
| $128^2$ | ~20 | ~42 |
| $256^2$ | ~80 | ~160 |

Estimates assume double precision arrays and FFT workspaces; adjust for additional diagnostics. Use `htop` or `Get-Process` to monitor consumption during runs.
