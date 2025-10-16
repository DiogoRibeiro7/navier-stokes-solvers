---
layout: single
title: "API Reference"
permalink: /api-reference/
toc: true
---

# API Reference

This reference documents the public C API exposed by the finite-difference and spectral solvers. Types referenced below originate from [`include/common_types.h`](../include/common_types.h).

- See [Usage Guide](/usage-guide/) for practical walkthroughs.
- Mathematical notation is discussed in [Mathematical Formulation](mathematical_formulation.md).

## Common Types

### `SolverStatus`

Enum indicating the outcome of a solver step:

| Value | Meaning |
|-------|---------|
| `SOLVER_SUCCESS` | Step converged normally. |
| `SOLVER_MAX_ITER` | Newton iteration hit the maximum iteration count. |
| `SOLVER_DIVERGED` | Residual grew or stagnated. |
| `SOLVER_MEM_ERROR` | Memory allocation failed. |
| `SOLVER_INVALID_PARAM` | Illegal parameter detected (e.g., negative grid size). |

### `SolverStats`

```c
typedef struct {
    double energy;
    double enstrophy;
    double max_velocity;
    double max_divergence;
    double residual;
    int    iterations;
    double cpu_time;
} SolverStats;
```

Used to report diagnostics to calling code.

---

## Finite-Difference Module (`ns_fd_solver.h`)

### `NSFiniteDiffData* ns_fd_allocate(int nx, int ny, double L, double H, double Re);`
- **Description:** Allocates and initialises a finite-difference solver context for a rectangular domain.
- **Parameters:**
  - `nx`, `ny`: grid points in $x$ and $y$ (must be $\ge 3$).
  - `L`, `H`: domain lengths.
  - `Re`: Reynolds number (positive).
- **Returns:** Pointer to an owning `NSFiniteDiffData` or `NULL` on allocation failure.
- **Notes:** Arrays are zeroed. Default boundary conditions are no-slip with unit lid velocity.

### `void ns_fd_free(NSFiniteDiffData *data);`
Releases all memory associated with `data`. Safe to call with `NULL`.

### Initialisation Utilities

- `void ns_fd_initialize_lid_cavity(NSFiniteDiffData *data);`
- `void ns_fd_initialize_taylor_green(NSFiniteDiffData *data);`
- `void ns_fd_apply_boundary_conditions(NSFiniteDiffData *data);`

These routines set initial fields and impose boundary values. Call immediately after allocation or before restarting a run. `ns_fd_apply_boundary_conditions` should be invoked whenever ghost cells or boundary nodes need updating (e.g., after Newton updates).

### Newton and Linear Algebra

- `SolverStatus ns_fd_newton_raphson_step(NSFiniteDiffData *data);`
- `double ns_fd_compute_residual(NSFiniteDiffData *data);`
- `void ns_fd_assemble_jacobian(NSFiniteDiffData *data);`
- `void ns_fd_solve_linear_system(NSFiniteDiffData *data, double *delta);`

`ns_fd_newton_raphson_step` performs one Newton iteration internally calling the latter helpers. Use `ns_fd_compute_residual` to monitor convergence between iterations.

### Time Advancement

- `SolverStatus ns_fd_advance_timestep(NSFiniteDiffData *data);`
- `double ns_fd_adaptive_timestep(NSFiniteDiffData *data);`

`ns_fd_advance_timestep` wraps Newton solve, updates `data->t`, `data->dt`, and preserves previous states for convergence checks. `ns_fd_adaptive_timestep` implements the CFL/viscous heuristic; override `data->dt` if necessary before calling `ns_fd_advance_timestep`.

### Diagnostics and Output

- `void ns_fd_compute_stats(NSFiniteDiffData *data, SolverStats *stats);`
- `void ns_fd_analyze_convergence(NSFiniteDiffData *data, int step);`
- `void ns_fd_output_solution(NSFiniteDiffData *data, const char *filename);`

`ns_fd_output_solution` writes columns `x y u v p` to a text file. Example usage:

```c
NSFiniteDiffData *sim = ns_fd_allocate(128, 128, 1.0, 1.0, 1000.0);
ns_fd_initialize_lid_cavity(sim);
for (int step = 0; step < 5000; ++step) {
    if (ns_fd_advance_timestep(sim) != SOLVER_SUCCESS) break;
    if (step % 100 == 0) {
        SolverStats stats;
        ns_fd_compute_stats(sim, &stats);
        printf("t=%g, residual=%e\n", sim->t, stats.residual);
    }
}
ns_fd_output_solution(sim, "results/fd_snapshot.dat");
ns_fd_free(sim);
```

---

## Spectral Module (`ns_spectral_solver.h`)

### `NSSpectralData* ns_spectral_allocate(int nx, int ny, double Lx, double Ly, double Re);`
Allocates spectral solver structures. `nx` should be even; `ny` must match FFT plan requirements.

### `void ns_spectral_free(NSSpectralData *data);`
Destroys FFT plans, frees arrays, and calls `fftw_cleanup()`.

### Grid and Transforms

- `void ns_spectral_initialize_grids(NSSpectralData *data);`
- `void ns_spectral_initialize_fft_plans(NSSpectralData *data);`
- `void ns_spectral_initialize_taylor_green(NSSpectralData *data);`
- `void ns_spectral_physical_to_spectral(NSSpectralData *data, double *field, fftw_complex *field_hat);`
- `void ns_spectral_spectral_to_physical(NSSpectralData *data, fftw_complex *field_hat, double *field);`
- `void ns_spectral_apply_dealiasing(NSSpectralData *data, fftw_complex *field_hat);`

The transform utilities expect arrays of size `nx*ny` (real) or `nkx*nky` (complex). Always call `ns_spectral_initialize_grids` and `ns_spectral_initialize_fft_plans` once after allocation.

### Core Operators

- `void ns_spectral_compute_derivatives(NSSpectralData *data);`
- `void ns_spectral_compute_velocity_from_vorticity(NSSpectralData *data);`
- `void ns_spectral_compute_nonlinear_terms(NSSpectralData *data);`
- `void ns_spectral_compute_velocity_from_vorticity(NSSpectralData *data);`

These routines populate `u`, `v`, `omega`, and `nonlinear_hat` in-place. `ns_spectral_compute_derivatives` internally updates velocity and vorticity fields consistent with the latest spectral coefficients.

### Time Integration

- `void ns_spectral_advance_rk4(NSSpectralData *data);`
- `double ns_spectral_compute_cfl_timestep(NSSpectralData *data);`

`ns_spectral_advance_rk4` expects `data->dt` to be set and advances `omega_hat`, updating `t`. Call `ns_spectral_compute_cfl_timestep` beforehand to choose a stable `dt`.

### Analysis and Output

- `double ns_spectral_compute_energy(NSSpectralData *data);`
- `double ns_spectral_compute_enstrophy(NSSpectralData *data);`
- `void ns_spectral_analyze_convergence(NSSpectralData *data, int step);`
- `void ns_spectral_analyze_spectrum(NSSpectralData *data);`
- `int ns_spectral_check_resolution(NSSpectralData *data);`
- `void ns_spectral_output_solution(NSSpectralData *data, const char *filename);`

Example usage:

```c
NSSpectralData *spec = ns_spectral_allocate(256, 256, 2*M_PI, 2*M_PI, 4000.0);
ns_spectral_initialize_grids(spec);
ns_spectral_initialize_fft_plans(spec);
ns_spectral_initialize_taylor_green(spec);
for (int step = 0; step < 2000; ++step) {
    spec->dt = ns_spectral_compute_cfl_timestep(spec);
    ns_spectral_advance_rk4(spec);
    if (step % 100 == 0) {
        printf("t=%g, E=%g\n", spec->t, ns_spectral_compute_energy(spec));
    }
}
ns_spectral_output_solution(spec, "results/spectral_final_solution.dat");
ns_spectral_free(spec);
```

---

## Time Integration Utilities (`time_integration.h`)

These helper functions operate on scalar values and can be reused by external modules.

### `double compute_cfl_condition(double max_u, double max_v, double dx, double dy, double cfl_factor);`
Computes a stability-limited timestep $\Delta t$ based on advective speeds.

### `double compute_viscous_timestep(double dx, double dy, double Re);`
Returns the viscous limit $0.5\,\mathrm{Re}\,\Delta x\Delta y/(\Delta x + \Delta y)$.

### `double adaptive_timestep(double max_vel, double dx, double dy, double Re, AdaptiveTimeParams *params);`
Combines advective and viscous limits with user-defined safety factors (see struct `AdaptiveTimeParams`).

**Example:**

```c
AdaptiveTimeParams params = {.cfl_factor = 0.5, .dt_min = 1e-4, .dt_max = 1e-1,
                             .safety_factor = 0.9, .adapt_enabled = 1};
double dt = adaptive_timestep(max_vel, dx, dy, Re, &params);
```

---

## Notes & Caveats

- Functions returning pointers transfer ownership to the caller; invoke the corresponding `*_free` routine to avoid leaks.
- File output routines overwrite existing files. Use unique filenames per snapshot if accumulation is desired.
- Spectral routines require FFTW initialisation; link with `-lfftw3` and call `fftw_cleanup()` only once the simulation is finished.
- Jacobian assembly assumes structured Cartesian grids; extending to non-uniform meshes requires re-deriving difference operators.

For further reading, consult [Benchmark Results](/benchmark-results/) to see how these APIs were used to produce validation data.
