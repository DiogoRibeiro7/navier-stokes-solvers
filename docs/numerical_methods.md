---
layout: single
title: "Numerical Methods"
permalink: /numerical-methods/
toc: true
---

# Numerical Methods

> Complementary reading: [Mathematical Formulation](mathematical_formulation.md) and [Theory Background](/theory-background/).

## Finite-Difference Discretization

### Spatial Operators

We employ second-order centered differences on a uniform Cartesian grid $(x_i, y_j)$ with spacings $\Delta x, \Delta y$:

$$
\left.\frac{\partial u}{\partial x}\right|_{i,j} \approx \frac{u_{i+1,j} - u_{i-1,j}}{2\Delta x}, \qquad
\left.\frac{\partial^2 u}{\partial x^2}\right|_{i,j} \approx \frac{u_{i+1,j} - 2u_{i,j} + u_{i-1,j}}{\Delta x^2}.
\tag{1}
$$

Figure reference: see `docs/images/stencil_illustration.txt`.

### Temporal Discretization (Backward Euler)

The semi-discrete momentum equations are advanced implicitly via backward Euler:

$$
\frac{\boldsymbol{u}^{n+1} - \boldsymbol{u}^n}{\Delta t} + \mathcal{N}(\boldsymbol{u}^{n+1}) = -\nabla p^{n+1} + \frac{1}{\mathrm{Re}} \nabla^2 \boldsymbol{u}^{n+1},
\tag{2}
$$

where $\mathcal{N}$ denotes the nonlinear advection terms. The residual is solved by Newton-Raphson.

### Newton-Raphson Linearisation

Given an iterate $\boldsymbol{u}^{(k)}$, the correction $\delta \boldsymbol{u}$ solves

$$
\mathbf{J}(\boldsymbol{u}^{(k)})\, \delta \boldsymbol{u} = -\mathbf{R}(\boldsymbol{u}^{(k)}),
\tag{3}
$$

with sparse Jacobian $\mathbf{J}$ assembled in compressed-row format (`J`, `ia`, `ja`). The Newton iteration terminates when $\|\mathbf{R}\|_\infty < 10^{-8}$ or after `MAX_ITER` steps.

**Pseudocode:**

```pseudo
function advance_fd_step(state):
    copy_previous_state(state)
    state.dt = adaptive_timestep(state)
    for k in 0 .. MAX_ITER:
        residual = assemble_residual(state)
        if norm_inf(residual) < tol:
            return SUCCESS
        jacobian = assemble_jacobian(state)
        delta = solve_sparse_system(jacobian, residual)
        update_state(state, delta)
    return MAX_ITER_REACHED
```

Sparse linear systems are solved with a bi-conjugate gradient stabilized (BiCGSTAB) routine; see Ferziger & Perić (2002) for similar strategies.

## Spectral (Fourier) Methods

### Basis and Transforms

Spectral simulations assume periodic domains and expand vorticity in Fourier modes:

$$
\omega(x,y,t) = \sum_{k_x=-N_x/2}^{N_x/2} \sum_{k_y=-N_y/2}^{N_y/2} \hat{\omega}_{k_x,k_y}(t) e^{i(k_x x + k_y y)}.
\tag{4}
$$

Fast Fourier Transforms (FFTs) and their inverses are executed via FFTW3 plans (`fftw_plan_dft_r2c_2d`, `fftw_plan_dft_c2r_2d`). Reference: Frigo & Johnson (2005).

### Pseudospectral Evaluation

Nonlinear products are computed in physical space to avoid expensive convolutions:

1. Transform $\omega^n$ to spectral space.
2. Recover velocity via $\hat{\psi} = -\hat{\omega}/k^2$ and (3).
3. Transform back to physical space to form Jacobian $J(\psi,\omega)$.
4. Transform the product to spectral space for time integration.

The 2/3 dealiasing rule zeroes modal coefficients with $|k| > 2N/3$ to suppress aliasing (Orszag, 1971). See placeholder diagram `docs/images/spectral_basis_functions.txt`.

### Spectral Poisson Solver

In Fourier space the Poisson problem reduces to algebraic inversion:

$$
\hat{\psi}_{k_x,k_y} = -\frac{\hat{\omega}_{k_x,k_y}}{k_x^2 + k_y^2},
\tag{5}
$$

with zero mean mode handled separately ($\hat{\psi}_{0,0} = 0$).

## Time Integration

### Fourth-Order Runge–Kutta (RK4)

For vorticity transport we use classical RK4:

```pseudo
function rk4_step(omega, dt):
    k1 = rhs(omega)
    k2 = rhs(omega + 0.5 * dt * k1)
    k3 = rhs(omega + 0.5 * dt * k2)
    k4 = rhs(omega + dt * k3)
    return omega + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
```

Derivation follows Butcher (2008). The `rhs` operator internally performs FFTs and diffusion terms.

### Adaptive Timestep Control

Both solvers estimate the permissible timestep using a CFL-like condition:

$$
\Delta t = \mathrm{CFL}\times \min\left(\frac{\Delta x}{|u|_{\max}+\epsilon}, \frac{\Delta y}{|v|_{\max}+\epsilon}\right).
\tag{6}
$$

Viscous constraints and global bounds (`MIN_DT`, `MAX_DT`) further limit $\Delta t$ (see `ns_fd_adaptive_timestep` and `ns_spectral_compute_cfl_timestep`).

### Stability Considerations

- **FD backward Euler** is A-stable, so stiffness from viscous terms is handled implicitly, but Newton linearisation must remain convergent; a good initial guess from the previous time level is essential.
- **Spectral RK4** is conditionally stable; the CFL factor 0.3 used in `spectral_memory.c` keeps the advective Courant number below 1.

Von Neumann analysis for linearised advection-diffusion suggests stability if $\Delta t \le 2/(|\boldsymbol{u}|k + \nu k^2)$, consistent with (6).

## References

- Butcher, J. C. (2008). *Numerical Methods for Ordinary Differential Equations* (2nd ed.). Wiley.
- Canuto, C., Hussaini, M., Quarteroni, A., & Zang, T. (2006). *Spectral Methods: Fundamentals in Single Domains*. Springer.
- Ferziger, J. H., & Perić, M. (2002). *Computational Methods for Fluid Dynamics*. Springer.
- Frigo, M., & Johnson, S. G. (2005). The design and implementation of FFTW3. *Proceedings of the IEEE*, 93(2), 216–231.
- Orszag, S. A. (1971). On the elimination of aliasing in finite-difference schemes by filtering high-wavenumber components. *J. Atmos. Sci.*, 28, 1074–1074.
