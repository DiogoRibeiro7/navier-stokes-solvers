---
layout: single
title: "Theory Background"
permalink: /theory-background/
toc: true
---

# Theory Background

This primer targets advanced undergraduate students seeking intuition for the phenomena simulated by the code. For derivations consult [Mathematical Formulation](mathematical_formulation.md).

## Fluid Dynamics Essentials

- **Continuum assumption:** Fluids are modelled as continuous media when molecular mean free paths are much smaller than the domain size.
- **Mass conservation:** Enforced through the incompressibility constraint $\nabla \cdot \boldsymbol{u} = 0$.
- **Momentum balance:** Newton's second law applied to a fluid element yields the Navier–Stokes equations.

ASCII schematic of a control volume:

```
+-----------+
|           |  u(x+dx,y)
|   CV      |<-----------
|           |
+-----------+
^           ^
u(x,y)    
```

See `docs/images/grid_structure_diagram.txt` for an illustration slot.

## Reynolds Number

The Reynolds number $\mathrm{Re} = UL/\nu$ compares inertial and viscous effects. Flow regimes:

- $\mathrm{Re} \lesssim 10^2$: laminar, viscous forces dominate.
- $\mathrm{Re} \sim 10^3$–$10^4$: transition and coherent vortices appear (e.g., lid-driven cavity at 1000).
- $\mathrm{Re} \gg 10^4$: turbulence with wide range of scales.

## Vorticity Dynamics

Vorticity $\omega = \nabla \times \boldsymbol{u}$ measures local rotation. In 2D incompressible flow the vorticity equation simplifies to

$$
\frac{\partial \omega}{\partial t} + \boldsymbol{u}\cdot\nabla \omega = \nu \nabla^2 \omega.
\tag{1}
$$

Key processes:

- **Advection:** Vorticity is transported by the flow.
- **Diffusion:** Viscosity smooths gradients.
- **Boundary generation:** No-slip walls create vorticity layers that diffuse into the interior.

Figure placeholder: `docs/images/example_velocity_fields.txt` (add contour plots of $\omega$).

## Turbulence Basics

- **Energy cascade:** In 3D turbulence energy flows from large to small scales; in 2D, enstrophy cascades to small scales while energy moves to large scales (Kraichnan, 1967).
- **Kolmogorov spectrum:** $E(k) \sim k^{-5/3}$ in inertial range. Spectral solver diagnostics (see [Benchmark Results](benchmark_results.md#spectral-energy-spectra)) reproduce this slope when adequately resolved.
- **Coherent structures:** Kelvin–Helmholtz rollers, Taylor–Green vortices, and lid-driven corner eddies are examples simulated here.

## Energy Cascades

A qualitative depiction of the dual cascade in 2D:

```
Large eddies --> (inverse cascade) --> Larger coherent structures
          |
          +--> (forward enstrophy cascade) --> Small-scale vorticity filaments
```

Spectral energy spectra stored in `docs/images/spectral_basis_functions.txt` should visualise these trends.

## Suggested Reading

- White, F. M. (2006). *Viscous Fluid Flow* (3rd ed.). McGraw-Hill.
- Pope, S. B. (2000). *Turbulent Flows*. Cambridge University Press.
- Davidson, P. A. (2015). *Turbulence: An Introduction for Scientists and Engineers*. Oxford University Press.
- Batchelor, G. K. (1967). *An Introduction to Fluid Dynamics*. Cambridge University Press.

For implementation details revisit [Numerical Methods](/numerical-methods/).
