---
layout: single
title: "Mathematical Formulation"
permalink: /mathematical-formulation/
toc: true
---

# Mathematical Formulation

> See also: [Numerical Methods](/numerical-methods/) and [Theory Background](/theory-background/).

## Governing Equations

We solve the incompressible Navier-Stokes equations in two space dimensions for velocity $\boldsymbol{u} = (u, v)$ and pressure $p$:

$$
\begin{aligned}
\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} &= -\frac{1}{\rho}\frac{\partial p}{\partial x} + \nu\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right), \\
\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} &= -\frac{1}{\rho}\frac{\partial p}{\partial y} + \nu\left(\frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2}\right), \\
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} &= 0.
\end{aligned}
\tag{1}
$$

Here $\rho$ is the fluid density and $\nu$ the kinematic viscosity. Equation (1) matches the presentation in standard texts such as Panton (2013) and Ferziger & Perić (2002).

## Stream Function–Vorticity Formulation

For the spectral solver we adopt the vorticity $\omega$ and stream function $\psi$ variables:

$$
\omega = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}, \qquad \nabla^2 \psi = -\omega,
\tag{2}
$$

with velocity recoveries

$$
 u = \frac{\partial \psi}{\partial y}, \qquad v = -\frac{\partial \psi}{\partial x}. 
\tag{3}
$$

The transport equation for $\omega$ becomes

$$
\frac{\partial \omega}{\partial t} + J(\psi, \omega) = \nu \nabla^2 \omega,
\tag{4}
$$

where $J(\psi, \omega) = \frac{\partial \psi}{\partial x} \frac{\partial \omega}{\partial y} - \frac{\partial \psi}{\partial y} \frac{\partial \omega}{\partial x}$ is the Jacobian. This formulation eliminates pressure and is particularly convenient for Fourier-based solvers (Boyd, 2001; Canuto et al., 2006).

## Non-dimensionalisation and Reynolds Number

We choose characteristic length $L$, velocity $U$, and define nondimensional variables:

$$
\tilde{\boldsymbol{x}} = \frac{\boldsymbol{x}}{L}, \quad \tilde{t} = \frac{tU}{L}, \quad \tilde{\boldsymbol{u}} = \frac{\boldsymbol{u}}{U}, \quad \tilde{p} = \frac{p}{\rho U^2}.
\tag{5}
$$

Substituting (5) into (1) yields the dimensionless Navier-Stokes equations with a single parameter, the Reynolds number

$$
\mathrm{Re} = \frac{UL}{\nu},
\tag{6}
$$

measuring the ratio of inertial to viscous forces. The nondimensional system is

$$
\begin{aligned}
\frac{\partial \tilde{\boldsymbol{u}}}{\partial \tilde{t}} + (\tilde{\boldsymbol{u}} \cdot \tilde{\nabla})\tilde{\boldsymbol{u}} &= -\tilde{\nabla}\tilde{p} + \frac{1}{\mathrm{Re}} \tilde{\nabla}^2 \tilde{\boldsymbol{u}}, \\
\tilde{\nabla} \cdot \tilde{\boldsymbol{u}} &= 0.
\end{aligned}
\tag{7}
$$

All discretisations in this project use the nondimensional variables.

## Boundary Conditions

Multiple boundary-condition types are supported via `BoundaryCondition` enumerations:

- **No-slip ($\boldsymbol{u}=0$):** Used on solid walls in lid-driven cavities except for the moving lid.
- **Moving lid:** Imposes prescribed tangential velocity $u = U_{\text{lid}}$, $v = 0$.
- **Neumann pressure:** Zero normal pressure gradient $\partial p/\partial n = 0$ closes the pressure Poisson problem in FD runs.
- **Periodic:** Spectral simulations assume doubly-periodic domains, enforcing $\boldsymbol{u}(x+L,y) = \boldsymbol{u}(x,y)$ and similar in $y$.
- **Free-slip and mixed options** may be added by extending the boundary handling routines.

Boundary operators must be applied consistently to velocity, vorticity, and stream function (see [Numerical Methods](numerical_methods.md#finite-difference-discretization)).

## Initial Conditions

Typical initial fields include:

- **Rest state:** $\boldsymbol{u}=0$, $p=0$, used before spinning up lid-driven cavities.
- **Taylor–Green vortex:** Analytical velocity/vorticity given by
  $$u(x, y, 0) = \sin x\, \cos y, \quad v(x, y, 0) = -\cos x\, \sin y,$$
  $$\omega(x, y, 0) = 2 \cos x \cos y,$$
  scaled to desired wavelength and Reynolds number.
- **Shear layer / vortex merger:** Smooth hyperbolic tangent profiles with localized perturbations to trigger instability.

The initialization routines in `fd_initialization.c` and `spectral_initialization.c` implement these states and call `ns_fd_apply_boundary_conditions` or spectral transforms as appropriate.

## References

- Boyd, J. P. (2001). *Chebyshev and Fourier Spectral Methods*. Dover.
- Canuto, C., Hussaini, M., Quarteroni, A., & Zang, T. (2006). *Spectral Methods: Fundamentals in Single Domains*. Springer.
- Ferziger, J. H., & Perić, M. (2002). *Computational Methods for Fluid Dynamics*. Springer.
- Panton, R. L. (2013). *Incompressible Flow* (4th ed.). Wiley.
- Pope, S. B. (2000). *Turbulent Flows*. Cambridge University Press.
