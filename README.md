# High-Order Navier-Stokes Solvers

Advanced numerical implementations of the 2D incompressible Navier-Stokes equations featuring Newton-Raphson iteration, spectral methods, and comprehensive convergence analysis.

## Author

**Diogo Ribeiro**  
ESMAD - Escola Superior de Média Arte e Design  
Lead Data Scientist at Mysense.ai  
Email: dfr@esmad.ipp.pt  
ORCID: [0009-0001-2022-7072](https://orcid.org/0009-0001-2022-7072)  
GitHub: [@diogoribeiro7](https://github.com/diogoribeiro7)

## Features

### Finite Difference Solver
- **Nonlinear Newton-Raphson iteration** for coupled momentum-continuity equations
- **Adaptive time stepping** based on CFL and viscous stability conditions
- **Sparse Jacobian assembly** (CSR format) for efficient linear algebra
- **Convergence diagnostics** with real-time residual monitoring
- **Lid-driven cavity** and **Taylor-Green vortex** test cases

### Spectral Solver
- **Fourier pseudospectral method** with FFT (exponential convergence)
- **Vorticity-stream function formulation** ensuring incompressibility
- **2/3 dealiasing rule** for accurate nonlinear term computation
- **4th-order Runge-Kutta** time integration
- **Energy cascade analysis** for turbulence diagnostics
- **Spectral resolution monitoring** with adaptive indicators

## Mathematical Background

Both solvers solve the incompressible Navier-Stokes equations:

```
∂u/∂t + (u·∇)u = -∇p + (1/Re)∇²u
∇·u = 0
```

Where:
- **u** = velocity field
- **p** = pressure
- **Re** = Reynolds number

### Finite Difference Method
- 2nd-order central differences for spatial discretization
- Backward Euler for time integration
- Newton-Raphson for nonlinear solver

### Spectral Method
- Fourier basis with exact derivative computation: ∂/∂x → ik_x
- Vorticity transport equation: ∂ω/∂t + u·∇ω = (1/Re)∇²ω
- Poisson equation for stream function: ∇²ψ = -ω

## Requirements

### Compiler
- GCC or Clang with C11 support
- OpenMP support (optional, for parallelization)

### Libraries
- **Standard**: math library (`-lm`)
- **For spectral solver**: FFTW3 library

### Installation of Dependencies

**Ubuntu/Debian:**
```bash
sudo apt-get install gcc libfftw3-dev libfftw3-3
```

**macOS (Homebrew):**
```bash
brew install gcc fftw
```

**Fedora/RHEL:**
```bash
sudo dnf install gcc fftw-devel
```

## Building

Clone and build:
```bash
git clone https://github.com/diogoribeiro7/navier-stokes-solvers.git
cd navier-stokes-solvers
make all
```

Build individual solvers:
```bash
make fd        # Finite difference solver only
make spectral  # Spectral solver only
```

## Usage

### Finite Difference Solver
```bash
./bin/ns_fd_solver
```

**Default parameters:**
- Grid: 32×32
- Reynolds number: 100
- Domain: 1×1 unit square
- Lid-driven cavity flow

**Output files:**
- `fd_solution_step_XXXX.dat` - Periodic snapshots
- `fd_final_solution.dat` - Final solution

### Spectral Solver
```bash
./bin/ns_spectral_solver
```

**Default parameters:**
- Resolution: 128×128 modes
- Reynolds number: 1000
- Domain: 2π×2π (periodic)
- Taylor-Green vortex

**Output files:**
- `spectral_solution_step_XXXX.dat` - Periodic snapshots
- `spectral_final_solution.dat` - Final solution
- Energy spectrum analysis to stdout

## File Structure

```
navier-stokes-solvers/
├── README.md
├── Makefile
├── include/
│   ├── common_types.h         # Common data structures
│   └── time_integration.h     # Time stepping utilities
├── src/
│   ├── finite_difference/
│   │   ├── ns_fd_solver.h
│   │   ├── fd_memory.c
│   │   ├── fd_initialization.c
│   │   ├── fd_newton_raphson.c
│   │   ├── fd_time_advance.c
│   │   ├── fd_analysis.c
│   │   └── main_fd.c
│   └── spectral/
│       ├── ns_spectral_solver.h
│       ├── spectral_memory.c
│       ├── spectral_initialization.c
│       ├── spectral_transforms.c
│       ├── spectral_core.c
│       ├── spectral_time_integration.c
│       ├── spectral_analysis.c
│       ├── spectral_output.c
│       └── main_spectral.c
├── bin/                        # Compiled executables
└── obj/                        # Object files
```

## Output Format

Both solvers produce space-delimited text files:

```
# x y u v [omega] [psi] [p]
0.000000 0.000000 0.0000e+00 0.0000e+00 ...
...
```

Columns:
- **x, y**: Grid coordinates
- **u, v**: Velocity components
- **omega**: Vorticity (spectral solver)
- **psi**: Stream function (spectral solver)
- **p**: Pressure (FD solver)

## Visualization

Use the provided Python scripts (requires matplotlib):

```bash
python scripts/plot_solution.py fd_final_solution.dat
python scripts/visualize_spectrum.py spectral_final_solution.dat
```

Or use gnuplot:
```gnuplot
set pm3d map
splot 'final_solution.dat' using 1:2:3
```

## Performance Benchmarks

Typical performance on Intel Core i7 (single thread):

| Solver | Grid/Modes | Steps/sec | Time to t=1.0 |
|--------|-----------|-----------|---------------|
| FD     | 32×32     | ~50       | ~20s          |
| FD     | 64×64     | ~12       | ~80s          |
| Spectral | 128×128 | ~15       | ~65s          |
| Spectral | 256×256 | ~3        | ~330s         |

## Validation

### Finite Difference Solver
- **Lid-driven cavity** benchmark (Re=100, 400, 1000)
- Comparison with Ghia et al. (1982) reference solutions
- Mass conservation: max(∇·u) < 1e-8

### Spectral Solver
- **Taylor-Green vortex** exact solution decay
- Enstrophy dissipation rate validation
- Spectral accuracy: exponential decay of high-k modes
- Energy spectrum: -5/3 slope in inertial range (high Re)

## Mathematical Techniques Demonstrated

1. **Nonlinear PDE Solving**: Newton-Raphson iteration with Jacobian assembly
2. **Spectral Methods**: Fourier pseudospectral with FFT
3. **Adaptive Algorithms**: CFL-based time step control
4. **Dealiasing**: 2/3 rule for nonlinear terms
5. **Vorticity Formulation**: Automatic incompressibility satisfaction
6. **Convergence Analysis**: Real-time monitoring and diagnostics

## References

1. **Ghia, U., Ghia, K. N., & Shin, C. T.** (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. *Journal of Computational Physics*, 48(3), 387-411.

2. **Canuto, C., Hussaini, M. Y., Quarteroni, A., & Zang, T. A.** (2007). *Spectral Methods: Fundamentals in Single Domains*. Springer.

3. **Peyret, R.** (2002). *Spectral Methods for Incompressible Viscous Flow*. Springer.

4. **Pope, S. B.** (2000). *Turbulent Flows*. Cambridge University Press.

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Submit a pull request

## Contact

For questions or collaborations:
- Email: dfr@esmad.ipp.pt
- GitHub: [@diogoribeiro7](https://github.com/diogoribeiro7)

## Acknowledgments

This work showcases advanced numerical methods for computational fluid dynamics, combining mathematical rigor with efficient C implementation. The solvers demonstrate both classical finite difference techniques and modern spectral methods used in high-performance CFD codes.

---

*Developed as a demonstration of advanced numerical analysis and scientific computing expertise.*
