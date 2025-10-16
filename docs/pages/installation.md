---
layout: single
title: "Installation"
permalink: /installation/
toc: true
---

This guide walks through installing Navier-Stokes Solvers on Linux, macOS, and Windows (via WSL), as well as setting up a containerised workflow. Verify your environment with the steps at the end of this page.

## Prerequisites

- C compiler (GCC ≥ 10 or Clang ≥ 13)
- Make
- FFTW3 (double precision)
- Python 3.9 or newer
- Git

Optional: CUDA toolkit (for future GPU extensions), Visual Studio Code for editing, and `gnuplot` or `paraview` for additional visualisation.

## Ubuntu / Debian

```bash
sudo apt update
sudo apt install build-essential git make python3 python3-venv python3-pip libfftw3-dev
```

1. Clone the repository and build the solvers:
   ```bash
   git clone https://github.com/your-org/navier-stokes-solvers.git
   cd navier-stokes-solvers
   make all
   ```
2. Install Python utilities:
   ```bash
   pip install -r scripts/requirements.txt
   ```

## macOS (Homebrew)

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install gcc fftw python git
```

> **Note:** On Apple Silicon, ensure Homebrew installs ARM64 versions (`/opt/homebrew`). Update your `PATH` accordingly.

Build and verify:
```bash
git clone https://github.com/your-org/navier-stokes-solvers.git
cd navier-stokes-solvers
make spectral  # builds FFT-dependent solver
python3 -m pip install -r scripts/requirements.txt
```

## Windows (WSL 2 Recommended)

1. Enable Windows Subsystem for Linux and install Ubuntu from the Microsoft Store.
2. Launch WSL and follow the Ubuntu instructions above.
3. For native Windows builds use MSYS2/MinGW and install packages:
   ```bash
   pacman -S --needed base-devel mingw-w64-x86_64-gcc mingw-w64-x86_64-fftw python
   ```
   Update `Makefile` variables `CC` and `PREFIX` if necessary.

## Docker Workflow

A reproducible Docker environment is provided.

```bash
git clone https://github.com/your-org/navier-stokes-solvers.git
cd navier-stokes-solvers
# Build image
docker build -t navier-stokes-solvers .
# Run container mounting local results
docker run --rm -it -v $(pwd)/results:/workspace/results navier-stokes-solvers ./bin/ns_fd_solver
```

Docker Compose includes service definitions for FD and spectral solvers; run `docker-compose up fd-solver` to start.

## Python Environment (Optional but recommended)

Use a virtual environment to keep dependencies isolated:

```bash
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scriptsctivate
pip install --upgrade pip
pip install -r scripts/requirements.txt
```

## Verification Checklist

1. **Build succeeds:** `make all` completes without errors.
2. **Run smoke tests:**
   ```bash
   ./bin/ns_fd_solver --help || ./bin/ns_fd_solver
   ./bin/ns_spectral_solver --help || ./bin/ns_spectral_solver
   ```
3. **Python tooling:**
   ```bash
   python -m compileall scripts
   python scripts/plot_solution.py --help
   ```
4. **Post-processing:** generate a figure from bundled sample data to confirm Matplotlib is functional:
   ```bash
   python scripts/plot_solution.py tests/data/fd_sample.dat -o figures/quick_check
   ```

If any step fails, consult the [FAQ](/faq/) or open an issue with platform details and logs.
