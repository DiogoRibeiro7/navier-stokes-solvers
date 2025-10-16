---
layout: single
title: "FAQ"
permalink: /faq/
toc: true
---

# Frequently Asked Questions

> For setup steps refer to the [Usage Guide](/usage-guide/); for theoretical insight see [Theory Background](/theory-background/).

## Installation

**Q:** *What dependencies do I need?*

**A:** GCC/Clang, Make, FFTW3 (for spectral solver), Python 3.9+ with packages listed in `scripts/requirements.txt`.

**Q:** *FFTW is missing on my system. What do I do?*

**A:** Install via package manager (e.g., `sudo apt install libfftw3-dev` on Ubuntu) or build from source with `--enable-shared`.

## Compilation

**Q:** *`make spectral` fails with linker errors.*

**A:** Ensure `-lfftw3` is discoverable. Set `LIBRARY_PATH` or edit `Makefile` to add the correct library directory.

**Q:** *Can I use MSVC?*

**A:** Official build system targets POSIX toolchains. On Windows, use WSL or MinGW/MSYS2 with GCC.

## Runtime Issues

**Q:** *`SOLVER_DIVERGED` during FD runs.*

**A:** Reduce `MAX_DT`, refine the grid, or lower Reynolds number. Check that boundary conditions are consistent.

**Q:** *Spectral run produces `nan` values.*

**A:** CFL condition violated. Inspect `ns_spectral_compute_cfl_timestep` output and cap `dt`. Also verify that FFTW wisdom is not corrupted.

## Results Interpretation

**Q:** *Where are velocity and pressure fields stored?*

**A:** Final solutions are written to `results/*.dat` with columns documented in [Usage Guide](usage_guide.md#understanding-output-files).

**Q:** *How do I compare against benchmark data?*

**A:** Use `scripts/compare_methods.py` to compute $L_2$ and $L_\infty$ metrics. See [Benchmark Results](/benchmark-results/) for reference values.

## Performance

**Q:** *How do I speed up spectral runs?*

**A:** Enable FFTW threading, compile with `-O3 -march=native`, and reduce snapshot frequency.

**Q:** *Memory usage spikes on large grids.*

**A:** Pre-allocate arrays once, avoid copying snapshots unnecessarily, and consult memory estimates in [Usage Guide](usage_guide.md#memory-estimates).

## Miscellaneous

**Q:** *Can I run 3D simulations?*

**A:** Not yet; the current codebase targets 2D flows. Contributions are welcome (see [Contributing](/contributing/)).

**Q:** *How do I cite this project?*

**A:** Cite the repository URL and relevant benchmark sections from [Benchmark Results](/benchmark-results/).

Still stuck? Open an issue on GitHub with logs, parameter files, and system details.
