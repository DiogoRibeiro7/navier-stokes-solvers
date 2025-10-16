#!/usr/bin/env python3
"""Compute and plot the 1D energy spectrum from spectral solver output.

This script performs a 2D Fourier transform of the velocity field, derives the
isotropic energy spectrum E(k), and highlights potential inertial-range scaling
regions with a reference -5/3 slope.

Example
-------
.. code-block:: bash

    python scripts/energy_spectrum.py results/spectral_final_solution.dat \
        -o figures/energy_spectrum.png
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy import fft as spfft
from tqdm import tqdm

from plot_solution import SolutionGrid, configure_matplotlib, load_solution

LOGGER = logging.getLogger(__name__)


def compute_energy_spectrum(solution: SolutionGrid) -> Tuple[np.ndarray, np.ndarray]:
    """Compute isotropic 1D energy spectrum from velocity field."""
    ny, nx = solution.shape
    if nx < 2 or ny < 2:
        raise ValueError("Need at least a 2x2 grid to compute energy spectrum.")

    dx = np.mean(np.diff(solution.x_coords))
    dy = np.mean(np.diff(solution.y_coords))

    u_hat = spfft.fft2(solution.velocity_u, norm="forward")
    v_hat = spfft.fft2(solution.velocity_v, norm="forward")
    energy_density = 0.5 * (np.abs(u_hat) ** 2 + np.abs(v_hat) ** 2)

    kx = spfft.fftfreq(nx, d=dx) * 2.0 * np.pi
    ky = spfft.fftfreq(ny, d=dy) * 2.0 * np.pi
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    k_magnitude = np.sqrt(kx_grid**2 + ky_grid**2)

    k_flat = k_magnitude.ravel()
    energy_flat = energy_density.ravel()

    dkx = np.min(np.abs(np.diff(np.sort(kx)))) if nx > 1 else 0.0
    dky = np.min(np.abs(np.diff(np.sort(ky)))) if ny > 1 else 0.0
    dk_candidates = [value for value in (dkx, dky) if value > 0]
    dk = min(dk_candidates) if dk_candidates else 1.0
    max_k = k_flat.max()
    bins = np.arange(0.0, max_k + dk, dk)

    bin_indices = np.digitize(k_flat, bins) - 1
    energy_sum = np.bincount(bin_indices, weights=energy_flat, minlength=len(bins))
    counts = np.bincount(bin_indices, minlength=len(bins))

    valid = counts > 0
    k_centers = 0.5 * (bins[:-1] + bins[1:])
    spectrum = np.zeros_like(k_centers)
    spectrum[valid[:-1]] = energy_sum[:-1][valid[:-1]] / counts[:-1][valid[:-1]]

    mask = k_centers > 0
    return k_centers[mask], spectrum[mask]


def detect_inertial_range(
    k: np.ndarray,
    spectrum: np.ndarray,
    slope_target: float,
    tolerance: float,
    window: int,
    show_progress: bool,
) -> Optional[Tuple[int, int, float]]:
    """Identify a contiguous wavenumber interval exhibiting target slope."""
    if len(k) < window:
        return None

    logk = np.log(k)
    logE = np.log(spectrum)

    iterator = range(len(k) - window)
    if show_progress:
        iterator = tqdm(iterator, desc="Scanning inertial range", leave=False)
    for start in iterator:
        stop = start + window
        slope, _ = np.polyfit(logk[start:stop], logE[start:stop], deg=1)
        if abs(slope - slope_target) < tolerance:
            return start, stop, slope
    return None


def plot_spectrum(
    k: np.ndarray,
    spectrum: np.ndarray,
    inertial: Optional[Tuple[int, int, float]],
    output_path: Path,
    dpi: int,
) -> None:
    """Render log-log plot of energy spectrum with reference slope."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.loglog(k, spectrum, "o-", label="E(k)")

    if inertial:
        start, stop, slope = inertial
        ax.fill_between(
            k[start:stop],
            spectrum[start:stop] * 0.9,
            spectrum[start:stop] * 1.1,
            color="tab:orange",
            alpha=0.3,
            label=f"Inertial range (slope ~= {slope:.2f})",
        )
        anchor_k = k[start]
        anchor_E = spectrum[start]
    else:
        anchor_k = k[len(k) // 3]
        anchor_E = spectrum[len(k) // 3]

    reference = anchor_E * (k / anchor_k) ** (-5.0 / 3.0)
    ax.loglog(k, reference, "--", label="Reference slope -5/3")

    ax.set_xlabel("Wavenumber k")
    ax.set_ylabel("E(k)")
    ax.set_title("Energy Spectrum")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend()

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    """CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Compute isotropic energy spectrum from spectral solutions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Spectral solution file containing velocity components.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("energy_spectrum.png"),
        help="Output figure path.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=400,
        help="Figure resolution.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress indicator for inertial range detection.",
    )
    return parser


def main(args: Optional[Sequence[str]] = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    parsed = parser.parse_args(args=args)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    configure_matplotlib(font_size=11.0, dpi=parsed.dpi)

    try:
        solution = load_solution(parsed.input_file, solver_type="spectral")
    except (FileNotFoundError, ValueError) as exc:
        LOGGER.error("%s", exc)
        return 1

    try:
        k, spectrum = compute_energy_spectrum(solution)
    except ValueError as exc:
        LOGGER.error("Spectrum computation failed: %s", exc)
        return 1

    inertial = detect_inertial_range(
        k,
        spectrum,
        slope_target=-5.0 / 3.0,
        tolerance=0.4,
        window=max(5, len(k) // 12),
        show_progress=not parsed.no_progress,
    )

    output_path = parsed.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_spectrum(k, spectrum, inertial, output_path, parsed.dpi)
    LOGGER.info("Saved spectrum plot to %s", output_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
