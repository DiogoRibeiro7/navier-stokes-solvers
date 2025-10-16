#!/usr/bin/env python3
"""Compare finite difference and spectral Navier-Stokes solver outputs.

The script aligns both solutions on a common grid, evaluates error metrics, and
produces publication-ready figures illustrating similarities and differences.
Optional convergence plots quantify how the FD solution approaches the spectral
reference as the mesh is refined.

Example
-------
.. code-block:: bash

    python scripts/compare_methods.py results/fd_final_solution.dat \
        results/spectral_final_solution.dat \
        -o figures/fd_vs_spectral --fd-series results/fd_*.dat
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import linregress
from tqdm import tqdm

from plot_solution import SolutionGrid, configure_matplotlib, load_solution

LOGGER = logging.getLogger(__name__)


def build_interpolator(field: np.ndarray, solution: SolutionGrid) -> RegularGridInterpolator:
    """Create interpolator for a scalar field defined on solution grid."""
    return RegularGridInterpolator(
        (solution.y_coords, solution.x_coords),
        field,
        bounds_error=False,
        fill_value=None,
    )


def map_solution(
    source: SolutionGrid,
    target: SolutionGrid,
) -> Dict[str, np.ndarray]:
    """Resample source fields onto the target grid."""
    points = np.column_stack(
        (target.grid_y.ravel(), target.grid_x.ravel())
    )

    u_interp = build_interpolator(source.velocity_u, source)
    v_interp = build_interpolator(source.velocity_v, source)

    mapped_u = u_interp(points).reshape(target.shape)
    mapped_v = v_interp(points).reshape(target.shape)

    mapped_pressure = None
    if source.pressure is not None:
        p_interp = build_interpolator(source.pressure, source)
        mapped_pressure = p_interp(points).reshape(target.shape)

    mapped_vorticity = None
    try:
        vorticity = source.ensure_vorticity()
        omega_interp = build_interpolator(vorticity, source)
        mapped_vorticity = omega_interp(points).reshape(target.shape)
    except ValueError:
        mapped_vorticity = None

    return {
        "u": mapped_u,
        "v": mapped_v,
        "pressure": mapped_pressure,
        "vorticity": mapped_vorticity,
    }


def compute_errors(
    reference: SolutionGrid,
    comparison: Dict[str, np.ndarray],
) -> Dict[str, float]:
    """Compute error metrics between reference solution and comparison fields."""
    metrics: Dict[str, float] = {}

    ref_mag = reference.velocity_magnitude
    cmp_mag = np.hypot(comparison["u"], comparison["v"])
    diff_mag = cmp_mag - ref_mag

    metrics["velocity_l2"] = np.linalg.norm(diff_mag) / np.linalg.norm(ref_mag)
    metrics["velocity_linf"] = np.max(np.abs(diff_mag))

    du = comparison["u"] - reference.velocity_u
    dv = comparison["v"] - reference.velocity_v
    metrics["u_l2"] = np.linalg.norm(du) / np.linalg.norm(reference.velocity_u)
    metrics["v_l2"] = np.linalg.norm(dv) / np.linalg.norm(reference.velocity_v)

    if reference.pressure is not None and comparison["pressure"] is not None:
        dp = comparison["pressure"] - reference.pressure
        metrics["pressure_l2"] = np.linalg.norm(dp) / np.linalg.norm(reference.pressure)
        metrics["pressure_linf"] = np.max(np.abs(dp))

    if comparison["vorticity"] is not None:
        ref_omega = reference.ensure_vorticity()
        d_omega = comparison["vorticity"] - ref_omega
        metrics["vorticity_l2"] = np.linalg.norm(d_omega) / np.linalg.norm(ref_omega)

    return metrics


def format_metrics_table(metrics: Dict[str, float]) -> pd.DataFrame:
    """Return metrics DataFrame for display and export."""
    items = sorted(metrics.items())
    table = pd.DataFrame(items, columns=["metric", "value"])
    return table


def render_velocity_comparison(
    fd_on_ref: Dict[str, np.ndarray],
    reference: SolutionGrid,
    fd_original: SolutionGrid,
    prefix: Path,
    contour_levels: int,
    cmap: str,
    error_cmap: str,
) -> Path:
    """Create side-by-side velocity magnitude comparison."""
    fd_mag = np.hypot(fd_on_ref["u"], fd_on_ref["v"])
    spec_mag = reference.velocity_magnitude
    diff = fd_mag - spec_mag

    vmax = max(fd_mag.max(), spec_mag.max())
    vmin = min(fd_mag.min(), spec_mag.min())
    symmetric = max(abs(diff.min()), abs(diff.max()))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    im0 = axes[0].contourf(
        fd_original.grid_x,
        fd_original.grid_y,
        fd_original.velocity_magnitude,
        levels=contour_levels,
        cmap=cmap,
    )
    axes[0].set_title("Finite Difference |u|")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].contourf(
        reference.grid_x,
        reference.grid_y,
        spec_mag,
        levels=contour_levels,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    axes[1].set_title("Spectral |u|")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    norm = TwoSlopeNorm(vcenter=0.0, vmin=-symmetric, vmax=symmetric)
    im2 = axes[2].contourf(
        reference.grid_x,
        reference.grid_y,
        diff,
        levels=contour_levels,
        cmap=error_cmap,
        norm=norm,
    )
    axes[2].set_title("FD - Spectral |u| Difference")
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04, label="Difference")

    output = prefix.with_name(f"{prefix.name}_velocity_comparison.png")
    fig.savefig(output, dpi=plt.rcParams["savefig.dpi"], bbox_inches="tight")
    plt.close(fig)
    return output


def render_vorticity_comparison(
    fd_on_ref: Dict[str, np.ndarray],
    reference: SolutionGrid,
    prefix: Path,
    contour_levels: int,
    cmap: str,
) -> Optional[Path]:
    """Create vorticity comparison if data are available."""
    if fd_on_ref["vorticity"] is None:
        LOGGER.info("Skipping vorticity comparison (data unavailable).")
        return None

    fd_vorticity = fd_on_ref["vorticity"]
    ref_vorticity = reference.ensure_vorticity()
    diff = fd_vorticity - ref_vorticity
    vmax = np.max(np.abs(diff))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    levels = contour_levels
    im0 = axes[0].contourf(
        reference.grid_x,
        reference.grid_y,
        fd_vorticity,
        levels=levels,
        cmap=cmap,
    )
    axes[0].set_title("Finite Difference omega")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].contourf(
        reference.grid_x,
        reference.grid_y,
        ref_vorticity,
        levels=levels,
        cmap=cmap,
    )
    axes[1].set_title("Spectral omega")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    im2 = axes[2].contourf(
        reference.grid_x,
        reference.grid_y,
        diff,
        levels=levels,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
    )
    axes[2].set_title("Difference omega")
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04, label="Difference")

    output = prefix.with_name(f"{prefix.name}_vorticity_comparison.png")
    fig.savefig(output, dpi=plt.rcParams["savefig.dpi"], bbox_inches="tight")
    plt.close(fig)
    return output


def render_convergence_plot(
    fd_series: Sequence[Path],
    reference: SolutionGrid,
    solver_type: str,
    prefix: Path,
    show_progress: bool,
) -> Optional[Path]:
    """Generate convergence plot from FD series compared to reference."""
    if not fd_series:
        return None

    entries: List[Dict[str, float]] = []
    iterator: Iterable[Path]
    if show_progress:
        iterator = tqdm(fd_series, desc="Evaluating series", unit="file")
    else:
        iterator = fd_series

    for path in iterator:
        candidate = load_solution(path, solver_type=solver_type)
        mapped = map_solution(candidate, reference)
        metrics = compute_errors(reference, mapped)
        resolution = candidate.shape[1]
        entries.append(
            {
                "resolution": resolution,
                "grid_spacing": (reference.x_coords.max() - reference.x_coords.min())
                / (candidate.shape[1] - 1),
                "velocity_l2": metrics["velocity_l2"],
            }
        )

    df = pd.DataFrame(entries).sort_values("resolution")
    if df.empty:
        return None

    slope, intercept, *_ = linregress(
        np.log(df["grid_spacing"]),
        np.log(df["velocity_l2"]),
    )
    fitted = np.exp(intercept) * df["grid_spacing"] ** slope

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.loglog(df["grid_spacing"], df["velocity_l2"], "o-", label="FD error")
    ax.loglog(df["grid_spacing"], fitted, "--", label=f"Fit slope ~= {slope:.2f}")
    ax.set_xlabel("Grid spacing h")
    ax.set_ylabel("Relative L2 error in |u|")
    ax.set_title("Convergence of FD towards Spectral Reference")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend()

    output = prefix.with_name(f"{prefix.name}_convergence.png")
    fig.savefig(output, dpi=plt.rcParams["savefig.dpi"], bbox_inches="tight")
    plt.close(fig)
    return output


def build_parser() -> argparse.ArgumentParser:
    """Create CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Compare finite difference and spectral solver outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fd_file", type=Path, help="Finite difference result file.")
    parser.add_argument("spectral_file", type=Path, help="Spectral result file.")
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=Path,
        default=Path("comparison"),
        help="File prefix for generated figures.",
    )
    parser.add_argument(
        "--fd-series",
        type=Path,
        nargs="*",
        default=[],
        help="Optional list of FD outputs at increasing resolution.",
    )
    parser.add_argument(
        "--solver-type",
        choices=("auto", "fd", "spectral", "hybrid"),
        default="auto",
        help="Hint for parsing FD series files.",
    )
    parser.add_argument(
        "--contour-levels",
        type=int,
        default=30,
        help="Contour levels for surface plots.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=400,
        help="Output figure resolution.",
    )
    parser.add_argument(
        "--colormap",
        default="viridis",
        help="Colormap for magnitude plots.",
    )
    parser.add_argument(
        "--error-colormap",
        default="coolwarm",
        help="Colormap for difference plots.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars for batch operations.",
    )
    return parser


def main(args: Optional[Sequence[str]] = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    parsed = parser.parse_args(args=args)

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    output_prefix = parsed.output_prefix.with_suffix("")
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    configure_matplotlib(font_size=11.0, dpi=parsed.dpi)

    try:
        fd_solution = load_solution(parsed.fd_file, solver_type="fd")
        spectral_solution = load_solution(parsed.spectral_file, solver_type="spectral")
    except (FileNotFoundError, ValueError) as exc:
        LOGGER.error("%s", exc)
        return 1

    fd_on_reference = map_solution(fd_solution, spectral_solution)
    metrics = compute_errors(spectral_solution, fd_on_reference)
    table = format_metrics_table(metrics)
    LOGGER.info("Error metrics:\n%s", table.to_string(index=False, float_format="{:.4e}".format))

    saved_files: List[Path] = []
    saved_files.append(
        render_velocity_comparison(
            fd_on_reference,
            spectral_solution,
            fd_solution,
            output_prefix,
            parsed.contour_levels,
            parsed.colormap,
            parsed.error_colormap,
        )
    )

    vort_file = render_vorticity_comparison(
        fd_on_reference,
        spectral_solution,
        output_prefix,
        parsed.contour_levels,
        parsed.error_colormap,
    )
    if vort_file:
        saved_files.append(vort_file)

    convergence = render_convergence_plot(
        parsed.fd_series,
        spectral_solution,
        solver_type=parsed.solver_type,
        prefix=output_prefix,
        show_progress=not parsed.no_progress,
    )
    if convergence:
        saved_files.append(convergence)

    table_path = output_prefix.with_name(f"{output_prefix.name}_metrics.csv")
    table.to_csv(table_path, index=False)
    saved_files.append(table_path)

    LOGGER.info("Generated %d artefact(s).", len(saved_files))
    for file in saved_files:
        LOGGER.info("  - %s", file)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
