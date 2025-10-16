#!/usr/bin/env python3
"""Advanced plotting utilities for Navier-Stokes solver outputs.

This script generates 2D and 3D visualisations from either the finite
difference (FD) or spectral solver results, producing publication-quality
figures for velocity, pressure, and vorticity fields.

Example
-------
.. code-block:: bash

    python scripts/plot_solution.py results/fd_final_solution.dat \
        --output-prefix figures/fd_final --dpi 400
    python scripts/plot_solution.py results/spectral_final_solution.dat \
        --solver-type spectral --surface-view 35 -55
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from tqdm import tqdm

LOGGER = logging.getLogger(__name__)


@dataclass
class SolutionGrid:
    """Structured representation of solver output on a 2D grid."""

    grid_x: np.ndarray
    grid_y: np.ndarray
    velocity_u: np.ndarray
    velocity_v: np.ndarray
    vorticity: Optional[np.ndarray]
    stream_function: Optional[np.ndarray]
    pressure: Optional[np.ndarray]
    solver_type: str

    @property
    def shape(self) -> Tuple[int, int]:
        """Return grid shape as (ny, nx)."""
        return self.grid_x.shape

    @property
    def velocity_magnitude(self) -> np.ndarray:
        """Compute velocity magnitude."""
        return np.hypot(self.velocity_u, self.velocity_v)

    @property
    def x_coords(self) -> np.ndarray:
        """Return unique x coordinates."""
        return self.grid_x[0, :]

    @property
    def y_coords(self) -> np.ndarray:
        """Return unique y coordinates."""
        return self.grid_y[:, 0]

    def ensure_vorticity(self) -> np.ndarray:
        """Return vorticity field, computing it if absent."""
        if self.vorticity is not None:
            return self.vorticity

        if min(self.shape) < 2:
            raise ValueError("Cannot compute vorticity on a degenerate grid.")

        try:
            dvy_dy, dvy_dx = np.gradient(
                self.velocity_v,
                self.y_coords,
                self.x_coords,
                edge_order=2,
            )
            du_dy, du_dx = np.gradient(
                self.velocity_u,
                self.y_coords,
                self.x_coords,
                edge_order=2,
            )
        except ValueError:
            dvy_dy, dvy_dx = np.gradient(self.velocity_v)
            du_dy, du_dx = np.gradient(self.velocity_u)

        self.vorticity = dvy_dx - du_dy
        return self.vorticity


def configure_matplotlib(font_size: float, dpi: int) -> None:
    """Apply styling suitable for publication-quality figures."""
    mpl.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "font.size": font_size,
            "axes.titlesize": font_size + 2,
            "axes.labelsize": font_size + 1,
            "legend.fontsize": font_size - 1,
            "xtick.labelsize": font_size - 1,
            "ytick.labelsize": font_size - 1,
            "image.cmap": "viridis",
        }
    )
    plt.style.use("seaborn-v0_8-paper")


def reshape_grid(
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    nx: int,
    ny: int,
) -> np.ndarray:
    """Reshape a flattened field onto a structured grid."""
    order = np.lexsort((x, y))
    return values[order].reshape(ny, nx)


def load_solution(
    filename: Path,
    solver_type: str = "auto",
    grid_shape: Optional[Tuple[int, int]] = None,
) -> SolutionGrid:
    """Load solver output and reshape fields onto a 2D grid."""
    if not filename.exists():
        raise FileNotFoundError(f"Input file not found: {filename}")

    try:
        data = np.loadtxt(filename, comments="#")
    except ValueError as exc:
        raise ValueError(f"Failed to parse data from {filename}: {exc}") from exc

    if data.ndim != 2 or data.shape[1] < 4:
        raise ValueError(
            "Expected at least four columns (x, y, u, v) in solution file."
        )

    x = data[:, 0]
    y = data[:, 1]
    u = data[:, 2]
    v = data[:, 3]

    if grid_shape:
        nx, ny = grid_shape
        if nx * ny != x.size:
            raise ValueError("Provided grid dimensions do not match data size.")
    else:
        unique_x = np.unique(x)
        unique_y = np.unique(y)
        nx = unique_x.size
        ny = unique_y.size
        if nx * ny != x.size:
            raise ValueError(
                "Unable to infer structured grid dimensions from data."
            )

    grid_x = reshape_grid(x, y, x, nx, ny)
    grid_y = reshape_grid(x, y, y, nx, ny)
    velocity_u = reshape_grid(x, y, u, nx, ny)
    velocity_v = reshape_grid(x, y, v, nx, ny)

    omega = reshape_grid(x, y, data[:, 4], nx, ny) if data.shape[1] > 4 else None
    psi = reshape_grid(x, y, data[:, 5], nx, ny) if data.shape[1] > 5 else None
    pressure = (
        reshape_grid(x, y, data[:, 6], nx, ny) if data.shape[1] > 6 else None
    )

    if solver_type == "auto":
        if omega is not None and pressure is None:
            solver = "spectral"
        elif pressure is not None and omega is None:
            solver = "fd"
        elif omega is not None and pressure is not None:
            solver = "hybrid"
        else:
            solver = "fd"
    else:
        solver = solver_type

    return SolutionGrid(
        grid_x=grid_x,
        grid_y=grid_y,
        velocity_u=velocity_u,
        velocity_v=velocity_v,
        vorticity=omega,
        stream_function=psi,
        pressure=pressure,
        solver_type=solver,
    )


def add_colorbar(ax: plt.Axes, im: mpl.cm.ScalarMappable, label: str) -> None:
    """Attach a neatly sized colorbar to an axes."""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.08)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(label)
    cbar.ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))


def plot_velocity_field(
    solution: SolutionGrid,
    output_path: Path,
    contour_levels: int,
    colormap: str,
    quiver_step: int,
    streamline_density: float,
) -> None:
    """Plot 2D velocity magnitude with vectors and streamlines."""
    magnitude = solution.velocity_magnitude
    fig, ax = plt.subplots(figsize=(7, 6))

    im = ax.contourf(
        solution.grid_x,
        solution.grid_y,
        magnitude,
        levels=contour_levels,
        cmap=colormap,
    )

    step = max(1, quiver_step)
    ax.quiver(
        solution.grid_x[::step, ::step],
        solution.grid_y[::step, ::step],
        solution.velocity_u[::step, ::step],
        solution.velocity_v[::step, ::step],
        color="white",
        alpha=0.75,
        linewidth=0.5,
        scale_units="xy",
        scale=None,
    )

    ax.streamplot(
        solution.x_coords,
        solution.y_coords,
        solution.velocity_u,
        solution.velocity_v,
        color="k",
        density=streamline_density,
        linewidth=0.8,
        arrowsize=1.0,
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Velocity Magnitude and Streamlines")
    ax.set_aspect("equal")
    add_colorbar(ax, im, label=r"$|\mathbf{u}|$")

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_velocity_surface(
    solution: SolutionGrid,
    output_path: Path,
    colormap: str,
    surface_view: Tuple[float, float],
) -> None:
    """Generate a 3D surface plot of the velocity magnitude."""
    magnitude = solution.velocity_magnitude
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        solution.grid_x,
        solution.grid_y,
        magnitude,
        cmap=colormap,
        linewidth=0,
        antialiased=True,
        shade=True,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel(r"$|\mathbf{u}|$")
    ax.view_init(elev=surface_view[0], azim=surface_view[1])
    fig.colorbar(surf, shrink=0.7, aspect=12, pad=0.08, label=r"$|\mathbf{u}|$")
    ax.set_title("Velocity Magnitude Surface")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_pressure_contours(
    solution: SolutionGrid,
    output_path: Path,
    contour_levels: int,
    colormap: str,
) -> None:
    """Plot filled contours of the pressure field."""
    if solution.pressure is None:
        raise ValueError("Pressure data not available for this solution.")

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.contourf(
        solution.grid_x,
        solution.grid_y,
        solution.pressure,
        levels=contour_levels,
        cmap=colormap,
    )
    ax.contour(
        solution.grid_x,
        solution.grid_y,
        solution.pressure,
        levels=10,
        colors="k",
        linewidths=0.6,
        alpha=0.6,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Pressure Contours")
    ax.set_aspect("equal")
    add_colorbar(ax, im, label="Pressure")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_vorticity_isolines(
    solution: SolutionGrid,
    output_path: Path,
    contour_levels: int,
    colormap: str,
) -> None:
    """Plot vorticity field with isolines."""
    vorticity = solution.ensure_vorticity()
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.contourf(
        solution.grid_x,
        solution.grid_y,
        vorticity,
        levels=contour_levels,
        cmap=colormap,
    )
    ax.contour(
        solution.grid_x,
        solution.grid_y,
        vorticity,
        levels=8,
        colors="k",
        linewidths=0.6,
        alpha=0.7,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Vorticity Isolines")
    ax.set_aspect("equal")
    add_colorbar(ax, im, label=r"$\omega$")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_centerline_profiles(
    solution: SolutionGrid,
    output_path: Path,
) -> None:
    """Plot velocity profiles along central lines of the domain."""
    ny, nx = solution.shape
    if nx < 2 or ny < 2:
        raise ValueError("Insufficient points for centerline profiles.")

    mid_x = nx // 2
    mid_y = ny // 2

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].plot(
        solution.velocity_u[:, mid_x],
        solution.y_coords,
        label="u(x=mid)",
        linewidth=2.0,
    )
    axes[0].set_xlabel("u")
    axes[0].set_ylabel("y")
    axes[0].set_title("Vertical Centerline")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(
        solution.x_coords,
        solution.velocity_v[mid_y, :],
        label="v(y=mid)",
        color="tab:red",
        linewidth=2.0,
    )
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("v")
    axes[1].set_title("Horizontal Centerline")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def build_argument_parser() -> argparse.ArgumentParser:
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate advanced plots from Navier-Stokes solver outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Path to solver output file (.dat).",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=Path,
        help="Output filename prefix. Defaults to input stem.",
    )
    parser.add_argument(
        "--solver-type",
        choices=("auto", "fd", "spectral", "hybrid"),
        default="auto",
        help="Solver type used to interpret optional fields.",
    )
    parser.add_argument(
        "--grid",
        type=int,
        nargs=2,
        metavar=("NX", "NY"),
        help="Override grid dimensions when automatic inference fails.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=400,
        help="Dots per inch for saved figures.",
    )
    parser.add_argument(
        "--font-size",
        type=float,
        default=11.0,
        help="Base font size for figure text.",
    )
    parser.add_argument(
        "--colormap",
        default="viridis",
        help="Colormap for velocity magnitude visualisations.",
    )
    parser.add_argument(
        "--pressure-colormap",
        default="magma",
        help="Colormap for pressure contours.",
    )
    parser.add_argument(
        "--vorticity-colormap",
        default="coolwarm",
        help="Colormap for vorticity isolines.",
    )
    parser.add_argument(
        "--contour-levels",
        type=int,
        default=30,
        help="Number of contour levels for filled contour plots.",
    )
    parser.add_argument(
        "--quiver-step",
        type=int,
        default=3,
        help="Subsampling factor for quiver arrows (>=1).",
    )
    parser.add_argument(
        "--streamline-density",
        type=float,
        default=1.8,
        help="Streamline density parameter for matplotlib.streamplot.",
    )
    parser.add_argument(
        "--surface-view",
        type=float,
        nargs=2,
        metavar=("ELEV", "AZIM"),
        default=(30.0, -60.0),
        help="Elevation and azimuth angles for the 3D surface plot.",
    )
    parser.add_argument(
        "--skip-profiles",
        action="store_true",
        help="Skip centerline profile plots.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bar output.",
    )
    return parser


def main(args: Optional[Sequence[str]] = None) -> int:
    """CLI entry point."""
    parser = build_argument_parser()
    parsed = parser.parse_args(args=args)

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    output_prefix = (
        parsed.output_prefix
        if parsed.output_prefix is not None
        else parsed.input_file.with_suffix("")
    )
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    configure_matplotlib(parsed.font_size, parsed.dpi)

    try:
        solution = load_solution(
            parsed.input_file,
            solver_type=parsed.solver_type,
            grid_shape=tuple(parsed.grid) if parsed.grid else None,
        )
    except (FileNotFoundError, ValueError) as exc:
        LOGGER.error("%s", exc)
        return 1

    LOGGER.info(
        "Loaded %s data on grid %s",
        solution.solver_type.upper(),
        solution.shape,
    )

    tasks = []
    tasks.append(
        (
            "velocity field",
            lambda: plot_velocity_field(
                solution=solution,
                output_path=output_prefix.with_name(
                    f"{output_prefix.name}_velocity.png"
                ),
                contour_levels=parsed.contour_levels,
                colormap=parsed.colormap,
                quiver_step=parsed.quiver_step,
                streamline_density=parsed.streamline_density,
            ),
        )
    )
    tasks.append(
        (
            "velocity surface",
            lambda: plot_velocity_surface(
                solution=solution,
                output_path=output_prefix.with_name(
                    f"{output_prefix.name}_velocity_surface.png"
                ),
                colormap=parsed.colormap,
                surface_view=tuple(parsed.surface_view),
            ),
        )
    )
    if solution.pressure is not None:
        tasks.append(
            (
                "pressure contours",
                lambda: plot_pressure_contours(
                    solution=solution,
                    output_path=output_prefix.with_name(
                        f"{output_prefix.name}_pressure.png"
                    ),
                    contour_levels=parsed.contour_levels,
                    colormap=parsed.pressure_colormap,
                ),
            )
        )
    else:
        LOGGER.info("Pressure data not found: skipping pressure contours.")

    try:
        solution.ensure_vorticity()
        tasks.append(
            (
                "vorticity isolines",
                lambda: plot_vorticity_isolines(
                    solution=solution,
                    output_path=output_prefix.with_name(
                        f"{output_prefix.name}_vorticity.png"
                    ),
                    contour_levels=parsed.contour_levels,
                    colormap=parsed.vorticity_colormap,
                ),
            )
        )
    except ValueError as exc:
        LOGGER.info("Skipping vorticity isolines: %s", exc)

    if not parsed.skip_profiles:
        tasks.append(
            (
                "centerline profiles",
                lambda: plot_centerline_profiles(
                    solution=solution,
                    output_path=output_prefix.with_name(
                        f"{output_prefix.name}_profiles.png"
                    ),
                ),
            )
        )

    iterator = tasks
    reporter = LOGGER.info
    if not parsed.no_progress:
        iterator = tqdm(tasks, desc="Generating figures", unit="figure")
        reporter = iterator.write

    saved_files = []
    for label, task in iterator:
        try:
            task()
            out_file = output_prefix.with_name(
                f"{output_prefix.name}_{label.replace(' ', '_')}.png"
            )
            saved_files.append(out_file)
            reporter(f"[OK] Saved {label}: {out_file}")
        except Exception as exc:  # pylint: disable=broad-except
            reporter(f"[WARN] Failed to create {label}: {exc}")

    if not parsed.no_progress and isinstance(iterator, tqdm):
        iterator.close()

    if saved_files:
        LOGGER.info("Generated %d figure(s).", len(saved_files))
    else:
        LOGGER.warning("No figures were generated.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
