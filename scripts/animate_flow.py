#!/usr/bin/env python3
"""Create animations of Navier-Stokes flow fields from time-series data.

The animation combines a velocity magnitude colormap, quiver vectors, and
Lagrangian particle tracers advected by the velocity field. Frames are exported
to GIF or MP4 formats with optional time annotations.

Example
-------
.. code-block:: bash

    python scripts/animate_flow.py results/spectral_step_*.dat \
        -o animations/taylor_green.mp4 --fps 24 --particle-count 200
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm

from plot_solution import SolutionGrid, load_solution

LOGGER = logging.getLogger(__name__)


def expand_inputs(patterns: Iterable[str]) -> List[Path]:
    """Resolve glob patterns and return a sorted list of frame files."""
    files: List[Path] = []
    for pattern in patterns:
        matches = sorted(Path().glob(pattern)) if any(
            char in pattern for char in "*?["
        ) else [Path(pattern)]
        if not matches:
            LOGGER.warning("Pattern produced no matches: %s", pattern)
        files.extend(match for match in matches if match.is_file())
    unique_files = sorted(set(files))
    if not unique_files:
        raise FileNotFoundError("No input files matched the provided patterns.")
    return unique_files


def load_frames(
    paths: Sequence[Path],
    solver_type: str,
    show_progress: bool,
) -> List[SolutionGrid]:
    """Load solution frames, preserving order."""
    iterator: Iterable[Path]
    if show_progress:
        iterator = tqdm(paths, desc="Loading frames", unit="frame")
    else:
        iterator = paths

    frames: List[SolutionGrid] = []
    for frame_path in iterator:
        try:
            frames.append(load_solution(frame_path, solver_type=solver_type))
        except ValueError as exc:
            LOGGER.error("Failed to load %s: %s", frame_path, exc)
            raise
    return frames


def initialise_matplotlib(dpi: int, font_size: float) -> None:
    """Set plotting defaults suitable for animations."""
    mpl.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "font.size": font_size,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size + 1,
            "xtick.labelsize": font_size - 1,
            "ytick.labelsize": font_size - 1,
        }
    )
    plt.style.use("seaborn-v0_8-darkgrid")


def create_interpolators(solution: SolutionGrid) -> Tuple[
    RegularGridInterpolator,
    RegularGridInterpolator,
]:
    """Build interpolators for velocity components."""
    grid_y = solution.y_coords
    grid_x = solution.x_coords
    return (
        RegularGridInterpolator((grid_y, grid_x), solution.velocity_u),
        RegularGridInterpolator((grid_y, grid_x), solution.velocity_v),
    )


def wrap_positions(
    positions: np.ndarray,
    x_bounds: Tuple[float, float],
    y_bounds: Tuple[float, float],
) -> np.ndarray:
    """Wrap particle positions into the periodic simulation domain."""
    width = x_bounds[1] - x_bounds[0]
    height = y_bounds[1] - y_bounds[0]
    positions[:, 0] = x_bounds[0] + np.mod(positions[:, 0] - x_bounds[0], width)
    positions[:, 1] = y_bounds[0] + np.mod(positions[:, 1] - y_bounds[0], height)
    return positions


def clamp_positions(
    positions: np.ndarray,
    x_bounds: Tuple[float, float],
    y_bounds: Tuple[float, float],
) -> np.ndarray:
    """Clamp particle positions to the simulation domain boundaries."""
    positions[:, 0] = np.clip(positions[:, 0], *x_bounds)
    positions[:, 1] = np.clip(positions[:, 1], *y_bounds)
    return positions


def advect_particles(
    positions: np.ndarray,
    interpolators: Tuple[RegularGridInterpolator, RegularGridInterpolator],
    dt: float,
    x_bounds: Tuple[float, float],
    y_bounds: Tuple[float, float],
    periodic: bool,
) -> np.ndarray:
    """Advance particle positions using Euler integration."""
    u_interp, v_interp = interpolators
    velocity = np.column_stack(
        (u_interp(positions), v_interp(positions))
    )
    new_positions = positions + dt * velocity
    if periodic:
        return wrap_positions(new_positions, x_bounds, y_bounds)
    return clamp_positions(new_positions, x_bounds, y_bounds)


def build_animation(
    frames: Sequence[SolutionGrid],
    output_path: Path,
    fps: int,
    dpi: int,
    particle_count: int,
    dt: float,
    time_start: float,
    cmap: str,
    quiver_step: int,
    tracer_size: float,
    periodic: bool,
    show_progress: bool,
    bitrate: Optional[int],
) -> None:
    """Create and save animation from loaded frames."""
    first = frames[0]
    initialise_matplotlib(dpi=dpi, font_size=11.0)

    extent = (
        float(first.x_coords.min()),
        float(first.x_coords.max()),
        float(first.y_coords.min()),
        float(first.y_coords.max()),
    )

    vmin = min(frame.velocity_magnitude.min() for frame in frames)
    vmax = max(frame.velocity_magnitude.max() for frame in frames)

    fig, ax = plt.subplots(figsize=(7.5, 6.0))

    image = ax.imshow(
        first.velocity_magnitude,
        origin="lower",
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
    )
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$|\mathbf{u}|$")

    step = max(1, quiver_step)
    quiver = ax.quiver(
        first.grid_x[::step, ::step],
        first.grid_y[::step, ::step],
        first.velocity_u[::step, ::step],
        first.velocity_v[::step, ::step],
        color="white",
        scale_units="xy",
        scale=None,
        alpha=0.8,
        linewidth=0.5,
    )

    rng = np.random.default_rng(seed=42)
    particles = np.column_stack(
        (
            rng.uniform(extent[0], extent[1], particle_count),
            rng.uniform(extent[2], extent[3], particle_count),
        )
    )
    scatter = ax.scatter(
        particles[:, 0],
        particles[:, 1],
        s=tracer_size,
        c="white",
        edgecolor="black",
        linewidth=0.3,
        alpha=0.85,
    )

    times = time_start + dt * np.arange(len(frames))
    time_text = ax.text(
        0.02,
        0.95,
        "",
        transform=ax.transAxes,
        color="white",
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.4),
    )

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Velocity Field Evolution with Particle Tracers")

    interpolators = [create_interpolators(frame) for frame in frames]
    x_bounds = (extent[0], extent[1])
    y_bounds = (extent[2], extent[3])

    def init_func():
        image.set_data(first.velocity_magnitude)
        quiver.set_UVC(
            first.velocity_u[::step, ::step],
            first.velocity_v[::step, ::step],
        )
        scatter.set_offsets(particles)
        time_text.set_text(f"t = {times[0]:.3f}")
        return image, quiver, scatter, time_text

    def update(frame_index: int):
        frame = frames[frame_index]
        image.set_data(frame.velocity_magnitude)
        quiver.set_UVC(
            frame.velocity_u[::step, ::step],
            frame.velocity_v[::step, ::step],
        )
        scatter.set_offsets(
            advect_particles(
                scatter.get_offsets(),
                interpolators[frame_index],
                dt,
                x_bounds,
                y_bounds,
                periodic,
            )
        )
        time_text.set_text(f"t = {times[frame_index]:.3f}")
        return image, quiver, scatter, time_text

    anim = animation.FuncAnimation(
        fig,
        update,
        init_func=init_func,
        frames=len(frames),
        blit=False,
        interval=1000.0 / fps,
    )

    writer: animation.AbstractMovieWriter
    suffix = output_path.suffix.lower()
    if suffix in {".gif"}:
        writer = animation.PillowWriter(fps=fps)
    else:
        writer = animation.FFMpegWriter(fps=fps, bitrate=bitrate)

    progress_bar: Optional[tqdm] = None
    progress_callback = None
    if show_progress:
        progress_bar = tqdm(total=len(frames), desc="Rendering", unit="frame")

        def _callback(frame_number: int, total: int) -> None:
            if progress_bar is not None:
                progress_bar.update()

        progress_callback = _callback

    try:
        anim.save(
            output_path,
            writer=writer,
            dpi=dpi,
            progress_callback=progress_callback,
        )
        LOGGER.info("Saved animation to %s", output_path)
    except RuntimeError as exc:
        LOGGER.error("Failed to save animation: %s", exc)
        raise
    finally:
        if progress_bar is not None:
            progress_bar.close()
        plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    """Construct CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Animate Navier-Stokes solver outputs (FD or spectral).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input files or glob patterns sorted lexicographically.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output animation file (.mp4 or .gif).",
    )
    parser.add_argument(
        "--solver-type",
        default="auto",
        choices=("auto", "fd", "spectral", "hybrid"),
        help="Solver type hint for optional columns.",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=20,
        help="Frames per second for animation.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Resolution of output frames.",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=0.05,
        help="Time increment between frames for annotation.",
    )
    parser.add_argument(
        "--time-start",
        type=float,
        default=0.0,
        help="Starting time for the animation overlay.",
    )
    parser.add_argument(
        "--colormap",
        default="plasma",
        help="Colormap for velocity magnitude.",
    )
    parser.add_argument(
        "--quiver-step",
        type=int,
        default=4,
        help="Subsampling factor for quiver arrows (>=1).",
    )
    parser.add_argument(
        "--particle-count",
        type=int,
        default=150,
        help="Number of tracer particles.",
    )
    parser.add_argument(
        "--particle-size",
        type=float,
        default=25.0,
        help="Marker size for tracer particles.",
    )
    parser.add_argument(
        "--periodic",
        action="store_true",
        help="Wrap tracer particles periodically instead of clamping.",
    )
    parser.add_argument(
        "--bitrate",
        type=int,
        default=1800,
        help="Bitrate for MP4 encoding (ignored for GIF).",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars during loading and rendering.",
    )
    return parser


def main(args: Optional[Sequence[str]] = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    parsed = parser.parse_args(args=args)

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    try:
        paths = expand_inputs(parsed.inputs)
    except FileNotFoundError as exc:
        LOGGER.error("%s", exc)
        return 1

    show_progress = not parsed.no_progress
    try:
        frames = load_frames(paths, parsed.solver_type, show_progress)
    except ValueError:
        return 1

    output_path = parsed.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        build_animation(
            frames=frames,
            output_path=output_path,
            fps=parsed.fps,
            dpi=parsed.dpi,
            particle_count=parsed.particle_count,
            dt=parsed.dt,
            time_start=parsed.time_start,
            cmap=parsed.colormap,
            quiver_step=parsed.quiver_step,
            tracer_size=parsed.particle_size,
            periodic=parsed.periodic,
            show_progress=show_progress,
            bitrate=parsed.bitrate,
        )
    except RuntimeError:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
