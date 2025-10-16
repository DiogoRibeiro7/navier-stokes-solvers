#!/usr/bin/env python3
"""Automate convergence studies for Navier-Stokes solver outputs.

The script can launch solver executables across multiple grid resolutions or
ingest precomputed results, compute relative L2 errors against a reference
solution, and produce log-log convergence plots alongside tabulated metrics.

Example
-------
.. code-block:: bash

    python scripts/convergence_study.py --resolutions 32 64 128 \
        --fd-command "../bin/ns_fd_solver --nx {nx} --ny {ny} --output {output}" \
        --spectral-results-pattern "results/spectral_final_{nx:04d}.dat" \
        --output-dir studies/taylor_green
"""

from __future__ import annotations

import argparse
import logging
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

from plot_solution import SolutionGrid, configure_matplotlib, load_solution
from compare_methods import compute_errors, map_solution

LOGGER = logging.getLogger(__name__)


@dataclass
class StudyConfig:
    """Configuration for convergence study."""

    resolutions: List[int]
    fd_command: Optional[str]
    spectral_command: Optional[str]
    fd_results_pattern: Optional[str]
    spectral_results_pattern: Optional[str]
    fd_output_template: str
    spectral_output_template: str
    output_dir: Path
    reference_path: Optional[Path]
    solver_type: str
    dpi: int
    no_progress: bool
    repo_root: Path


def run_command(command: str, cwd: Path) -> None:
    """Execute shell command with logging."""
    LOGGER.info("Running: %s", command)
    try:
        subprocess.run(
            shlex.split(command),
            cwd=cwd,
            check=True,
            capture_output=False,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        LOGGER.error("Command failed: %s", exc)
        raise


def resolve_result_path(
    resolution: int,
    pattern: Optional[str],
    output_template: str,
    output_dir: Path,
) -> Path:
    """Determine filesystem path of solver output for a resolution."""
    if pattern:
        return Path(pattern.format(nx=resolution, ny=resolution)).expanduser().resolve()
    return (output_dir / output_template.format(nx=resolution, ny=resolution)).resolve()


def obtain_solution(
    resolution: int,
    command_template: Optional[str],
    pattern: Optional[str],
    output_template: str,
    output_dir: Path,
    solver_label: str,
    workdir: Path,
) -> Path:
    """Run solver command if provided or locate existing results."""
    result_path = resolve_result_path(resolution, pattern, output_template, output_dir)
    if result_path.exists():
        LOGGER.info("%s result found: %s", solver_label, result_path)
        return result_path

    if command_template is None:
        raise FileNotFoundError(
            f"No {solver_label} result for resolution {resolution} and no "
            "command template provided to generate it."
        )

    result_path.parent.mkdir(parents=True, exist_ok=True)
    command = command_template.format(nx=resolution, ny=resolution, output=result_path)
    run_command(command, cwd=workdir)

    if not result_path.exists():
        raise FileNotFoundError(
            f"{solver_label} command completed but output not found: {result_path}"
        )
    return result_path


def build_reference_solution(
    config: StudyConfig,
) -> SolutionGrid:
    """Load or compute the reference solution used for error metrics."""
    if config.reference_path:
        LOGGER.info("Loading reference solution: %s", config.reference_path)
        return load_solution(config.reference_path, solver_type="spectral")

    if not config.resolutions:
        raise ValueError("No resolutions provided for reference selection.")

    max_resolution = max(config.resolutions)
    spectral_path = resolve_result_path(
        max_resolution,
        config.spectral_results_pattern,
        config.spectral_output_template,
        config.output_dir,
    )
    if not spectral_path.exists():
        if config.spectral_command:
            command = config.spectral_command.format(
                nx=max_resolution,
                ny=max_resolution,
                output=spectral_path,
            )
            spectral_path.parent.mkdir(parents=True, exist_ok=True)
            run_command(command, cwd=config.repo_root)
        else:
            raise FileNotFoundError(
                "Reference solution not provided and spectral result for the "
                f"largest resolution {max_resolution} not found: {spectral_path}"
            )

    if not spectral_path.exists():
        raise FileNotFoundError(
            f"Spectral reference file expected at {spectral_path} but not found "
            "after attempting to generate it."
        )

    LOGGER.info(
        "Using spectral solution at resolution %d as reference.", max_resolution
    )
    return load_solution(spectral_path, solver_type="spectral")


def perform_study(config: StudyConfig) -> pd.DataFrame:
    """Execute convergence study and return metrics table."""
    reference = build_reference_solution(config)

    records: List[Dict[str, float]] = []
    iterator = config.resolutions
    if not config.no_progress:
        iterator = tqdm(config.resolutions, desc="Resolutions", unit="grid")

    for nx in iterator:
        fd_path = obtain_solution(
            nx,
            config.fd_command,
            config.fd_results_pattern,
            config.fd_output_template,
            config.output_dir,
            "FD",
            config.repo_root,
        )

        fd_solution = load_solution(fd_path, solver_type=config.solver_type)
        mapped = map_solution(fd_solution, reference)
        metrics = compute_errors(reference, mapped)

        grid_spacing = (
            reference.x_coords.max() - reference.x_coords.min()
        ) / (fd_solution.shape[1] - 1)

        records.append(
            {
                "nx": nx,
                "ny": fd_solution.shape[0],
                "grid_spacing": grid_spacing,
                "velocity_l2": metrics["velocity_l2"],
                "velocity_linf": metrics["velocity_linf"],
                "u_l2": metrics["u_l2"],
                "v_l2": metrics["v_l2"],
                "pressure_l2": metrics.get("pressure_l2", np.nan),
            }
        )

    df = pd.DataFrame(records).sort_values("nx")
    if df.empty:
        raise RuntimeError("No convergence data collected.")

    slope, intercept = np.polyfit(
        np.log(df["grid_spacing"]), np.log(df["velocity_l2"]), deg=1
    )
    df["fitted_error"] = np.exp(intercept) * df["grid_spacing"] ** slope
    df.attrs["observed_order"] = -slope

    LOGGER.info("Observed convergence order ~= %.2f", -slope)
    return df


def plot_convergence(df: pd.DataFrame, output_dir: Path, dpi: int) -> Path:
    """Plot log-log convergence of velocity error."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.loglog(df["grid_spacing"], df["velocity_l2"], "o-", label="Computed error")
    ax.loglog(df["grid_spacing"], df["fitted_error"], "--", label="Best fit")
    ax.set_xlabel("Grid spacing h")
    ax.set_ylabel("Relative L2 error in |u|")
    order = df.attrs.get("observed_order", np.nan)
    ax.set_title(f"Convergence Study (observed order ~= {order:.2f})")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend()

    output_path = output_dir / "convergence_loglog.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_error_bars(df: pd.DataFrame, output_dir: Path, dpi: int) -> Path:
    """Plot bar chart of key error metrics."""
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.bar(df["nx"], df["velocity_l2"], width=5, label="Velocity L2")
    ax.set_xlabel("Grid resolution (nx)")
    ax.set_ylabel("Relative L2 error")
    ax.set_title("Error metrics across resolutions")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

    output_path = output_dir / "convergence_errors.png"
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output_path


def save_results(df: pd.DataFrame, output_dir: Path) -> Path:
    """Write metrics table to CSV."""
    table_path = output_dir / "convergence_metrics.csv"
    df.to_csv(table_path, index=False)
    return table_path


def build_parser() -> argparse.ArgumentParser:
    """Construct CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Run or analyse multi-resolution convergence studies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--resolutions",
        type=int,
        nargs="+",
        required=True,
        help="Grid resolutions (nx values) to include.",
    )
    parser.add_argument(
        "--fd-command",
        type=str,
        help="Command template for FD solver with placeholders {nx}, {ny}, {output}.",
    )
    parser.add_argument(
        "--spectral-command",
        type=str,
        help="Command template for spectral solver (optional).",
    )
    parser.add_argument(
        "--fd-results-pattern",
        type=str,
        help="Pattern for existing FD results with placeholders {nx}, {ny}.",
    )
    parser.add_argument(
        "--spectral-results-pattern",
        type=str,
        help="Pattern for spectral results with placeholders {nx}, {ny}.",
    )
    parser.add_argument(
        "--fd-output-template",
        type=str,
        default="fd_result_{nx:04d}.dat",
        help="Filename template (relative to output dir) for FD outputs.",
    )
    parser.add_argument(
        "--spectral-output-template",
        type=str,
        default="spectral_result_{nx:04d}.dat",
        help="Filename template for spectral outputs.",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        help="Optional reference solution (defaults to highest spectral result).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("convergence_study"),
        help="Directory where outputs and generated data are stored.",
    )
    parser.add_argument(
        "--solver-type",
        choices=("auto", "fd", "spectral", "hybrid"),
        default="fd",
        help="Solver hint for FD solutions.",
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
        help="Disable progress bar output.",
    )
    return parser


def main(args: Optional[Sequence[str]] = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    parsed = parser.parse_args(args=args)

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    config = StudyConfig(
        resolutions=parsed.resolutions,
        fd_command=parsed.fd_command,
        spectral_command=parsed.spectral_command,
        fd_results_pattern=parsed.fd_results_pattern,
        spectral_results_pattern=parsed.spectral_results_pattern,
        fd_output_template=parsed.fd_output_template,
        spectral_output_template=parsed.spectral_output_template,
        output_dir=parsed.output_dir.resolve(),
        reference_path=parsed.reference,
        solver_type=parsed.solver_type,
        dpi=parsed.dpi,
        no_progress=parsed.no_progress,
        repo_root=Path(__file__).resolve().parent.parent,
    )

    config.output_dir.mkdir(parents=True, exist_ok=True)
    configure_matplotlib(font_size=11.0, dpi=config.dpi)

    try:
        df = perform_study(config)
    except Exception as exc:  # pylint: disable=broad-except
        LOGGER.error("Convergence study failed: %s", exc)
        return 1

    saved: List[Path] = []
    saved.append(save_results(df, config.output_dir))
    saved.append(plot_convergence(df, config.output_dir, config.dpi))
    saved.append(plot_error_bars(df, config.output_dir, config.dpi))

    LOGGER.info("Outputs generated in %s:", config.output_dir)
    for path in saved:
        LOGGER.info("  - %s", path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
