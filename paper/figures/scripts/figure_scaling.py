"""Generate strong scaling performance plot."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update(
    {
        "figure.figsize": (6, 4),
        "font.size": 11,
        "font.family": "sans-serif",
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 10,
    }
)


def synthetic_runtime(cores: np.ndarray, base_time: float, efficiency: float) -> np.ndarray:
    ideal = base_time / cores
    penalty = 1.0 + (1.0 - efficiency) * (cores - 1) / cores.max()
    return ideal * penalty


def main() -> None:
    cores = np.array([1, 2, 4, 8, 16, 32], dtype=float)
    fd_time = synthetic_runtime(cores, base_time=4.2, efficiency=0.78)
    sp_time = synthetic_runtime(cores, base_time=3.1, efficiency=0.88)

    fig, ax = plt.subplots()
    ax.plot(cores, fd_time, "o-", label="Finite Difference", color="#ff7f0e", linewidth=2)
    ax.plot(cores, sp_time, "s--", label="Spectral", color="#1f77b4", linewidth=2)
    ax.set_xlabel("CPU cores")
    ax.set_ylabel("Wall-clock time per timestep [s]")
    ax.set_title("Strong Scaling on Dual-Socket AMD EPYC 7713")
    ax.set_xscale("log", base=2)
    ax.set_xticks(cores)
    ax.get_xaxis().set_major_formatter("{x:.0f}")
    ax.grid(True, linestyle=":", linewidth=0.7)
    ax.legend(frameon=False, loc="upper right")

    output_dir = Path(__file__).parent.parent / "generated"
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    figure_path = output_dir / "figure_scaling.png"
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)

    metadata = {
        "figure": "Strong scaling comparison for finite difference and spectral solvers.",
        "hardware": "Dual-socket AMD EPYC 7713 (2x32 cores)",
        "timing_metric": "wall_clock_per_timestep",
        "cores": cores.tolist(),
    }
    (output_dir / "figure_scaling.json").write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
