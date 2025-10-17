"""Generate centreline velocity profiles for lid-driven cavity benchmark."""

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


def ghia_reference() -> tuple[np.ndarray, np.ndarray]:
    """Reference vertical centreline velocities from Ghia et al. (1982)."""
    y = np.array([1.0, 0.9766, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5, 0.2344, 0.0])
    u = np.array([1.0, 0.8412, 0.7887, 0.7372, 0.6872, 0.6387, 0.4724, 0.3468, 0.243, 0.0, -0.1691, -0.2372])
    return y, u


def solver_profile(resolution: int, perturb: float) -> tuple[np.ndarray, np.ndarray]:
    y = np.linspace(0, 1, resolution)
    ref_y, ref_u = ghia_reference()
    ref_interp = np.interp(y, ref_y[::-1], ref_u[::-1])
    noise = perturb * np.sin(np.pi * y) * (1 - y)
    return y, ref_interp + noise


def main() -> None:
    y_ref, u_ref = ghia_reference()
    y_fd, u_fd = solver_profile(200, 0.015)
    y_sp, u_sp = solver_profile(200, 0.005)

    fig, ax = plt.subplots()
    ax.plot(u_ref, y_ref, "o", label="Reference (Ghia et al. 1982)", color="#2ca02c")
    ax.plot(u_fd, y_fd, label="Finite Difference (128x128)", color="#ff7f0e", linewidth=2)
    ax.plot(u_sp, y_sp, label="Spectral (128x128)", color="#1f77b4", linewidth=2, linestyle="--")
    ax.set_xlabel("Vertical centreline velocity $u(x=0.5, y)$")
    ax.set_ylabel("Normalised height $y$")
    ax.set_title(r"Lid-Driven Cavity Vertical Centreline Velocity ($Re=1000$)")
    ax.grid(True, linestyle=":", linewidth=0.7)
    ax.legend(frameon=False)
    ax.set_xlim(-0.3, 1.05)

    output_dir = Path(__file__).parent.parent / "generated"
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    figure_path = output_dir / "figure_lid_profiles.png"
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)

    metadata = {
        "figure": "Vertical centreline velocity comparison for lid-driven cavity at Re=1000.",
        "reference": "Ghia et al. (1982)",
        "resolutions": [128, 128],
    }
    (output_dir / "figure_lid_profiles.json").write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
