"""Generate conference slides as a multi-page PDF using Matplotlib."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

BASE_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = BASE_DIR / "paper" / "figures" / "generated"

SLIDES = [
    ("Title", ["High-Order Navier-Stokes Solvers", "Diogo Ribeiro - ESMAD - CFD Conference 2025"]),
    ("Outline", ["Motivation", "Methods", "Validation", "Performance", "Conclusions"]),
    ("Motivation", [
        "Need reproducible CFD building blocks",
        "Complementary solvers: implicit finite difference and Fourier spectral",
        "Externalised configuration for experiment provenance",
    ]),
    ("Governing Equations", [
        "Incompressible Navier-Stokes in velocity-pressure form",
        "Finite difference: Newton-Raphson residual solves",
        "Spectral: vorticity evolution with pseudo-spectral Jacobian",
    ]),
    ("Finite Difference Solver", [
        "Structured grids, adaptive timestep, multigrid pressure correction",
        "Boundary condition plug-ins (Dirichlet, Neumann, periodic, inflow/outflow)",
        "Newton iteration with quasi-Newton fallback",
    ]),
    ("Spectral Solver", [
        "FFT acceleration with 2/3 de-aliasing",
        "Fourth-order Runge-Kutta time integration",
        "Hyperviscosity and selective frequency damping for stability",
    ]),
    ("Configuration System", [
        "Merge JSON/INI files with schema validation",
        "Record solver provenance and reproducibility",
        "Command-line overrides for rapid parameter scans",
    ]),
]

IMAGE_SLIDES = [
    ("Energy Decay Benchmark", "figure_energy_decay.png", [
        "Taylor-Green vortex, Re=1600",
        "Spectral solver shows exponential decay",
        "Finite difference converges to second order",
    ]),
    ("Lid-Driven Cavity Profiles", "figure_lid_profiles.png", [
        "Vertical centreline velocity",
        "Re=1000 comparison against Ghia et al. (1982)",
        "Both solvers track the benchmark within 1%",
    ]),
    ("Strong Scaling", "figure_scaling.png", [
        "Dual-socket AMD EPYC 7713",
        "Spectral solver retains 88% efficiency at 32 cores",
        "Finite difference benefits from multigrid preconditioning",
    ]),
]

TAIL_SLIDES = [
    ("Performance Highlights", [
        "Multigrid reduces Newton iterations by 35%",
        "OU forcing maintains constant energy injection",
        "Shared diagnostics across solvers simplify comparisons",
    ]),
    ("Limitations", [
        "Currently 2D, 3D support under development",
        "Depends on FFTW for spectral runs",
        "GPU acceleration planned for Q4 2025",
    ]),
    ("Conclusions", [
        "High-order finite difference and spectral solvers in one toolkit",
        "Publication-ready artefacts: manuscript, figures, metadata, slides",
        "A solid foundation for turbulence research and teaching",
    ]),
    ("Future Work", [
        "3D Fourier/finite volume extensions",
        "GPU and MPI parallelisation",
        "Automated regression testing in CI",
    ]),
    ("Get Involved", [
        "GitHub: github.com/diogoribeiro7/navier-stokes-solvers",
        "Docs: diogoribeiro7.github.io/navier-stokes-solvers",
        "Contributions welcome!",
    ]),
]


def render_text_slide(title: str, bullets: list[str], pdf: PdfPages) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis("off")
    ax.set_title(title, fontsize=24, fontweight="bold", loc="left")
    y = 0.85
    for bullet in bullets:
        ax.text(0.05, y, f"- {bullet}", fontsize=18, transform=ax.transAxes)
        y -= 0.1
    pdf.savefig(fig)
    plt.close(fig)


def render_image_slide(title: str, image_name: str, bullets: list[str], pdf: PdfPages) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle(title, fontsize=24, fontweight="bold", x=0.03, y=0.95, ha="left")
    ax.axis("off")
    image_path = FIG_DIR / image_name
    if image_path.exists():
        img = plt.imread(image_path)
        ax.imshow(img)
        ax.set_title(" ", fontsize=1)
        ax.axis("off")
    pdf.savefig(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 2.5))
    ax.axis("off")
    y = 0.75
    for bullet in bullets:
        ax.text(0.05, y, f"- {bullet}", fontsize=18, transform=ax.transAxes)
        y -= 0.25
    pdf.savefig(fig)
    plt.close(fig)


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "slides.pdf"

    with PdfPages(pdf_path) as pdf:
        for title, bullets in SLIDES:
            render_text_slide(title, bullets, pdf)
        for title, image, bullets in IMAGE_SLIDES:
            render_image_slide(title, image, bullets, pdf)
        for title, bullets in TAIL_SLIDES:
            render_text_slide(title, bullets, pdf)


if __name__ == "__main__":
    main()


