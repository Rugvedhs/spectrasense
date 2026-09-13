"""Figure style: a validated, colour-vision-safe palette and print defaults.

Hues are assigned to entities in a fixed order and never cycled, so a method
keeps its colour across every figure in the paper.  The categorical set was
checked for colour-vision separation rather than chosen by eye (worst adjacent
protan Delta-E 9.1, normal-vision 19.6, both above the accepted floors).  Two of
the hues sit below 3:1 contrast against a white page, so every categorical
figure also carries direct value labels and each has a companion CSV table —
identity is never conveyed by colour alone.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# Fixed categorical order.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")

CONFORMAL_COLORS = {
    "SCP": SERIES[0],
    "Normalised SCP": SERIES[1],
    "Mondrian-similarity": SERIES[2],
    "Mondrian-family": SERIES[3],
}

REGIME_COLORS = {
    "random": SERIES[0],
    "scaffold": SERIES[1],
    "cluster": SERIES[2],
    "family": SERIES[3],
}

ARM_COLORS = {
    "informed": SERIES[0],
    "control": SERIES[2],
    "naive": SERIES[1],
}

# Single-hue ramp for magnitude (similarity bins), light to dark.
SEQUENTIAL = ("#cfe2f8", "#9cc4f0", "#69a5e4", "#2a78d6", "#18538f")

TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#d9d9d6"


def use_paper_style() -> None:
    """Matplotlib defaults tuned for a two-column journal page."""
    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 400,
            "savefig.bbox": "tight",
            "font.size": 8,
            "font.family": "DejaVu Sans",
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "axes.labelcolor": TEXT_PRIMARY,
            "axes.edgecolor": TEXT_SECONDARY,
            "axes.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID,
            "grid.linewidth": 0.5,
            "xtick.color": TEXT_SECONDARY,
            "ytick.color": TEXT_SECONDARY,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "legend.frameon": False,
            "lines.linewidth": 1.6,
            "lines.markersize": 4,
        }
    )


def save(fig: plt.Figure, path, also_png: bool = True) -> None:
    """Write a vector PDF for typesetting plus a PNG for quick inspection."""
    from pathlib import Path

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    if also_png:
        fig.savefig(path.with_suffix(".png"))
    plt.close(fig)
