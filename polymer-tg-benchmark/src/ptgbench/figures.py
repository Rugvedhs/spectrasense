"""Figure style: journal conventions, not dashboard conventions.

The rules enforced here are the ones that separate a published figure from a
slide:

* **No titles inside the axes.** The caption carries the description; panels are
  identified by a bold letter. A sentence sitting above the plot area is a
  presentation habit.
* **Serif type matching the body text**, at the size it will actually be printed,
  so the figure reads as part of the page rather than as an insert.
* **Thin marks and hairline spines.** Chunky bars and heavy lines are the most
  reliable visual tell of an unedited default.
* **Selective labelling.** A number on every mark is clutter; label the ones that
  carry the argument and let the companion table hold the rest.
* **Direct labels over legend boxes** wherever a line ends in free space.

The palette is Okabe–Ito, the standard colour-vision-safe qualitative set for
scientific figures. It was checked with a validator rather than by eye, and
passes on both the adjacent and the all-pairs criterion (worst all-pairs deuteran
ΔE 11.0, normal-vision ΔE 15.6), so it is safe for scatter and small-multiple
forms as well as bars and lines. Okabe–Ito orange sits below 3:1 contrast on
white, so every figure using it also carries direct labels or a companion table.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe-Ito, fixed order. A series keeps its hue across every figure.
BLUE = "#0072B2"
VERMILLION = "#D55E00"
GREEN = "#009E73"
ORANGE = "#E69F00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"

SERIES = (BLUE, VERMILLION, GREEN, ORANGE, PURPLE, SKY)

CONFORMAL_COLORS = {
    "SCP": BLUE,
    "Normalised SCP": VERMILLION,
    "Mondrian-family": ORANGE,
    "Mondrian-similarity": GREEN,
}

REGIME_COLORS = {
    "random": BLUE,
    "scaffold": VERMILLION,
    "cluster": ORANGE,
    "family": GREEN,
}

ARM_COLORS = {"informed": BLUE, "control": GREEN, "naive": VERMILLION}

MODEL_COLORS = {
    "hist_gbr": BLUE,
    "extra_trees": VERMILLION,
    "random_forest": GREEN,
    "svr": ORANGE,
    "median": "#8c8c8c",
}

INK = "#1a1a1a"
MUTED = "#5c5c5c"
HAIRLINE = "#9a9a9a"
GRID = "#e2e2e2"
FILL = "#d9d9d9"
REFERENCE = "#8c8c8c"

# Printed widths. Journals set figures to a column measure; drawing at that size
# means the type in the figure matches the type on the page.
SINGLE_COLUMN = 3.46   # 88 mm
ONE_HALF_COLUMN = 5.5  # 140 mm
DOUBLE_COLUMN = 7.09   # 180 mm


def use_paper_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "font.family": "serif",
            "font.serif": ["DejaVu Serif"],
            "font.size": 7.5,
            "axes.titlesize": 7.5,
            "axes.labelsize": 7.5,
            "axes.labelcolor": INK,
            "axes.edgecolor": HAIRLINE,
            "axes.linewidth": 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "axes.axisbelow": True,
            "grid.color": GRID,
            "grid.linewidth": 0.4,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.fontsize": 6.8,
            "legend.frameon": False,
            "legend.handlelength": 1.4,
            "legend.handletextpad": 0.5,
            "legend.labelspacing": 0.35,
            "legend.borderaxespad": 0.2,
            "lines.linewidth": 1.1,
            "lines.markersize": 3.2,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel(ax, letter: str, dx: float = -0.02, dy: float = 1.04) -> None:
    """Bold panel identifier, placed outside the plot area."""
    ax.text(
        dx, dy, f"({letter})", transform=ax.transAxes, fontsize=8,
        fontweight="bold", color=INK, ha="left", va="bottom",
    )


def yardstick(ax, axis: str = "y") -> None:
    """A single faint grid direction, for reading values off long axes."""
    ax.grid(True, axis=axis, color=GRID, linewidth=0.4)
    ax.set_axisbelow(True)


def reference_line(ax, value: float, label: str | None = None,
                   orientation: str = "h", where: float = 0.012) -> None:
    """A dashed neutral rule for a nominal level, labelled in the margin."""
    if orientation == "h":
        ax.axhline(value, color=REFERENCE, linewidth=0.7, linestyle=(0, (4, 3)),
                   zorder=1)
        if label:
            ax.text(where, value, label, transform=ax.get_yaxis_transform(),
                    fontsize=6, color=REFERENCE, ha="left", va="bottom")
    else:
        ax.axvline(value, color=REFERENCE, linewidth=0.7, linestyle=(0, (4, 3)),
                   zorder=1)


def end_label(ax, x, y, text: str, color: str, dx: float = 0.12,
              fontsize: float = 6.6) -> None:
    """Name a line where it ends, so the reader never consults a legend box."""
    ax.text(x + dx, y, text, color=color, fontsize=fontsize, va="center",
            ha="left", fontweight="medium")


def spread_labels(values, minimum_gap: float):
    """Nudge end-label positions apart while keeping their order.

    Two series that finish within a hair of each other would otherwise print
    their names on top of one another. This moves the labels, never the data.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    placed = list(values)
    for rank, index in enumerate(order):
        if rank == 0:
            continue
        previous = placed[order[rank - 1]]
        if placed[index] - previous < minimum_gap:
            placed[index] = previous + minimum_gap
    return placed


def save(fig: plt.Figure, path) -> None:
    """Vector PDF for typesetting, PNG for inspection."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"))
    plt.close(fig)
