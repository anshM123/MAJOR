"""Shared matplotlib style for the figures of the Letter."""
import matplotlib as mpl

BLUE, ORANGE, AQUA, GRAY = "#257de0", "#f47230", "#00b894", "#8a8a8a"
INK, INK2, GRID = "#191919", "#666666", "#d9d9d9"
RED, VIOLET, GREEN = "#ec4647", "#6a4ce2", "#008a17"

def apply():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.0,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
    })
