# SPDX-License-Identifier: GPL-2.0-only
# Copyright (c) 2026 Muhammad Uzair
"""One style for every data figure of the TAES paper (taken over from the T-TST paper).

Fonts: Times-like serif at final size (8 pt labels, 7 pt ticks and legends), STIX math,
TrueType embedding. Palette: the CVD-checked categorical set of the ieee-transactions
figures.md, in fixed order and never cycled. Every method keeps one color, one marker and
one dash pattern in every figure, so the plots also read in grayscale.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GREEN, GOLD = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
PURPLE, TEAL, RED, GREY = "#8a5cd1", "#0f9ba8", "#d64545", "#7b8794"
INK, INK2 = "#1f2933", "#52606d"
COL, PAGE = 3.5, 7.16

# method identity: (color, marker, dash)
DASH = {"solid": "-", "dash": (0, (4, 1.6)), "dot": (0, (1, 1.2)), "dashdot": (0, (5, 1.5, 1, 1.5)),
        "longdash": (0, (7, 2)), "dashdotdot": (0, (4, 1.2, 1, 1.2, 1, 1.2))}
METHOD = {
    "SnapFast": (INK, "s", DASH["dash"]),          # snapshot, 100 ms (Hypatia)
    "SnapSlow": (GREY, "s", DASH["dot"]),          # snapshot, 1 s
    "OdZero": (PURPLE, "^", DASH["dashdot"]),      # on-demand, no setup delay
    "OdOne": (RED, "v", DASH["dot"]),              # on-demand, 1 s setup
    "Asr": (GOLD, "D", DASH["longdash"]),          # merged snapshots, 30 s
    "LifeOne": (GREEN, "o", DASH["solid"]),        # lifetime-aware, eps 0.01
    "LifeFive": (BLUE, "o", DASH["solid"]),        # lifetime-aware, eps 0.05
    "LifeTen": (ORANGE, "o", DASH["solid"]),       # lifetime-aware, eps 0.10
}
THRESH = dict(color=RED, lw=0.7, ls=(0, (1, 1.5)))


def setup():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Nimbus Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "legend.frameon": False,
        "legend.handlelength": 2.2, "legend.columnspacing": 1.0, "legend.handletextpad": 0.4,
        "axes.linewidth": 0.7, "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "axes.grid": True, "grid.color": "#c4ccd6",
        "grid.linestyle": ":", "grid.linewidth": 0.5, "axes.axisbelow": True,
        "lines.linewidth": 1.0, "lines.markersize": 3.8, "lines.markeredgewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    })


def style(key, hollow=False):
    c, m, ls = METHOD[key]
    return dict(color=c, marker=m, ls=ls, mfc="white" if hollow else c, mec=c)
