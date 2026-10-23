# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "scipy", "matplotlib"]
# ///
"""Render the six-panel reproductions of Barwicz Figs. 11 and 12.

Reads output/example3/stage3_meshes.mat, written by
examples/example3_reproduce_fig11_fig12.m, and writes the two composite PNGs
next to it.  This is the renderer that produced the images in figures/.

Layout mirrors the paper: rows = SiON (a,b), Si3N4 (c,d), Si (e,f);
columns = TE-like, TM-like.  Contour levels are Tymon Barwicz's own, taken
from his plotting scripts.

Usage (from the repository root):

    uv run examples/plot_fig11_fig12.py
    # or, with numpy/scipy/matplotlib already installed:
    python examples/plot_fig11_fig12.py
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np
from scipy.io import loadmat

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MESHES = ROOT / "output" / "example3" / "stage3_meshes.mat"
OUT = MESHES.parent

lv = np.arange

# name, axis limits (nm), Fig-12 TE/TM levels, Fig-11 TE/TM levels, title, panel tags
SYSTEMS = [
    (
        "SiON", (750, 4000),
        lv(0, 0.00111, 0.00002), lv(0, 0.00111, 0.00004),
        np.concatenate([lv(0.20, 0.401, 0.05), lv(0.41, 0.481, 0.01),
                        lv(0.50, 0.751, 0.05), lv(0.78, 1.201, 0.02)]),
        np.concatenate([lv(0.20, 0.401, 0.05), lv(0.41, 0.481, 0.01),
                        lv(0.50, 0.751, 0.05), lv(0.78, 1.201, 0.02)]),
        "SiON/SiO$_2$ ($n_{core}$=1.50)", ("(a)", "(b)"),
    ),
    (
        "Si3N4", (200, 1000),
        np.concatenate([[0, 0.005], lv(0.01, 0.231, 0.01)]),
        np.concatenate([[0, 0.005], lv(0.01, 0.231, 0.01)]),
        lv(0.10, 0.701, 0.05), lv(0.10, 0.701, 0.05),
        "Si$_3$N$_4$/SiO$_2$ ($n_{core}$=2.00)", ("(c)", "(d)"),
    ),
    (
        "Si", (120, 600),
        np.concatenate([lv(0, 0.81, 0.2), lv(1, 6.51, 0.5), [7, 12]]),
        np.concatenate([[0, 0.02, 0.05, 0.1], lv(0.2, 0.81, 0.2),
                        lv(1, 6.51, 0.5), [7, 12]]),
        np.concatenate([lv(0, 0.051, 0.01), lv(0.06, 0.801, 0.02)]),
        np.concatenate([lv(0, 0.061, 0.01), [0.08], lv(0.10, 0.801, 0.02)]),
        "Si/SiO$_2$ ($n_{core}$=3.50)", ("(e)", "(f)"),
    ),
]

FIGURES = [
    ("te_factor", "tm_factor",
     "Radiation efficiency P$_{3D}$/P$_{2D}$ — reproduction of Barwicz Fig. 11",
     "fig11_reproduction.png"),
    ("te_total", "tm_total",
     "Scattering losses normalized to roughness variance (dB/cm)/nm$^2$"
     " — reproduction of Barwicz Fig. 12",
     "fig12_reproduction.png"),
]


def main() -> int:
    if not MESHES.exists():
        print(f"missing {MESHES}", file=sys.stderr)
        print("Run examples/example3_reproduce_fig11_fig12.m first.", file=sys.stderr)
        return 1

    M = loadmat(MESHES)

    for key_te, key_tm, suptitle, fname in FIGURES:
        fig, axes = plt.subplots(3, 2, figsize=(12, 14))
        for r, (name, lims, l12te, l12tm, l11te, l11tm, title, tags) in enumerate(SYSTEMS):
            b = M[f"{name}_b"] * 1e9
            a = M[f"{name}_a"] * 1e9
            levels = (l11te, l11tm) if key_te.endswith("factor") else (l12te, l12tm)
            for c, (key, pol) in enumerate([(key_te, "TE-like"), (key_tm, "TM-like")]):
                ax = axes[r, c]
                z = M[f"{name}_{key}"]
                ax.contourf(b, a, z, levels=levels[c], cmap="gray_r", extend="max")
                cs = ax.contour(b, a, z, levels=levels[c], colors="k", linewidths=0.4)
                ax.clabel(cs, inline=True, fontsize=6, fmt="%g")
                ax.set_xlim(*lims)
                ax.set_ylim(*lims)
                ax.set_title(f"{tags[c]} {title} — {pol}", fontsize=10)
                ax.set_xlabel("Width of the core (nm)", fontsize=9)
                ax.set_ylabel("Height of the core (nm)", fontsize=9)
        fig.suptitle(suptitle, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(OUT / fname, dpi=150)
        plt.close(fig)
        print(f"wrote {OUT / fname}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
