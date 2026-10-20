#!/usr/bin/env python3
"""soapff paper figures: workflow, similarity histograms of the held-out frames, MD time series."""
import json, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "figs")
plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})

# 1. workflow
fig, ax = plt.subplots(figsize=(13, 4.2)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 4.2)
def box(x, y, w, h, title, body, c):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.15", fc=c, ec="#555", lw=1))
    ax.text(x + w / 2, y + h - 0.28, title, ha="center", va="top", fontsize=12, weight="bold")
    ax.text(x + w / 2, y + h - 0.75, body, ha="center", va="top", fontsize=10, linespacing=1.4)
def arrow(x0, y0, x1, y1):
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color="#333", lw=1.6))
box(0.1, 2.3, 3.0, 1.8, "reference system", "initial structure:\nLAMMPS data with charges\nand ClayFF type names", "#ECE8DF")
box(0.1, 0.1, 3.0, 1.8, "new structure", "elements, positions, cell\n(xyz, cif, POSCAR, data)\nno charges, no types", "#ECE8DF")
box(3.9, 2.3, 2.8, 1.8, "library", "per element:\nSOAP vectors X\n+ charges q + types", "#dbe6f3")
box(3.9, 0.1, 2.8, 1.8, "SOAP of every atom", "clay, water, Na\nr_cut 5 A, l_max 3\nn_max 4, sigma 0.4 A", "#dbe6f3")
box(7.5, 0.6, 2.6, 3.0, "nearest entry", "S = F X^T\nj* = argmax_j S_ij\ncharge = q_j*\ntype = type_j*\n(similarity < 0.95:\nClayFF fallback)", "#f3dfd5")
box(10.8, 0.6, 2.1, 3.0, "LAMMPS input", "system.data\n(charges, bonds,\nangles from types)\nsystem.settings\nin.lammps\ncharges.csv", "#e3efe1")
arrow(3.1, 3.2, 3.9, 3.2); arrow(3.1, 1.0, 3.9, 1.0); arrow(6.7, 3.2, 7.5, 2.6); arrow(6.7, 1.0, 7.5, 1.6); arrow(10.1, 2.1, 10.8, 2.1)
fig.tight_layout(); fig.savefig(os.path.join(F, "workflow.png"), dpi=170); plt.close(fig)

# 3. MD time series
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
for s, c, lab in (("initial", "#A8461F", "force field generated from the initial structure"),):
    txt = open(os.path.join(HERE, f"md_{s}", "log.lammps")).read()
    blk = txt.split("fix md all nvt")[1].split("Loop time")[0]
    rows = [l.split() for l in blk.splitlines() if re.match(r"^\s+\d+\s+[-\d.]", l)]
    a = np.array([[float(x) for x in r[:4]] for r in rows])
    axs[0].plot(a[:, 0] * 1e-3, a[:, 1], color=c, lw=1.5, label=lab)
    axs[1].plot(a[:, 0] * 1e-3, a[:, 3] / 2502, color=c, lw=1.5, label=lab)
axs[0].axhline(300, color="#5C676C", ls=":", lw=1); axs[0].set_ylabel("temperature (K)")
axs[1].set_ylabel("potential energy per atom (kcal/mol)")
for ax in axs:
    ax.set_xlabel("time (ps)"); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(F, "md.png"), dpi=170); plt.close(fig)
print("ok")
