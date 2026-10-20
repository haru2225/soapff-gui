#!/usr/bin/env python3
"""soapff paper: OVITO rendering of the whole test system (clay, water, Na), coloured by the ClayFF type soapff assigned.
Uses the helpers of ../../test8/ovito_figs.py.  Run with the OVITO Python environment."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "/Users/harutokono/clay/test8")
import ovito_figs as of
from ovito.modifiers import SliceModifier
HERE = os.path.dirname(os.path.abspath(__file__))
W, H = 1900, 1500
vp = of.viewport(15.5, 12.7, 14.3)
ymid = float(np.median(of.pos[:, 1]))
img = of.render(of.pipeline([SliceModifier(normal=(0, 1, 0), distance=ymid, slab_width=8.0)]), vp, (W, H), "system_full.png")
ext = of.calibrate(vp, (W, H), of.pick_two((np.abs(of.pos[:, 1] - ymid) < 4.0) & ~np.isin(of.lab, ["OW", "HW", "Na"])))
fig, ax = plt.subplots(figsize=(11, 8.5))
ax.imshow(img, extent=ext, aspect="equal")
ax.set_xlim(-1.5, 32.5); ax.set_ylim(-1.5, 26.5); ax.set_xlabel("x (A)"); ax.set_ylabel("z (A)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=of.COL[t], markeredgecolor="#333", markersize=8, label=of.NAME[t])
                   for t in ("st", "ao", "mgo", "ob", "obos", "oh", "ohs", "ho", "Na", "OW", "HW")],
          fontsize=9, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", "system.png"), dpi=160, bbox_inches="tight")
print("ok")
