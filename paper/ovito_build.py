#!/usr/bin/env python3
"""soapff paper: the four stages of building the quartz (101)-water interface (OVITO, side view).
Input: extxyz files of the stages (a cut slab, b hydroxylated, c water placed, d relaxed) written from build_ase.py;
added OH/H atoms (array 'added') drawn larger, water (index >= n_slab) in light colours.  Run with the OVITO Python."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from ovito.io import import_file
from ovito.modifiers import AssignColorModifier, ComputePropertyModifier, CreateBondsModifier, ExpressionSelectionModifier
from ovito.vis import TachyonRenderer, Viewport
HERE = os.path.dirname(os.path.abspath(__file__))
D = sys.argv[1]
N_SLAB = 252
EL = {"Si": ("#e8c547", 0.42), "O": ("#d9412b", 0.36), "H": ("#8c8c8c", 0.2)}
TITLES = {"a": "(a) 切り出したスラブ", "b": "(b) OH と H で終端", "c": "(c) 水を配置", "d": "(d) 緩和後（参照）"}
plt.rcParams.update({"font.family": "Hiragino Sans"})
fig, axs = plt.subplots(1, 4, figsize=(16, 9))
for ax, k in zip(axs, "abcd"):
    p = import_file(os.path.join(D, f"{k}.xyz"))
    data = p.source.data
    data.cell_.vis.enabled = False
    for t in data.particles_.particle_types_.types_:
        if t.name in EL:
            t.color = to_rgb(EL[t.name][0]); t.radius = EL[t.name][1]
    tid = {t.name: t.id for t in data.particles.particle_types.types}
    if k in "cd":
        p.modifiers.append(ExpressionSelectionModifier(expression=f"ParticleIndex >= {N_SLAB} && ParticleType == {tid['O']}"))
        p.modifiers.append(AssignColorModifier(color=to_rgb("#7fb2e0")))
        p.modifiers.append(ExpressionSelectionModifier(expression=f"ParticleIndex >= {N_SLAB} && ParticleType == {tid['H']}"))
        p.modifiers.append(AssignColorModifier(color=to_rgb("#dcebf7")))
    if k == "b":
        p.modifiers.append(ComputePropertyModifier(output_property="Radius", expressions=[f"added ? (ParticleType=={tid['H']} ? 0.32 : 0.55) : 0"]))
    b = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
    b.set_pairwise_cutoff("Si", "O", 2.0); b.set_pairwise_cutoff("O", "H", 1.15)
    b.vis.width = 0.2
    p.modifiers.append(b)
    p.add_to_scene()
    cell = np.array(data.cell[:, :3]); c = cell.sum(0) / 2
    vp = Viewport(type=Viewport.Type.Ortho, camera_pos=(c[0], -80, c[2]), camera_dir=(0, 1, 0), fov=float(cell[2, 2]) / 2 + 1)
    path = os.path.join(HERE, "figs", f"_build_{k}.png")
    vp.render_image(size=(700, 1700), filename=path, background=(1, 1, 1), renderer=TachyonRenderer(ambient_occlusion=True))
    p.remove_from_scene()
    ax.imshow(plt.imread(path)); ax.axis("off"); ax.set_title(TITLES[k], loc="left", fontsize=14)
    os.remove(path)
h = [Line2D([], [], marker="o", ls="", color=c, markeredgecolor="#333", markersize=9, label=l) for c, l in
     (("#e8c547", "Si"), ("#d9412b", "O（石英）"), ("#8c8c8c", "H（OH）"), ("#7fb2e0", "O（水）"), ("#dcebf7", "H（水）"))]
fig.legend(handles=h, loc="lower center", ncol=5, fontsize=12, frameon=False)
fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.06, wspace=0.03); fig.savefig(os.path.join(HERE, "figs", "build.png"), dpi=150)
print("ok")
