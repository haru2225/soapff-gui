#!/usr/bin/env python3
"""soapff paper, Fig. 2: how a bare structure file becomes a force field, drawn with real data (the bulk brucite crystal of the library as the
reference, a larger brucite supercell as the new structure; real SOAP vectors and similarities).
    python make_workflow.py -> figs/workflow.png     (run with the soapff Python environment)"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Rectangle
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui")
from soapff_cui.core import read_lammps_data, read_structure, soap_vectors
HERE = os.path.dirname(os.path.abspath(__file__))
EX = "/Users/harutokono/clay/soapff-cui/examples"
plt.rcParams.update({"font.family": "Hiragino Sans", "font.size": 11})
COL = {"mgh": "#3fa34d", "oh": "#f08c3c", "ho": "#8c8c8c", "ow": "#7fb2e0", "hw": "#cfe3f3"}
RAD = {"mgh": 0.55, "oh": 0.42, "ho": 0.25, "ow": 0.42, "hw": 0.25}
ALIAS = {"OW": "ow", "HW": "hw", "OH": "oh", "HO": "ho"}

lib = np.load(os.path.join(EX, "minerals", "library_clayff_bulk.npz"), allow_pickle=False)
X, Q, EL, TY, ORG = lib["X"], lib["q"], lib["element"], np.array([ALIAS.get(t, t) for t in lib["type"]]), lib["origin"]
species = list(lib["species"])

# reference: brucite (001) slab in water (from the reference data file), side view
ref, qref, tref = read_lammps_data(os.path.join(EX, "minerals", "refs_clayff", "brucite_bulk.data"))
tref = [ALIAS.get(t, t) for t in tref]
# new structure: dry brucite crystal, bare file (elements, positions, cell) + the types soapff gave it
CE = "/private/tmp/claude-501/-Users-harutokono-ScoreMD/64bceb4c-942b-429c-8731-e59ed7d791ef/scratchpad/clayff_eval"
new = read_structure(os.path.join(CE, "tests", "brucite_0.extxyz"))
rows = list(csv.DictReader(open(os.path.join(CE, "out", "brucite_0", "charges.csv"))))
tnew = [r["clayff_type"] for r in rows]; qnew = [float(r["charge"]) for r in rows]
Xn = soap_vectors(new, species, 5.0, 4, 3, 0.4)
# the atom shown: an H of the dry crystal near the middle of the slice
pos = new.positions
cand = [i for i in range(len(new)) if new[i].symbol == "H" and 3.0 < pos[i, 1] < 4.0]
a = cand[len(cand) // 2]
# best library entry of each H type -> similarity to the new atom
Hm = np.nonzero(EL == "H")[0]
sims = X[Hm].astype(float) @ Xn[a].astype(float)
groups = {}
for k, i in enumerate(Hm):
    key = (TY[i], round(float(Q[i]), 4), str(ORG[i]))
    if key not in groups or sims[k] > groups[key][0]:
        groups[key] = (sims[k], i)
top = sorted(groups.items(), key=lambda kv: -kv[1][0])[:4]
print("atom", a, [(k, round(v[0], 4)) for k, v in top])


def draw(ax, atoms, types, xlim, zlim, ymax, colour=True, labels=None, hi=None, bonds=True):
    p = atoms.positions
    m = (p[:, 1] < ymax) & (p[:, 0] > xlim[0]) & (p[:, 0] < xlim[1]) & (p[:, 2] > zlim[0]) & (p[:, 2] < zlim[1])
    idx = np.nonzero(m)[0]
    if bonds:
        for i in idx:
            for j in idx:
                if j > i and atoms[i].symbol != atoms[j].symbol and "H" in (atoms[i].symbol, atoms[j].symbol) and "O" in (atoms[i].symbol, atoms[j].symbol):
                    if np.linalg.norm(p[i] - p[j]) < 1.1:
                        ax.plot(p[[i, j], 0], p[[i, j], 2], color="#999", lw=1.2, zorder=1)
                elif j > i and {atoms[i].symbol, atoms[j].symbol} == {"Mg", "O"} and np.linalg.norm(p[i] - p[j]) < 2.3:
                    ax.plot(p[[i, j], 0], p[[i, j], 2], color="#bbb", lw=1.0, zorder=1)
    for i in idx[np.argsort(p[idx, 1])[::-1]]:
        c = COL.get(types[i], "#bbbbbb") if colour else "#c9c9c9"
        r = RAD.get(types[i], 0.4) if colour else (0.55 if atoms[i].symbol == "Mg" else 0.42 if atoms[i].symbol == "O" else 0.25)
        ax.add_patch(Circle((p[i, 0], p[i, 2]), r, fc=c, ec="#444", lw=0.6, zorder=2))
    ax.set_xlim(*xlim); ax.set_ylim(*zlim); ax.set_aspect("equal"); ax.axis("off")
    return idx


fig = plt.figure(figsize=(14, 8.2))
_add = fig.add_axes
def _top_axes(*a, **k):
    ax = _add(*a, **k); ax.set_zorder(5); return ax    # axes above the background cards
fig.add_axes = _top_axes
def card(x, y, w, h, fc, title):
    fig.patches.append(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.002,rounding_size=0.012", transform=fig.transFigure,
                                      fc=fc, ec="#666", lw=1, zorder=0))
    fig.text(x + 0.01, y + h - 0.012, title, ha="left", va="top", fontsize=13, weight="bold")
def arrow(x0, y0, x1, y1, txt=None):
    fig.patches.append(FancyArrowPatch((x0, y0), (x1, y1), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=22, lw=2.2, color="#333", zorder=5))
    if txt:
        fig.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.012, txt, ha="center", va="bottom", fontsize=10, color="#333")

# ---- developer row (once) --------------------------------------------------------------------------------------
fig.text(0.012, 0.995, "開発者（1 回だけ）", fontsize=14, weight="bold", color="#8a3b12", va="top")
card(0.01, 0.60, 0.30, 0.34, "#f4eee4", "① 参照の系（電荷と型が分かっている）")
axr = fig.add_axes([0.012, 0.650, 0.296, 0.235])
draw(axr, ref, tref, (0, 12.5), (1.4, 12.6), 4.0, True)
for typ, c in (("Mg (mgh)  +1.05", "mgh"), ("O (oh)  −0.95", "oh"), ("H (ho)  +0.425", "ho"), ("O (ow)  −0.82", "ow"), ("H (hw)  +0.41", "hw")):
    pass
fig.text(0.02, 0.62, "原子ごとに 型（色）と 電荷", fontsize=10, color="#444")
arrow(0.315, 0.77, 0.365, 0.77, "SOAP")
card(0.37, 0.60, 0.62, 0.34, "#e2ebf5", "② ライブラリ：全原子の 1 行 = （SOAP ベクトル，電荷，型）")
# a few library rows, real SOAP vectors as small heat strips
ex_rows = []
for want in (("mgh", None), ("oh", None), ("ho", None), ("ow", None), ("hw", None)):
    ii = np.nonzero((TY == want[0]) & (ORG == ("montmorillonite_2W.data" if want[0] in ("ow", "hw") else "brucite_bulk.data")))[0]
    ex_rows.append(ii[len(ii) // 3])
fig.text(0.385, 0.87, "型", fontsize=10, weight="bold"); fig.text(0.46, 0.87, "SOAP ベクトル（実際の値，1200 成分）", fontsize=10, weight="bold")
fig.text(0.90, 0.87, "電荷 q", fontsize=10, weight="bold")
for k, i in enumerate(ex_rows):
    y = 0.80 - k * 0.0375
    ax = fig.add_axes([0.46, y, 0.42, 0.030]); v = X[i].astype(float)
    ax.imshow(np.sqrt(np.abs(v[:1200])).reshape(20, 60), aspect="auto", cmap="viridis", vmin=0, vmax=np.sqrt(np.abs(X[:, :1200])).max() * 0.6); ax.axis("off")
    fig.patches.append(Circle((0.395, y + 0.015), 0.0085, transform=fig.transFigure, fc=COL[TY[i]], ec="#333", lw=0.7, zorder=6))
    fig.text(0.41, y + 0.015, TY[i], va="center", fontsize=11)
    fig.text(0.905, y + 0.015, f"{Q[i]:+.4f}", va="center", fontsize=11)
fig.text(0.68, 0.607, "…  5 構造で 3,432 行（Na モンモリロナイト，石英，アルバイト，ブルーサイト，パイロフィライト）", ha="center", fontsize=10, color="#333")

# ---- user row ---------------------------------------------------------------------------------------------------
fig.text(0.012, 0.572, "利用者（構造ファイルを渡すだけ）", fontsize=14, weight="bold", color="#1f5c8a", va="top")
card(0.01, 0.02, 0.215, 0.50, "#f4eee4", "③ 構造ファイル")
axn = fig.add_axes([0.012, 0.215, 0.211, 0.235])
idx = draw(axn, new, tnew, (0, 12.5), (1.4, 12.6), 4.0, False)
axn.add_patch(Circle((pos[a, 0], pos[a, 2]), 0.45, fc="#e8c547", ec="#b3260e", lw=2.2, zorder=4, clip_on=False))
axn.add_patch(Circle((pos[a, 0], pos[a, 2]), 5.0, fc="none", ec="#b3260e", lw=1.4, ls="--", zorder=4))
fig.text(0.02, 0.145, "元素・座標・箱だけ\n（電荷も型もなし）", fontsize=10.5, color="#444", va="top")
fig.text(0.02, 0.085, "赤丸：注目する H\n破線：SOAP のカットオフ 5 Å", fontsize=10, color="#b3260e", va="top")
arrow(0.23, 0.30, 0.265, 0.30)
card(0.27, 0.02, 0.45, 0.50, "#e6f0e4", "④ いちばん似た行を探す（同じ元素の中で）")
axa = fig.add_axes([0.28, 0.380, 0.17, 0.068]); axa.imshow(np.sqrt(np.abs(Xn[a][:1200])).reshape(20, 60), aspect="auto", cmap="viridis", vmin=0, vmax=np.sqrt(np.abs(X[:, :1200])).max() * 0.6); axa.axis("off")
fig.text(0.365, 0.458, "赤丸の H の SOAP ベクトル", ha="center", fontsize=10)
fig.text(0.465, 0.40, "と比べる\n（類似度 =\nコサイン類似度）", fontsize=10, va="center")
for k, ((typ, q, _o), (s, i)) in enumerate(top):
    y = 0.305 - k * 0.072; best = k == 0
    fig.patches.append(FancyBboxPatch((0.28, y - 0.004), 0.43, 0.066, boxstyle="round,pad=0.001,rounding_size=0.008", transform=fig.transFigure,
                                      fc="#fff3c4" if best else "white", ec="#b3260e" if best else "#aaa", lw=2.2 if best else 0.8, zorder=1))
    ax = fig.add_axes([0.29, y + 0.004, 0.17, 0.05]); ax.imshow(np.sqrt(np.abs(X[i][:1200].astype(float))).reshape(20, 60), aspect="auto", cmap="viridis", vmin=0, vmax=np.sqrt(np.abs(X[:, :1200])).max() * 0.6); ax.axis("off")
    fig.text(0.47, y + 0.03, f"{typ}（{ORG[i].replace('pyrophyllite_bulk.data', 'パイロフィライト結晶').replace('brucite_bulk.data', 'ブルーサイト結晶').replace('montmorillonite_2W.data', 'モンモリロナイト')}）", va="center", fontsize=10)
    fig.text(0.47, y + 0.008, f"q = {q:+.4f}", va="center", fontsize=10)
    fig.text(0.70, y + 0.03, f"{s:.4f}", ha="right", va="center", fontsize=13, weight="bold", color="#b3260e" if best else "#333")
fig.text(0.70, 0.374, "類似度", ha="right", fontsize=10, color="#333")
fig.text(0.495, 0.045, "最大の行の 電荷 と 型 をそのまま使う（0.95 未満なら ClayFF の型に戻す）", ha="center", fontsize=10, color="#1b5e20", weight="bold")
arrow(0.725, 0.30, 0.76, 0.30)
card(0.765, 0.02, 0.225, 0.50, "#fbe8dc", "⑤ LAMMPS の入力")
axo = fig.add_axes([0.767, 0.215, 0.221, 0.235])
draw(axo, new, tnew, (0, 12.5), (1.4, 12.6), 4.0, True)
fig.text(0.775, 0.17, "全原子に 電荷 と 型\n結合・角度は型から決まる", fontsize=10.5, color="#444", va="top")
fig.text(0.775, 0.105, "system.data  system.settings\nin.lammps  charges.csv", fontsize=10, color="#333", va="top", family="monospace")
fig.savefig(os.path.join(HERE, "figs", "workflow.png"), dpi=170); plt.close(fig)
print("ok")
