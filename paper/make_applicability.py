#!/usr/bin/env python3
"""soapff paper: applicability range.  Highest SOAP similarity of every atom of bare dry bulk crystals to the library
(4 systems), per element, with the 0.95 threshold.   python make_applicability.py -> figs/applicability.png"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
A = "/Users/harutokono/clay/soapff-cui/examples/applicability"
plt.rcParams.update({"font.family": "Hiragino Sans", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
MIN = [("quartz", "石英"), ("albite", "アルバイト"), ("brucite", "ブルーサイト"), ("pyrophyllite", "パイロフィライト"),
       ("kaolinite", "カオリナイト"), ("microcline", "マイクロクリン"), ("forsterite", "フォルステライト")]
EC = {"Si": "#e0b228", "Al": "#b06ad0", "Mg": "#3fa34d", "O": "#d9412b", "H": "#7d7d7d", "Na": "#3f6fd8", "K": "#7a4a2a"}
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.axhspan(0.95, 1.01, color="#e3efe1", zorder=0); ax.axhspan(0.70, 0.95, color="#fbe8dc", zorder=0)
rng = np.random.default_rng(0)
seen = set()
for k, (m, jp) in enumerate(MIN):
    rows = list(csv.DictReader(open(os.path.join(A, f"charges_{m}.csv"))))
    els = sorted({r["element"] for r in rows}, key=lambda e: list(EC).index(e))
    for n, e in enumerate(els):
        s = np.array([float(r["similarity"]) for r in rows if r["element"] == e])
        s = np.where(s < 0.70, 0.705, s)        # K (no K in the library, similarity 0) drawn at the bottom
        x = k + (n - (len(els) - 1) / 2) * 0.14 + rng.uniform(-0.04, 0.04, len(s))
        ax.scatter(x, s, s=9, color=EC[e], alpha=0.7, lw=0, label=e if e not in seen else None, zorder=3); seen.add(e)
ax.axhline(0.95, color="#333", ls="--", lw=1.2); ax.text(6.55, 0.952, "しきい値 0.95", ha="right", va="bottom", fontsize=10)
ax.text(-0.45, 0.995, "SOAP で決まる（ライブラリの電荷と型）", fontsize=10, color="#1b5e20", va="top")
ax.text(-0.45, 0.715, "ClayFF の型に戻す（距離による型付け）", fontsize=10, color="#8a3b12", va="bottom")
ax.text(5.0, 0.712, "K：ライブラリにない（類似度 0）", fontsize=9, color="#7a4a2a", ha="center", va="bottom")
ax.set_xticks(range(len(MIN))); ax.set_xticklabels([j for _, j in MIN]); ax.set_xlim(-0.5, len(MIN) - 0.5); ax.set_ylim(0.70, 1.005)
ax.set_ylabel("ライブラリとのいちばん高い類似度")
ax.legend(ncol=7, fontsize=10, loc="lower left", bbox_to_anchor=(0.0, 1.0), frameon=False, markerscale=2)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", "applicability.png"), dpi=170)
print("ok")
