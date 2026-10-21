"""Figure 4a of Teich-McGoldrick, Greathouse and Cygan, J. Phys. Chem. C 116, 15099 (2012) again with soapff-gui: g(r) of the
interlayer K around all O of muscovite 2M1 at 298, 500 and 800 K, 0.1 MPa (SOVA pair distribution function 'all:K-O' of
T*/analysis/sova_pdf.csv; 'K-ob' and 'K-obts' for the split by ClayFF type).
  fig4a.png          our curves in the style of the paper's figure
  fig4a_overlay.png  our curves (thin lines) on the paper's figure paper_fig4a.png (axes read from its frame: x 0 -> pixel 91.5,
                     10 -> 585; g 0 -> pixel 375, 8 -> 30)
  fig4a_split.png    K-ob and K-obts separately
    python make_fig4a.py"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
COL = {298: "black", 500: "red", 800: "green"}


def pdf(T, name):
    rows = list(csv.DictReader(open(f"{HERE}/T{T}/analysis/sova_pdf.csv")))
    return np.array([float(x["r_A"]) for x in rows]), np.array([float(x[name]) for x in rows])


# our figure
fig, ax = plt.subplots(figsize=(6.0, 4.4))
for T in (298, 500, 800):
    r, g = pdf(T, "all:K-O")
    ax.plot(r, g, color=COL[T], lw=1.6, label=f"T = {T} K")
ax.set_xlim(0, 10); ax.set_ylim(0, 8)
ax.set_xlabel("r (Å)"); ax.set_ylabel("g(r)  K – O")
ax.text(0.2, 7.2, "P = 1×10$^{-4}$ GPa", fontsize=12)
ax.legend(loc="upper right", bbox_to_anchor=(0.97, 0.93), frameon=False)
ax.set_title("a)  soapff-gui + SOVA (672 atoms, 100 ps)", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(f"{HERE}/fig4a.png", dpi=200); plt.close(fig)

# overlay on the paper's figure
im = plt.imread(f"{HERE}/paper_fig4a.png")
H, W = im.shape[:2]
x0, x1, y0, y1 = 91.5, 585.0, 375.0, 30.0
sx, sy = 10.0 / (x1 - x0), 8.0 / (y0 - y1)
fig, ax = plt.subplots(figsize=(7.0, 5.0))
ax.imshow(im, extent=(-x0 * sx, (W - x0) * sx, -(H - y0) * sy, y0 * sy), aspect="auto")
for T in (298, 500, 800):
    r, g = pdf(T, "all:K-O")
    ax.plot(r, g, color={298: "#1f77b4", 500: "#ff7f0e", 800: "#9467bd"}[T], lw=1.2, ls="--", label=f"soapff-gui {T} K")
ax.set_xlim(-x0 * sx, (W - x0) * sx); ax.set_ylim(-(H - y0) * sy, y0 * sy); ax.axis("off")
ax.legend(loc="upper right", bbox_to_anchor=(0.98, 0.78), fontsize=8, frameon=True)
fig.tight_layout(pad=0); fig.savefig(f"{HERE}/fig4a_overlay.png", dpi=200); plt.close(fig)

# split by type
fig, ax = plt.subplots(figsize=(6.2, 4.6))
for T in (298, 500, 800):
    for name, ls, lab in (("K-ob", "-", "Oob"), ("K-obts", "--", "Oobts")):
        r, g = pdf(T, name)
        ax.plot(r, g, ls, color=COL[T], lw=1.4, label=f"{T} K, {lab}")
ax.set_xlim(2, 8); ax.set_ylim(bottom=0); ax.set_xlabel("r (Å)"); ax.set_ylabel("g(r)  K–O"); ax.legend(fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(f"{HERE}/fig4a_split.png", dpi=200)

# peak heights and positions against the paper's (read from the figure)
PAPER = {298: (5.3, 5.1, 3.45), 500: (4.5, 4.15, 2.95), 800: (3.7, 3.4, 2.5)}
for T in (298, 500, 800):
    r, g = pdf(T, "all:K-O")
    out = []
    for (a, b), p in zip(((2.5, 3.8), (4.2, 5.4), (6.5, 7.6)), PAPER[T]):
        m = (r > a) & (r < b); i = np.argmax(g[m])
        out.append(f"{r[m][i]:.2f} Å  g {g[m][i]:.2f} (paper {p})")
    print(T, " | ".join(out))
