"""Figure like Fig. 4a of Teich-McGoldrick, Greathouse and Cygan, J. Phys. Chem. C 116, 15099 (2012): g(r) of the interlayer K
around the bridging O (ob) and the bridging O of a tetrahedral substitution (obts) at 298, 500 and 800 K, 1 atm.  The curves are
the SOVA pair distribution functions of T*/analysis/sova_pdf.csv (soapff-gui, 'K-ob K-obts')."""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
fig, ax = plt.subplots(figsize=(6.2, 4.6))
col = {298: "#1f77b4", 500: "#2ca02c", 800: "#d62728"}
for T in (298, 500, 800):
    rows = list(csv.DictReader(open(f"{HERE}/T{T}/analysis/sova_pdf.csv")))
    r = [float(x["r_A"]) for x in rows]
    for name, ls, lab in (("K-ob", "-", "Oob"), ("K-obts", "--", "Oobts")):
        ax.plot(r, [float(x[name]) for x in rows], ls, color=col[T], lw=1.5, label=f"{T} K, {lab}")
ax.set_xlim(2, 8); ax.set_ylim(bottom=0)
ax.set_xlabel("r (Å)"); ax.set_ylabel("g(r)  K–O")
ax.set_title("Muscovite 2M$_1$, 0.1 MPa", loc="left", fontsize=11)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(f"{HERE}/fig4a.png", dpi=200)
