#!/usr/bin/env python3
"""D_xy of water O and Na+ from the LAMMPS compute msd of the 4 segments (200 ps) at each temperature, 50-200 ps linear fit,
block mean and standard error, Arrhenius activation energy; Greathouse et al. (2016) J. Phys. Chem. C 120, 1640, Table 3
(Na-montmorillonite, -0.750 e/u.c., 2W; the layer charge and the water content of this system, 0.75 e and 8.0 H2O per O20(OH)4, are the same) for comparison.  python analyze_diffusion.py [FIGURE.png]  -> diffusion_summary.json, figure"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
DT_FS = 0.5; FIT = (50.0, 200.0); R = 8.314e-3
G = {"Ow": {300: 13.1, 330: 24.3, 366: 39.4}, "Na": {300: 5.1, 330: 10.4, 366: 17.0}}   # 1e-10 m2/s
G_EA = {"Ow": 15.2, "Na": 16.5}                                                         # kJ/mol, reported
res = {}
temps = [T for T in (300, 330, 366) if all(os.path.isfile(f"{HERE}/T{T}/seg{s}/msd_Ow.dat") and os.path.isfile(f"{HERE}/T{T}/seg{s}/final.data") for s in (1, 2, 3, 4))]
for g in ("Ow", "Na"):
    res[g] = {}
    for T in temps:
        Ds = []
        for s in (1, 2, 3, 4):
            d = np.loadtxt(f"{HERE}/T{T}/seg{s}/msd_{g}.dat", comments="#"); t = d[:, 0] * DT_FS / 1000.0
            m = (t >= FIT[0]) & (t <= FIT[1])
            Ds.append(np.polyfit(t[m], d[m, 1] + d[m, 2], 1)[0] / 4 * 100)         # A^2/ps -> 1e-10 m2/s
        Ds = np.array(Ds)
        res[g][T] = dict(segments=[round(float(x), 3) for x in Ds], mean=float(Ds.mean()), sem=float(Ds.std(ddof=1) / 2))
    if len(temps) == 3:
        Tk = np.array(temps, float); D = np.array([res[g][T]["mean"] for T in temps])
        p, cov = np.polyfit(1 / Tk, np.log(D), 1, cov=True)
        res[g]["Ea"] = float(-p[0] * R); res[g]["Ea_err"] = float(np.sqrt(cov[0, 0]) * R)
        res[g]["Ea_greathouse_fit"] = float(-np.polyfit(1 / Tk, np.log([G[g][T] for T in temps]), 1)[0] * R)
        res[g]["Ea_greathouse_reported"] = G_EA[g]
# layer spacing d001 = L_z (one layer per cell) over the last 50 ps of each equilibration/segment
lz = {}
for T in temps:
    v = []
    for s in (1, 2, 3, 4):
        rows = [l.split() for l in open(f"{HERE}/T{T}/seg{s}/log.lammps") if l.split()[:1] and l.split()[0].isdigit() and len(l.split()) == 11]
        v.append(np.mean([float(r[-1]) for r in rows[-25:]]))
    lz[T] = float(np.mean(v))
res["d001"] = lz
json.dump(res, open(f"{HERE}/diffusion_summary.json", "w"), indent=1)
for g in ("Ow", "Na"):
    for T in temps:
        r = res[g][T]; print(g, T, np.round(r["segments"], 2), "mean %.2f +- %.2f (Greathouse %.1f, ratio %.2f)" % (r["mean"], r["sem"], G[g][T], r["mean"] / G[g][T]))
    if "Ea" in res[g]:
        print(g, "Ea %.1f +- %.1f kJ/mol (Greathouse reported %.1f, fit %.1f)" % (res[g]["Ea"], res[g]["Ea_err"], G_EA[g], res[g]["Ea_greathouse_fit"]))
print("d001", {k: round(v, 2) for k, v in lz.items()})
if len(sys.argv) > 1:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "Hiragino Sans"
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); col = {300: "C0", 330: "C1", 366: "C3"}
    for k, (g, lab) in enumerate((("Ow", "water O"), ("Na", "Na$^+$"))):
        for T in temps:
            ms = []
            for s in (1, 2, 3, 4):
                d = np.loadtxt(f"{HERE}/T{T}/seg{s}/msd_{g}.dat", comments="#"); t = d[:, 0] * DT_FS / 1000; ms.append(d[:, 1] + d[:, 2])
            ms = np.mean(ms, 0); ax[k].plot(t, ms, color=col[T], label=f"{T} K")
            m = (t >= FIT[0]) & (t <= FIT[1]); p = np.polyfit(t[m], ms[m], 1); ax[k].plot(t[m], np.polyval(p, t[m]), "k:", lw=0.8)
        ax[k].set_xlabel("time (ps)"); ax[k].set_ylabel("in-plane MSD (Å$^2$)"); ax[k].set_title(f"({'ab'[k]}) {lab}"); ax[k].legend(fontsize=8)
    Tk = np.array(temps)
    for g, mk, lab, c in (("Ow", "o", "water", "C0"), ("Na", "s", "Na$^+$", "C2")):
        ax[2].errorbar(1000 / Tk, [res[g][T]["mean"] for T in temps], yerr=[res[g][T]["sem"] for T in temps], fmt=mk + "-", color=c, label=f"{lab} (this work)", capsize=3)
        ax[2].plot(1000 / Tk, [G[g][T] for T in temps], mk + "--", mfc="none", color="gray", label=f"{lab} (Greathouse et al.)")
    ax[2].set_yscale("log"); ax[2].set_xlabel("1000 / T (K$^{-1}$)"); ax[2].set_ylabel("$D_{xy}$ (10$^{-10}$ m$^2$/s)"); ax[2].set_title("(c) in-plane diffusion coefficient"); ax[2].legend(fontsize=7)
    plt.tight_layout(); plt.savefig(sys.argv[1], dpi=200)
