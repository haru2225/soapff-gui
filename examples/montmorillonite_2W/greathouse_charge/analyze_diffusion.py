#!/usr/bin/env python3
"""D_xy of water O and Na+ of the Greathouse-charge model (T300, T330, T366 made and run with soapff-gui): each folder has
4 segments of 200 ps (msd_s<k>_Ow.dat, msd_s<k>_Na.dat; compute msd defined at the start of each segment); 50-200 ps linear
fit of the in-plane MSD, mean and standard error of the 4 segments, Arrhenius activation energy; Greathouse et al. (2016)
J. Phys. Chem. C 120, 1640 (Na-montmorillonite, 0.75 e / O20(OH)4, 2W) for comparison.
    python analyze_diffusion.py [FIGURE.png] -> diffusion_summary.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
DT_FS = 0.5; FIT = (50.0, 200.0); R = 8.314e-3
G = {"Ow": {300: 13.1, 330: 24.3, 366: 39.4}, "Na": {300: 5.1, 330: 10.4, 366: 17.0}}   # 1e-10 m2/s
G_EA = {"Ow": 15.2, "Na": 16.5}
res = {}
for g in ("Ow", "Na"):
    res[g] = {}
    for T in (300, 330, 366):
        Ds = []
        for s in (1, 2, 3, 4):
            d = np.loadtxt(f"{HERE}/T{T}/msd_s{s}_{g}.dat", comments="#")
            t = (d[:, 0] - d[0, 0]) * DT_FS / 1000.0 + 0.1          # first output 200 steps (0.1 ps) after the origin
            m = (t >= FIT[0]) & (t <= FIT[1])
            Ds.append(np.polyfit(t[m], d[m, 1] + d[m, 2], 1)[0] / 4 * 100)
        Ds = np.array(Ds)
        res[g][T] = dict(segments=[round(float(x), 3) for x in Ds], mean=float(Ds.mean()), sem=float(Ds.std(ddof=1) / 2))
    Tk = np.array([300, 330, 366.0]); D = np.array([res[g][T]["mean"] for T in (300, 330, 366)])
    p, cov = np.polyfit(1 / Tk, np.log(D), 1, cov=True)
    res[g]["Ea"] = float(-p[0] * R); res[g]["Ea_err"] = float(np.sqrt(cov[0, 0]) * R)
    res[g]["Ea_greathouse_fit"] = float(-np.polyfit(1 / Tk, np.log([G[g][T] for T in (300, 330, 366)]), 1)[0] * R)
# layer spacing L_z / 2 over the production (log.lammps thermo, last column lz)
res["d001"] = {}
for T in (300, 330, 366):
    rows = [l.split() for l in open(f"{HERE}/T{T}/log.lammps") if l.split()[:1] and l.split()[0].isdigit() and len(l.split()) == 11]
    lz = np.array([float(r[10]) for r in rows[-400:]])
    res["d001"][T] = float(lz.mean() / 2)
json.dump(res, open(f"{HERE}/diffusion_summary.json", "w"), indent=1)
for g in ("Ow", "Na"):
    for T in (300, 330, 366):
        r = res[g][T]
        print(g, T, np.round(r["segments"], 2), "mean %.2f +- %.2f (Greathouse %.1f, ratio %.2f)" % (r["mean"], r["sem"], G[g][T], r["mean"] / G[g][T]))
    print(g, "Ea %.1f +- %.1f kJ/mol (Greathouse %.1f, fit %.1f)" % (res[g]["Ea"], res[g]["Ea_err"], G_EA[g], res[g]["Ea_greathouse_fit"]))
print("d001", {k: round(v, 2) for k, v in res["d001"].items()})
