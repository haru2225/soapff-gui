"""Summary of the NVT and NPT runs of the 11 bulk library minerals (each <name>/run, made with soapff-gui):
force field (fallback, type check), temperature, potential energy and volume of the NVT (50 ps) and the NPT (50 ps).
    python summarize.py   -> summary.csv and the table printed (as in README.md)"""
import csv
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def blocks(log):
    """thermo blocks of log.lammps: list of arrays (step temp press pe ke etotal vol density lx ly lz)."""
    out, cur = [], None
    for line in open(log):
        w = line.split()
        if w[:2] == ["Step", "Temp"]:
            cur = []
            out.append(cur)
        elif cur is not None and line.startswith("Loop time"):
            cur = None
        elif cur is not None and w and w[0].isdigit():
            cur.append([float(x) for x in w])
    return [np.array(b) for b in out]


rows = []
for d in sorted(glob.glob(f"{HERE}/*/run")):
    name = os.path.basename(os.path.dirname(d))
    rep = json.load(open(f"{d}/report.json"))
    b = blocks(f"{d}/log.lammps")
    nvt, npt = b[-2], b[-1]                              # minimisations, relaxation NVT, NVT 50 ps, NPT 50 ps
    n = rep["n_atoms"] if "n_atoms" in rep else sum(1 for _ in open(f"{d}/charges.csv")) - 1
    h = lambda a: a[len(a) // 2:]                         # second half
    r = dict(mineral=name, atoms=n, fallback=rep["fallback"], type_differ=rep["type_check"]["differ"],
             T_nvt=h(nvt)[:, 1].mean(), T_npt=h(npt)[:, 1].mean(),
             pe_drift_npt=(h(npt)[-len(h(npt)) // 2:, 3].mean() - h(npt)[:len(h(npt)) // 2, 3].mean()) / n,
             P_npt=h(npt)[:, 2].mean(), dV_npt_pct=(h(npt)[:, 6].mean() / nvt[0, 6] - 1) * 100,
             density_npt=h(npt)[:, 7].mean())
    rows.append(r)
with open(f"{HERE}/summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print("| mineral | atoms | fallback | type mismatches | T NVT (K) | T NPT (K) | P NPT (atm) | volume change (%) | density (g/cm³) | PE change (kcal/mol/atom) |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r['mineral']} | {r['atoms']} | {r['fallback']} | {r['type_differ']} | {r['T_nvt']:.1f} | {r['T_npt']:.1f} | "
          f"{r['P_npt']:.0f} | {r['dV_npt_pct']:+.2f} | {r['density_npt']:.3f} | {r['pe_drift_npt']:+.4f} |")
