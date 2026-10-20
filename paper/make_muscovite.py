"""figs/muscovite.png and muscovite.json: muscovite 2M1 (2x2x1, 672 atoms) with the force field made by soapff-gui
(ClayFF, SOAP charges, library of 13 structures). (a) 0 K elastic constants (LAMMPS examples/ELASTIC procedure) against
Brillouin scattering (Vaughan and Guggenheim 1986); (b) uniaxial tension along x, y, z and (c) shear xy, xz, yz at 300 K
(fix deform, 1e9 /s; tension to strain 0.16, shear to 0.2; run from soapff-gui on 8 cores), from the files written by the commands of soapff_cui.properties.
    python make_muscovite.py"""
import json
import sys

import numpy as np

sys.path.insert(0, "/Users/harutokono/clay/soapff-gui/src")
from soapff_cui import properties as P  # noqa: E402

RUN = "/Users/harutokono/clay/soapff-cui/examples/minerals/muscovite_mech/muscovite_2M1_mech_gui"   # GUI run (remote)
RUN0 = "/Users/harutokono/clay/soapff-cui/examples/minerals/muscovite_mech/e0test/e005_np4"         # 0 K, strain 0.005
EXP = dict(C11=181.0, C22=178.4, C33=58.6, C44=16.5, C55=19.5, C66=72.0, C23=21.2, C13=25.6, C12=48.8, C15=-14.2,
           C25=1.1, C35=1.0, C46=-5.2)                      # GPa, Vaughan and Guggenheim (1986)

marks = P.markers(RUN)
e0 = P.elastic0(RUN0, dict(P.markers(RUN0))["elastic0"])
marks = [x for x in marks if x[0] != "elastic0"]
C = e0["C"]
E = np.zeros((6, 6))
for k, x in EXP.items():
    i, j = int(k[1]) - 1, int(k[2]) - 1
    E[i, j] = E[j, i] = x
out = dict(C_clayff_0K=np.round(C, 2).tolist(), C_exp=E.tolist(), moduli_clayff_0K=P._moduli(C), moduli_exp=P._moduli(E),
           tests={})
curves = {}
marks = [(n, dict(m, _dir=RUN)) for n, m in marks]          # tension to strain 0.16, shear to 0.2 (failure)
for name, m in marks:
    if name not in ("tensile", "shear"):
        continue
    try:
        r = P.mechanical(m["_dir"], m, 0.5)
    except (OSError, IndexError, ValueError):
        continue
    key = f"tension_{m['dir']}" if name == "tensile" else f"shear_{m['plane']}"
    curves[key] = (r["strain"], r["stress"])
    out["tests"][key] = {k: v for k, v in r.items() if k not in ("strain", "stress")}
    ts = np.convolve(r["stress"], np.ones(10) / 10, mode="same")       # 0.5 ps moving average for the strength
    i = int(np.argmax(ts[10:-10])) + 10
    out["tests"][key].update(peak_avg=float(ts[i]), peak_avg_strain=float(r["strain"][i]))
json.dump(out, open("muscovite.json", "w"), indent=1)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
labels = ["C11", "C22", "C33", "C44", "C55", "C66", "C12", "C13", "C23"]
x = np.arange(len(labels))
ax[0].bar(x - 0.2, [C[int(k[1]) - 1, int(k[2]) - 1] for k in labels], 0.4, label="ClayFF (this work, 0 K)", color="#28788d")
ax[0].bar(x + 0.2, [EXP[k] for k in labels], 0.4, label="experiment (Brillouin)", color="#bbbbbb")
ax[0].set_xticks(x); ax[0].set_xticklabels(labels); ax[0].set_ylabel("elastic constant (GPa)")
ax[0].legend(fontsize=9); ax[0].set_title("(a) elastic constants", loc="left", fontsize=11)
for a, kind, title in ((ax[1], "tension", "(b) uniaxial tension, 300 K"), (ax[2], "shear", "(c) shear, 300 K")):
    for key, (s, t) in sorted(curves.items()):
        if key.startswith(kind):
            a.plot(s, t, lw=1.3, label=key.split("_")[1])
    a.set_xlabel("engineering strain" if kind == "tension" else "engineering shear strain")
    a.set_ylabel("stress (GPa)" if kind == "tension" else "shear stress (GPa)")
    a.legend(fontsize=9); a.set_title(title, loc="left", fontsize=11); a.axhline(0, color="#999", lw=0.6)
fig.tight_layout()
fig.savefig("figs/muscovite.png", dpi=200)
print(json.dumps(out["tests"], indent=1))
