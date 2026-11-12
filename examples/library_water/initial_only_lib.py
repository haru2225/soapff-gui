"""Whole-structure SOAP of the 16 initial structures (10 packmol interfaces, 6 water shapes) against the library of the initial structures only
(library_clayff_water.npz): wrong types, atoms below 0.95, minimum similarity per kind of atom -> initial_only_lib.json
    python initial_only_lib.py     (with the Python of soapff-gui)"""
import json
import os
import numpy as np
from soapff_cui import core
from initial_only import nearest, truth_rule, EX, HERE, THR
LIB = core.load_library(os.path.join(HERE, "library_clayff_water.npz"))
MIN = ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite", "kaolinite", "lizerdite", "talc", "brucite"]
cases = {f"interface_{m}": f"interfaces_001_water/{m}/structure.extxyz" for m in MIN}
cases.update({f"shape_{c}": f"interface_shapes/{c}/structure.extxyz" for c in ["sphere", "hemisphere", "cylinder", "cavity", "xygauss", "ions"]})
rows = {}
for name, rel in cases.items():
    at = core.read_structure(os.path.join(EX, rel))
    tr = truth_rule(at)
    t, s = nearest(at, LIB)
    water = np.isin(tr, ["ow", "hw"]); ion = np.isin(tr, ["Na", "K", "Cs", "Ca", "Ba", "Cl"]); mine = ~water & ~ion
    rows[name] = dict(atoms=len(at), water=int((tr == "ow").sum()), ions=int(ion.sum()), wrong=int((t != tr).sum()), low=int((s < THR).sum()),
                      wrong_types=[f"{a}->{b}" for a, b in zip(tr[t != tr], t[t != tr])],
                      min_water=float(s[water].min()), min_ion=float(s[ion].min()) if ion.any() else None, min_mineral=float(s[mine].min()) if mine.any() else None)
    print(name, rows[name], flush=True)
g = lambda pre: [v for k, v in rows.items() if k.startswith(pre)]
tot = {}
for pre in ("interface", "shape"):
    r = g(pre)
    tot[pre] = dict(structures=len(r), atoms=sum(v["atoms"] for v in r), water=sum(v["water"] for v in r), ions=sum(v["ions"] for v in r),
                    wrong=sum(v["wrong"] for v in r), low=sum(v["low"] for v in r),
                    min_water_ion=min(min(v["min_water"], v["min_ion"] or 1) for v in r), min_mineral=min(v["min_mineral"] for v in r))
print(tot)
json.dump(dict(cases=rows, total=tot), open(os.path.join(HERE, "initial_only_lib.json"), "w"), indent=1)
