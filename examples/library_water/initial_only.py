"""Initial structures only (no MD frame): the 10 packmol interfaces and the 6 water shapes, as built (unrelaxed).
  1. mineral atoms and ions of the 10 interfaces typed by their SOAP against (i) the water-free mineral library on the whole structure,
     (ii) the same on the structure without the water, (iii) the solution library on the whole structure  -> initial_only.json (inwater)
  2. all atoms of the 16 structures typed by the whole-structure SOAP against the solution library: wrong types, atoms below 0.95, and the
     minimum similarity per kind of atom (from soap_only_atoms.npz, scheme A)  -> initial_only.json (water)
    python initial_only.py     (with the Python of soapff-gui)"""
import json
import os

import numpy as np

from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..")
THR = 0.95
LIBM = core.load_library(os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz"))
LIBS = core.load_library(os.path.join(HERE, "library_clayff_water_md.npz"))
alias = lambda t: np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t])


def truth_rule(at):
    from soapff_cui._clayff.ff import ff_types
    role, _ = core.roles(at)
    s = at.copy()
    s.set_array("role", np.array(["water" if r == "water" else "ion" if r == "ion" else "bulk" for r in role], dtype="U12"))
    t = alias(ff_types(s))
    el = np.array(at.get_chemical_symbols())
    t = np.where(role == "water", np.where(el == "O", "ow", "hw"), t)
    return np.where(role == "ion", el, t)


def nearest(at, lib):
    el = np.array(at.get_chemical_symbols())
    keep = np.nonzero(np.isin(el, lib["species"]))[0]
    Fk = core.soap_vectors(at[keep], lib["species"], **lib["soap"]).astype(np.float32)
    F = np.zeros((len(at), Fk.shape[1]), np.float32); F[keep] = Fk
    X = lib["X"].astype(np.float32)
    ty, sim = np.full(len(at), "", dtype=object), np.zeros(len(at))
    for e in sorted(set(el[keep])):
        m = np.nonzero(el == e)[0]; r = np.nonzero(lib["element"] == e)[0]
        S = F[m] @ X[r].T
        j = S.argmax(1)
        ty[m] = lib["type"][r][j]; sim[m] = S[np.arange(len(m)), j]
    return alias(ty), sim


tot = {k: dict(n=0, low=0, wrong=0, min=1.0) for k in ("mineral_whole", "mineral_dry", "solution_whole")}
per = {}
for m in ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite", "kaolinite", "lizerdite", "talc", "brucite"]:
    at = core.read_structure(os.path.join(EX, "interfaces_001_water", m, "structure.extxyz"))
    tr = truth_rule(at)
    wat = np.isin(tr, ["ow", "hw"])
    keep = ~wat
    out = {}
    t1, s1 = nearest(at, LIBM); out["mineral_whole"] = (t1[keep], s1[keep])
    t2, s2 = nearest(at[keep], LIBM); out["mineral_dry"] = (t2, s2)
    t3, s3 = nearest(at, LIBS); out["solution_whole"] = (t3[keep], s3[keep])
    for k, (t, s) in out.items():
        w = int((t != tr[keep]).sum()); lo = int((s < THR).sum())
        tot[k]["n"] += int(keep.sum()); tot[k]["wrong"] += w; tot[k]["low"] += lo; tot[k]["min"] = min(tot[k]["min"], float(s.min()))
    per[m] = {k: dict(wrong=int((t != tr[keep]).sum()), low=int((s < THR).sum()), min=float(s.min())) for k, (t, s) in out.items()}
    print(m, int(keep.sum()), per[m], flush=True)
print(tot)

d = np.load(os.path.join(HERE, "soap_only_atoms.npz"), allow_pickle=True)
case = d["case"]; sel = np.array([str(c).startswith(("examples/interface_", "examples/shape_")) for c in case])
sel &= ~np.array([str(c).startswith("examples/interface_pyro") and False for c in case])
groups = {"interface": [], "shape": []}
res = {}
for g, pre in (("interface", "examples/interface_"), ("shape", "examples/shape_")):
    m = np.array([str(c).startswith(pre) for c in case])
    tr, ty, si, el = d["truth"][m], d["type_A"][m], d["sim_A"][m], d["element"][m]
    water = np.isin(tr, ["ow", "hw"]); ion = np.isin(tr, ["Na", "K", "Cs", "Ca", "Ba", "Cl"])
    res[g] = dict(structures=len(set(case[m])), atoms=int(m.sum()), ow=int((tr == "ow").sum()), ions=int(ion.sum()),
                  wrong=int((ty != tr).sum()), low=int((si < THR).sum()),
                  min_water_ion=float(si[water | ion].min()), min_mineral=float(si[~water & ~ion].min()))
    print(g, res[g])
json.dump(dict(inwater=dict(total=tot, per=per), water=res), open(os.path.join(HERE, "initial_only.json"), "w"), indent=1)
