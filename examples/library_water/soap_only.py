"""Typing by SOAP alone (no shape rule for the water, no marks), on the structures with water of test.py (the 18 held-out last frames, the 10
packmol interfaces, the 6 shapes, the hydrate after its NPT, the 10 hydration states; not the initial hydrate structure, which is a reference).
  A  one step: the SOAP of every atom of the whole structure against the solution library (library_clayff_water_md.npz).
  B  two steps, SOAP only: A finds the water (nearest row ow / hw); those atoms leave the structure; the other atoms are typed with the
     water-free mineral library on the structure without the water (the atoms whose best row is below 0.95 keep their type of A).
The truth is that of test.py.  Counted: wrong types, atoms below the similarity 0.95 (of the type taken), per kind of atom.

    python soap_only.py     (with the Python of soapff-gui)  ->  soap_only.json"""
import collections
import json
import os
import sys

import numpy as np
from ase.io import read

from soapff_cui import core
from soapff_cui._clayff.ff import ff_types

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..")
LIBM = core.load_library(os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz"))
LIBS = core.load_library(os.path.join(HERE, "library_clayff_water_md.npz"))
SYSTEMS = ["water", "NaCl", "KCl", "CsCl", "CaCl2", "BaCl2", "hydrate_2W"] + \
          ["if_" + m for m in ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
                               "pyrophyllite_2", "kaolinite", "lizerdite", "talc", "brucite"]]
alias = lambda t: np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t])
THR = 0.95


def truth_rule(at):
    role, _ = core.roles(at)
    s = at.copy()
    s.set_array("role", np.array(["water" if r == "water" else "ion" if r == "ion" else "bulk" for r in role], dtype="U12"))
    t = alias(ff_types(s))
    el = np.array(at.get_chemical_symbols())
    t = np.where(role == "water", np.where(el == "O", "ow", "hw"), t)
    return np.where(role == "ion", el, t)


def nearest(at, lib):
    """type and similarity of the most similar row (same element) for every atom of `at`, the SOAP of the whole `at`"""
    el = np.array(at.get_chemical_symbols())
    keep = np.nonzero(np.isin(el, lib["species"]))[0]          # the structure seen by the SOAP: the elements of the library
    Fk = core.soap_vectors(at[keep], lib["species"], **lib["soap"]).astype(np.float32)
    F = np.zeros((len(at), Fk.shape[1]), np.float32); F[keep] = Fk
    X = lib["X"].astype(np.float32)
    ty, sim = np.full(len(at), "", dtype=object), np.zeros(len(at))
    for e in sorted(set(el[keep])):
        m = np.nonzero(el == e)[0]; r = np.nonzero(lib["element"] == e)[0]
        for k in range(0, len(m), 4000):
            b = m[k:k + 4000]
            S = F[b] @ X[r].T
            j = S.argmax(1)
            ty[b] = lib["type"][r][j]; sim[b] = S[np.arange(len(b)), j]
    return alias(ty), sim


def kind(t):
    return "water" if t in ("ow", "hw") else "ion" if t in core.FORMAL_ION else "mineral"


ATOMS = []        # per-atom results of every case (soap_only_atoms.npz): element, true type, type and similarity of A, of B
CASES = []


def run(at, truth):
    tA, sA = nearest(at, LIBS)
    water = np.isin(tA, ["ow", "hw"])
    tB, sB = tA.copy(), sA.copy()
    dry = at[~water]
    t1, s1 = nearest(dry, LIBM)
    idx = np.nonzero(~water)[0]
    use = s1 >= THR
    tB[idx[use]] = t1[use]; sB[idx[use]] = s1[use]
    out = {}
    for name, t, s in (("A", tA, sA), ("B", tB, sB)):
        wrong = t != truth
        out[name] = dict(wrong=int(wrong.sum()), low=int((s < THR).sum()),
                         wrong_kinds=dict(collections.Counter(f"{a}->{b}" for a, b in zip(truth[wrong], t[wrong]))),
                         wrong_by_kind=dict(collections.Counter(kind(a) for a in truth[wrong])))
    ATOMS.append((np.array(at.get_chemical_symbols()), truth, tA, sA, tB, sB))
    out["atoms"] = len(at)
    out["water_missed"] = int((np.isin(truth, ["ow", "hw"]) & ~water).sum())          # true water not found as water by SOAP
    out["mineral_as_water"] = int((~np.isin(truth, ["ow", "hw"]) & water).sum())     # mineral or ion found as water
    return out


res = {}
for s in SYSTEMS:
    d = os.path.join(HERE, "refs", s)
    at0, q0, t0 = core.read_lammps_data(os.path.join(d, "system.data"))
    last = core.read_dump_frames(os.path.join(d, "traj.lammpstrj"), at0.get_chemical_symbols())[-1]
    res["held_out/" + s] = run(last, alias(t0)); print(s, res["held_out/" + s], flush=True)

cases = {f"interface_{m}": f"interfaces_001_water/{m}/structure.extxyz" for m in
         ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite", "kaolinite",
          "lizerdite", "talc", "brucite"]}
cases.update({f"shape_{c}": f"interface_shapes/{c}/structure.extxyz" for c in ["sphere", "hemisphere", "cylinder", "cavity", "xygauss", "ions"]})
for name, rel in cases.items():
    at = core.read_structure(os.path.join(EX, rel))
    res["examples/" + name] = run(at, truth_rule(at)); print(name, res["examples/" + name], flush=True)

d = os.path.join(EX, "montmorillonite_2W", "diffusion", "start")
at0, q0, t0 = core.read_lammps_data(os.path.join(d, "system.data"))
fin = read(os.path.join(d, "final.data"), format="lammps-data", style="full")
fin = fin[np.argsort(fin.arrays["id"])] if "id" in fin.arrays else fin
fin.pbc = True
res["examples/hydrate_2W_after_npt"] = run(fin, alias(t0)); print("hydrate_2W_after_npt", res["examples/hydrate_2W_after_npt"], flush=True)

d = os.path.join(EX, "montmorillonite_swelling", "gui_run")
D = open(os.path.join(d, "system.data")).read().splitlines()
k = next(i for i, l in enumerate(D) if l.strip() == "Masses") + 2
names = {}
while D[k].strip():
    p = D[k].split(); names[int(p[0])] = (D[k].split("#")[1].strip(), float(p[1])); k += 1
mass_el = {26.98154: "Al", 28.0855: "Si", 24.305: "Mg", 22.98977: "Na", 15.9994: "O", 1.00794: "H"}
for at in read(os.path.join(d, "swelling_states.lammpstrj"), ":", format="lammps-dump-text"):
    tid = np.asarray(at.arrays["type"], int)
    truth = np.array([names[i][0] for i in tid])
    at.set_chemical_symbols([mass_el[min(mass_el, key=lambda m: abs(m - names[i][1]))] for i in tid])
    at.pbc = True
    nw = int((truth == "ow").sum())
    res[f"examples/interlayer_{nw}"] = run(at, truth); print("interlayer", nw, res[f"examples/interlayer_{nw}"], flush=True)

tot = {k: sum(v["atoms"] if k == "atoms" else v[k]["wrong"] if k in "AB" else 0 for v in res.values()) for k in ("atoms", "A", "B")}
low = {k: sum(v[k]["low"] for v in res.values()) for k in "AB"}
print("TOTAL atoms", tot["atoms"], "wrong A", tot["A"], "B", tot["B"], "below 0.95 A", low["A"], "B", low["B"])
np.savez_compressed(os.path.join(HERE, "soap_only_atoms.npz"),
                    case=np.concatenate([[k] * len(a[0]) for k, a in zip(res, ATOMS)]), element=np.concatenate([a[0] for a in ATOMS]),
                    truth=np.concatenate([a[1] for a in ATOMS]).astype(str), type_A=np.concatenate([a[2] for a in ATOMS]).astype(str),
                    sim_A=np.concatenate([a[3] for a in ATOMS]), type_B=np.concatenate([a[4] for a in ATOMS]).astype(str),
                    sim_B=np.concatenate([a[5] for a in ATOMS]))
json.dump(dict(cases=res, total=dict(atoms=tot["atoms"], wrong_A=tot["A"], wrong_B=tot["B"], low_A=low["A"], low_B=low["B"])),
          open(os.path.join(HERE, "soap_only.json"), "w"), indent=1)
