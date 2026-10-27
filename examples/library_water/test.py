"""Tests of the two-step typing: the mineral with the water-free library (../../src/soapff_gui/data/library_clayff_bulk.npz) on
the structure without the water and the added ions, then the water and the ions with the solution library
(library_clayff_water_md.npz) on the whole structure (soapff_cui.core.make with solution_library, as soapff-gui does).

  (1) held-out frames: the last frame of the MD of every reference system of build_refs.py (not in the library); no marks,
      so the water is separated from the mineral by its shape;
  (2) the structures with water of the other examples: the 10 water-(0 0 1) interfaces (../interfaces_001_water/, packmol
      seed 1), the 6 packmol shapes (../interface_shapes/), the NaCl droplet of Meng et al. (../meng2022_droplet/paper/), the
      2-layer hydrate (../montmorillonite_2W/diffusion/start/, initial structure and after its 100 ps NPT) and the
      montmorillonite of the swelling curve (../montmorillonite_swelling/gui_run/).  The packmol structures are marked as the
      GUI marks them (per-atom array "added": every atom from the first water molecule on, the substrate coming first);
  (2b) the 10 hydration states of the Na montmorillonite of ../montmorillonite_swelling/gui_run/ (4-59 waters added to the interlayer,
      the last structure of the 300 K NPT of each state; the types of the MD input are the true types);
  (3) the dry structures: the 8 sheet-silicate supercells and their 300 K MD frames (../sheet_silicates/md_*/final.data).
Every atom's type (charges.csv) is compared with the reference types: those of the reference system (1); the bond-graph ClayFF
types with the water by its shape and the ions by their element (2, 3).  Reported per structure: atoms, wrong types, fallback,
how the added water and ions were typed (SOAP / rule / below the threshold), the smallest similarity per type.

    python test.py [case ...]     (with the Python of soapff-gui; the Meng droplet, 105,977 atoms, only when named)  ->  test.json"""
import collections
import csv
import json
import os
import sys
import tempfile

import numpy as np
from ase.io import read

from soapff_cui import core
from soapff_cui._clayff.ff import ff_types

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..")
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")
SOL = os.path.join(HERE, "library_clayff_water_md.npz")
SYSTEMS = ["water", "NaCl", "KCl", "CsCl", "CaCl2", "BaCl2", "hydrate_2W"] + \
          ["if_" + m for m in ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
                               "pyrophyllite_2", "kaolinite", "lizerdite", "talc", "brucite"]]
alias = lambda t: np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t])


def truth_rule(at):
    """bond-graph ClayFF types, the water by its shape and the ions by their element (core.roles)"""
    role, _ = core.roles(at)
    s = at.copy()
    s.set_array("role", np.array(["water" if r == "water" else "ion" if r == "ion" else "bulk" for r in role], dtype="U12"))
    t = alias(ff_types(s))
    el = np.array(at.get_chemical_symbols())
    t = np.where(role == "water", np.where(el == "O", "ow", "hw"), t)
    return np.where(role == "ion", el, t)


def mark_packmol(at):
    """the atoms packmol added: from the first water molecule on (the substrate comes first and has no water)"""
    role, _ = core.roles(at)
    w = np.nonzero(role == "water")[0]
    at.set_array("added", np.arange(len(at)) >= (w[0] if len(w) else len(at)))
    return at


ATOMS = []        # per-atom results of every case (test_atoms.npz): group, case, element, true type, assigned type, source, similarity


def run(at, truth, case=("", "")):
    with tempfile.TemporaryDirectory() as d:
        rep = core.make(at, d, LIB, solution_library=SOL, md=False, check=False, log=lambda *a: None)
        rows = list(csv.DictReader(open(os.path.join(d, "charges.csv"))))
    t = alias([r["clayff_type"] for r in rows])
    src = np.array([r["source"] for r in rows])
    sim = np.array([float(r["similarity"]) for r in rows])
    wrong = t != truth
    ATOMS.append((case[0], case[1], np.array(at.get_chemical_symbols()), np.asarray(truth), t, src, sim))
    soap = np.isin(src, ["soap", "soap-low"])
    return dict(atoms=len(at), wrong=int(wrong.sum()),
                wrong_kinds=dict(collections.Counter(f"{a}->{b}" for a, b in zip(truth[wrong], t[wrong]))),
                by_soap=int(soap.sum()), by_rule=int((src == "rule").sum()), fallback=int((src == "fallback").sum()),
                rule_types=dict(collections.Counter(t[src == "rule"])), solution=rep["solution"],
                water_check=rep["water_check"], type_check_differ=rep["type_check"]["differ"],
                consistency=rep["consistency"], charge_sum=round(rep["charge_sum"], 4),
                min_similarity={x: round(float(sim[(t == x) & soap].min()), 4) for x in sorted(set(t[soap]))})


def show(group, name, r):
    print(group, name, {k: r[k] for k in ("atoms", "wrong", "by_soap", "by_rule", "fallback", "rule_types")},
          "solution", r["solution"] and {k: r["solution"][k] for k in ("atoms", "soap", "below_threshold", "mineral_type")},
          flush=True)


only = set(sys.argv[1:])
res = {"held_out": {}, "examples": {}, "dry": {}}
out = os.path.join(HERE, "test.json")
if os.path.isfile(out) and only:
    res = json.load(open(out))
for s in SYSTEMS:
    if only and s not in only:
        continue
    d = os.path.join(HERE, "refs", s)
    at0, q0, t0 = core.read_lammps_data(os.path.join(d, "system.data"))
    last = core.read_dump_frames(os.path.join(d, "traj.lammpstrj"), at0.get_chemical_symbols())[-1]
    res["held_out"][s] = run(last, alias(t0), ("held_out", s))
    show("held out", s, res["held_out"][s])

cases = {f"interface_{m}": f"interfaces_001_water/{m}/structure.extxyz" for m in
         ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite", "kaolinite",
          "lizerdite", "talc", "brucite"]}
cases.update({f"shape_{c}": f"interface_shapes/{c}/structure.extxyz" for c in
              ["sphere", "hemisphere", "cylinder", "cavity", "xygauss", "ions"]})
cases["meng2022_droplet"] = "meng2022_droplet/paper/structure.extxyz"
cases["hydrate_2W_initial"] = "montmorillonite_2W/diffusion/start/system.data"
cases["swelling_gui_run"] = "montmorillonite_swelling/gui_run/system.data"
for name, rel in cases.items():
    if (only and name not in only) or (name == "meng2022_droplet" and name not in only):   # large: only when asked
        continue
    at = core.read_structure(os.path.join(EX, rel))
    if rel.endswith(".extxyz"):
        mark_packmol(at)
    res["examples"][name] = run(at, truth_rule(at), ("examples", name))
    show("example", name, res["examples"][name])

# the 2-layer hydrate after its 100 ps NPT: the types of its start run (same atom order)
if not only or "hydrate_2W_after_npt" in only:
    d = os.path.join(EX, "montmorillonite_2W", "diffusion", "start")
    at0, q0, t0 = core.read_lammps_data(os.path.join(d, "system.data"))
    fin = read(os.path.join(d, "final.data"), format="lammps-data", style="full")
    fin = fin[np.argsort(fin.arrays["id"])] if "id" in fin.arrays else fin
    fin.pbc = True
    res["examples"]["hydrate_2W_after_npt"] = run(fin, alias(t0), ("examples", "hydrate_2W_after_npt"))
    show("example", "hydrate_2W_after_npt", res["examples"]["hydrate_2W_after_npt"])

# the 10 hydration states of the Na montmorillonite (4-59 waters added to the interlayer, the last structure of the NPT of each; the types
# of the MD input are the true types)
if not only or "interlayer" in only:
    d = os.path.join(EX, "montmorillonite_swelling", "gui_run")
    D = open(os.path.join(d, "system.data")).read().splitlines()
    k = next(i for i, l in enumerate(D) if l.strip() == "Masses") + 2
    names = {}
    while D[k].strip():
        p = D[k].split()
        names[int(p[0])] = (D[k].split("#")[1].strip(), float(p[1]))
        k += 1
    mass_el = {26.98154: "Al", 28.0855: "Si", 24.305: "Mg", 22.98977: "Na", 15.9994: "O", 1.00794: "H"}
    for at in read(os.path.join(d, "swelling_states.lammpstrj"), ":", format="lammps-dump-text"):
        tid = np.asarray(at.arrays["type"], int)
        truth = np.array([names[i][0] for i in tid])
        at.set_chemical_symbols([mass_el[min(mass_el, key=lambda m: abs(m - names[i][1]))] for i in tid])
        at.pbc = True
        nw = int((truth == "ow").sum())
        res["examples"][f"interlayer_{nw}"] = run(at, truth, ("examples", f"interlayer_{nw}"))
        show("example", f"interlayer_{nw}", res["examples"][f"interlayer_{nw}"])

for m in ["muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "kaolinite", "lizerdite", "talc", "pyrophyllite_2"]:
    if only and m not in only:
        continue
    sc = read(os.path.join(EX, "sheet_silicates", m + ".extxyz"))
    t = truth_rule(sc)
    res["dry"][m] = run(sc, t, ("dry", m))
    a = read(os.path.join(EX, "sheet_silicates", "md_" + m, "final.data"), format="lammps-data", style="full")
    a = a[np.argsort(a.arrays["id"])] if "id" in a.arrays else a
    a.pbc = True
    res["dry"][m + "_md"] = run(a, t, ("dry", m + "_md"))
    show("dry", m, res["dry"][m]); show("dry", m + "_md", res["dry"][m + "_md"])
json.dump(res, open(out, "w"), indent=1)
if not only:
    np.savez_compressed(os.path.join(HERE, "test_atoms.npz"), group=np.concatenate([[a[0]] * len(a[2]) for a in ATOMS]),
                        case=np.concatenate([[a[1]] * len(a[2]) for a in ATOMS]), element=np.concatenate([a[2] for a in ATOMS]),
                        truth=np.concatenate([a[3] for a in ATOMS]), assigned=np.concatenate([a[4] for a in ATOMS]),
                        source=np.concatenate([a[5] for a in ATOMS]), similarity=np.concatenate([a[6] for a in ATOMS]))
