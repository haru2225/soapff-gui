"""The SOAP libraries of the paper, rebuilt from the reference structures in this folder, and their one-to-one checks.

  (a) the bundled default library (src/soapff_gui/data/library_clayff_bulk.npz): the 11 bulk layered minerals (refs_bulk/) and
      their basal (0 0 1) slabs (refs_basal/), water left out of the SOAP, the slab rows pruned at similarity 0.999;
  (b) the basal slabs cut again from the bulk references with soapff_cui.surfaces (pymatgen), compared with refs_basal/;
  (c) the optional library "Montmorillonite edge (Underwood & Bourg 2020)" (library_clayff_edge.npz): the particle with
      edges (refs_edge/) and its force field taken whole from its data and settings files;
  (d) the initial 12-structure library of the paper (Table "one-to-one"): the 2-layer hydrate (refs_hydrate/, water and
      Na+ in the SOAP), the dry Na montmorillonite and the 10 other bulk crystals, checked at l_max = 3, 4, 5, and the hydrate
      alone.
Each check: the nearest other row of the same element carries the same charge (core.check_library).

    python run.py            (with the Python of soapff-gui; about 10 minutes)  ->  results.json, *.npz"""
import json
import os
import sys

import numpy as np

from soapff_cui import core, surfaces

HERE = os.path.dirname(os.path.abspath(__file__))
BULK = ["brucite_bulk", "kaolinite_bulk", "lizerdite_bulk", "muscovite_1M_bulk", "muscovite_2M1_bulk", "nanpingite_bulk",
        "phlogopite_1M_bulk", "pyrophyllite_2_bulk", "pyrophyllite_bulk", "talc_bulk", "montmorillonite_dry_Na"]
BASAL = ["brucite", "kaolinite", "lizerdite", "muscovite_1M", "muscovite_2M1", "nanpingite", "phlogopite_1M", "pyrophyllite_2",
         "pyrophyllite", "talc", "montmorillonite_dry_Na"]
BUNDLED = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data")
quiet = lambda *a: None
res = {}


def same_library(a, b):
    """True when the two npz files hold the same rows (SOAP vectors to 1e-6, charges, types, elements)"""
    x, y = np.load(a), np.load(b)
    return bool(x["X"].shape == y["X"].shape and np.allclose(x["X"], y["X"], atol=1e-6) and np.allclose(x["q"], y["q"])
                and (x["type"] == y["type"]).all() and (x["element"] == y["element"]).all())


def check(path):
    r = core.check_library(path, log=quiet)
    return dict(rows=r["total"], same_charge=r["total"] - r["wrong"], confusions=r["confusions"])


# (a) the default library
refs = [os.path.join(HERE, "refs_bulk", b + ".data") for b in BULK]
slabs = [os.path.join(HERE, "refs_basal", m + "_001.data") for m in BASAL]
out = os.path.join(HERE, "library_clayff_bulk.npz")
core.build_library(refs + slabs, out=out, water_in_soap=False, prune=0.999, protect=len(refs), log=quiet)
res["default"] = dict(check(out), structures=len(refs) + len(slabs),
                      same_as_bundled=same_library(out, os.path.join(BUNDLED, "library_clayff_bulk.npz")))
print("default library", res["default"], flush=True)

# (b) the basal slabs cut again (pymatgen; the bulk reference without water, 16 A of mineral, 16 A of vacuum)
res["basal_recut"] = {}
for b, m in zip(BULK, BASAL):
    at, q, types = core.read_lammps_data(os.path.join(HERE, "refs_bulk", b + ".data"))
    keep = ~np.isin(types, ["ow", "hw"])
    slab, typ, info = surfaces.make_slab(at[keep], [t for t, k in zip(types, keep) if k], (0, 0, 1), min_thickness=16.0,
                                         vacuum=16.0)
    a0, q0, t0 = core.read_lammps_data(os.path.join(HERE, "refs_basal", m + "_001.data"))
    res["basal_recut"][m] = dict(atoms=len(slab), same_atoms=len(slab) == len(a0), same_types=sorted(typ) == sorted(t0),
                                 formal_charge=info["formal_charge"], ions_face=info["ions_face"])
print("basal slabs cut again", all(v["same_atoms"] and v["same_types"] for v in res["basal_recut"].values()), flush=True)

# (c) the optional edge library: the particle files only, with their force field
edge = [os.path.join(HERE, "refs_edge", e) for e in ("clay_montmorillonite.data", "clay_montmorillonite_2.data")]
ff = core.read_forcefield(edge[0], os.path.join(HERE, "refs_edge", "clay_montmorillonite.settings"), "Montmorillonite edge",
                          "T. R. Underwood, I. C. Bourg, J. Phys. Chem. C 124, 3702 (2020)")
out = os.path.join(HERE, "library_clayff_edge.npz")
core.build_library(edge, out=out, water_in_soap=False, prune=0.999, protect=1, forcefield=ff, log=quiet)
res["edge"] = dict(check(out), structures=2, same_as_bundled=same_library(out, os.path.join(BUNDLED, "library_clayff_edge.npz")),
                   forcefield=ff)
print("edge library", {k: v for k, v in res["edge"].items() if k != "forcefield"}, flush=True)

# (d) the initial 12-structure library (water in the SOAP), l_max = 3, 4, 5; and the hydrate alone
hydrate = os.path.join(HERE, "refs_hydrate", "montmorillonite_2W.data")
initial = [hydrate, os.path.join(HERE, "refs_bulk", "montmorillonite_dry_Na.data")] + \
          [os.path.join(HERE, "refs_bulk", b + ".data") for b in BULK if b != "montmorillonite_dry_Na"]
res["initial"] = {}
for lmax in (3, 4, 5):
    for tag, files in (("hydrate", [hydrate]), ("library_12", initial)):
        out = os.path.join(HERE, f"_initial_{tag}_l{lmax}.npz")
        core.build_library(files, out=out, water_in_soap=True, soap=dict(l_max=lmax), log=quiet)
        res["initial"][f"{tag}_lmax{lmax}"] = check(out)
        os.remove(out)
        print(tag, "l_max", lmax, res["initial"][f"{tag}_lmax{lmax}"], flush=True)
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
