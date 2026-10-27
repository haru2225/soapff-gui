"""The solution library of the edge library "Montmorillonite edge (Underwood & Bourg 2020)" (library_clayff_edge_water.npz): the two particle files of the
edge library (../library/refs_edge/), and the structures that packmol builds with water and cations, not relaxed, no MD frame: liquid water, the five
chloride solutions (NaCl, KCl, CsCl, CaCl2, BaCl2; refs/<name>/<name>.data of build_refs.py, ClayFF water and ions) and the particle with a
water film and 4 Na+ + 4 Cl- (../montmorillonite_edge/out_particle_water/system.data: the particle and its Na+ the edge model, the water and Cl- ClayFF).
The water and the ions are ClayFF in every system (the model of the particle has neither Cl- nor water).  Rows of the systems with water
are pruned at 0.99, the second particle file at 0.999 (as the edge library).

    python build_library_clayff_edge_water.py     (with the Python of soapff-gui)  ->  library_clayff_edge_water.npz, build_library_clayff_edge_water.json"""
import json
import os

from soapff_cui import core
from build_library import HERE, reference

LIBDIR = os.path.join(HERE, "..", "library", "refs_edge")
SYSTEMS = ["water", "NaCl", "KCl", "CsCl", "CaCl2", "BaCl2"]
refs = [os.path.join(LIBDIR, "clay_montmorillonite.data"), os.path.join(LIBDIR, "clay_montmorillonite_2.data")] + \
       [reference(s) for s in SYSTEMS] + [os.path.join(HERE, "..", "montmorillonite_edge", "out_particle_water", "system.data")]
prune = [None, 0.999] + [0.99] * (len(SYSTEMS) + 1)
out = os.path.join(HERE, "library_clayff_edge_water.npz")
core.build_library(refs, [], out=out, water_in_soap=True, prune=prune, protect=1)
lib = core.load_library(out)
r = core.check_library(out)
info = dict(rows=int(len(lib["q"])), species=list(lib["species"]), rows_per_type={t: int((lib["type"] == t).sum()) for t in sorted(set(lib["type"]))},
            one_to_one=dict(rows=r["total"], same_charge=r["total"] - r["wrong"], confusions=r["confusions"]))
json.dump(info, open(os.path.join(HERE, "build_library_clayff_edge_water.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
