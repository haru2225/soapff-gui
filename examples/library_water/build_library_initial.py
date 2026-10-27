"""The library with water and aqueous ions from the initial structures only (library_clayff_water.npz): the 11 bulk layered minerals and their
basal (0 0 1) slabs of the water-free library, plus the structures that packmol builds (build_refs.py; not relaxed, no MD frame): liquid water, the
five chloride solutions (NaCl, KCl, CsCl, CaCl2, BaCl2) and the 11 slabs with water packed above them.  The 2-layer hydrate is not included (its structure
comes from an NPT run).  Rows of the systems with water are pruned at 0.99, the slab rows at 0.999 (as in build_library.py).

    python build_library_initial.py     (with the Python of soapff-gui)  ->  library_clayff_water.npz, build_library_initial.json"""
import json
import os

from soapff_cui import core
from build_library import BASAL, BULK, LIBDIR, HERE, reference

SYSTEMS = ["water", "NaCl", "KCl", "CsCl", "CaCl2", "BaCl2"] + \
          ["if_" + m for m in ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
                               "pyrophyllite_2", "kaolinite", "lizerdite", "talc", "brucite"]]
refs = [os.path.join(LIBDIR, "refs_bulk", b + ".data") for b in BULK] + \
       [os.path.join(LIBDIR, "refs_basal", m + "_001.data") for m in BASAL] + [reference(s) for s in SYSTEMS]
prune = [None] * len(BULK) + [0.999] * len(BASAL) + [0.99] * len(SYSTEMS)
out = os.path.join(HERE, "library_clayff_water.npz")
core.build_library(refs, [], out=out, water_in_soap=True, prune=prune, protect=len(BULK))
lib = core.load_library(out)
r = core.check_library(out)
info = dict(rows=int(len(lib["q"])), species=list(lib["species"]), rows_per_type={t: int((lib["type"] == t).sum()) for t in sorted(set(lib["type"]))},
            one_to_one=dict(rows=r["total"], same_charge=r["total"] - r["wrong"], confusions=r["confusions"]))
json.dump(info, open(os.path.join(HERE, "build_library_initial.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
