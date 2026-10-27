"""The library with the water and the aqueous ions in the SOAP (library_clayff_water_md.npz): the 11 bulk layered minerals and their
basal (0 0 1) slabs of the water-free library (../library/refs_bulk/, refs_basal/), plus the reference systems of
build_refs.py: each initial structure (<name>.data, a copy of its system.data: packmol's water and ions, not relaxed) and the frames of its 300 K MD up to
the second-to-last one (the last frame is kept out, for test.py).  The bulk crystals are kept whole, the slab rows pruned at
similarity 0.999 (as in the water-free library), the rows of the systems with water at PRUNE_WATER: liquid water has no two
identical environments, so its rows are thinned to the ones that differ (rows of another type or charge are never dropped).

    python build_library.py [prune_water]     (with the Python of soapff-gui)  ->  library_clayff_water_md.npz, build_library.json"""
import json
import os
import shutil
import sys

from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
LIBDIR = os.path.join(HERE, "..", "library")
BULK = ["brucite_bulk", "kaolinite_bulk", "lizerdite_bulk", "muscovite_1M_bulk", "muscovite_2M1_bulk", "nanpingite_bulk",
        "phlogopite_1M_bulk", "pyrophyllite_2_bulk", "pyrophyllite_bulk", "talc_bulk", "montmorillonite_dry_Na"]
BASAL = ["brucite", "kaolinite", "lizerdite", "muscovite_1M", "muscovite_2M1", "nanpingite", "phlogopite_1M", "pyrophyllite_2",
         "pyrophyllite", "talc", "montmorillonite_dry_Na"]
SYSTEMS = ["water", "NaCl", "KCl", "CsCl", "CaCl2", "BaCl2", "hydrate_2W"] + \
          ["if_" + m for m in ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
                               "pyrophyllite_2", "kaolinite", "lizerdite", "talc", "brucite"]]
PRUNE_WATER = float(sys.argv[1]) if len(sys.argv) > 1 else 0.99


def split_dump(path):
    """the frames of a LAMMPS dump as text blocks"""
    txt = open(path).read()
    return ["ITEM: TIMESTEP" + b for b in txt.split("ITEM: TIMESTEP")[1:]]


def training_frames(name):
    """refs/<name>/<name>_train.lammpstrj: every frame of the MD but the last (test.py's)"""
    d = os.path.join(HERE, "refs", name)
    fr = split_dump(os.path.join(d, "traj.lammpstrj"))
    out = os.path.join(d, name + "_train.lammpstrj")
    open(out, "w").write("".join(fr[:-1]))
    return out, len(fr) - 1


def reference(name):
    """refs/<name>/<name>.data: a copy of its system.data (the rows of the library are named after their file)"""
    d = os.path.join(HERE, "refs", name)
    out = os.path.join(d, name + ".data")
    shutil.copyfile(os.path.join(d, "system.data"), out)
    return out


if __name__ == "__main__":
    refs = [os.path.join(LIBDIR, "refs_bulk", b + ".data") for b in BULK] + \
           [os.path.join(LIBDIR, "refs_basal", m + "_001.data") for m in BASAL] + \
           [reference(s) for s in SYSTEMS]
    frames, info = [], {}
    for k, s in enumerate(SYSTEMS):
        f, n = training_frames(s)
        frames.append(f"{len(BULK) + len(BASAL) + k + 1}:{f}")
        info[s] = dict(atoms=core.read_lammps_data(refs[len(BULK) + len(BASAL) + k])[0].__len__(), training_frames=n)
    prune = [None] * len(BULK) + [0.999] * len(BASAL) + [PRUNE_WATER] * len(SYSTEMS)
    out = os.path.join(HERE, "library_clayff_water_md.npz")
    core.build_library(refs, frames, out=out, water_in_soap=True, prune=prune, protect=len(BULK))
    lib = core.load_library(out)
    r = core.check_library(out)
    json.dump(dict(prune_water=PRUNE_WATER, systems=info, rows=int(len(lib["q"])), species=lib["species"],
                   rows_per_type={t: int((lib["type"] == t).sum()) for t in sorted(set(lib["type"]))},
                   one_to_one=dict(rows=r["total"], same_charge=r["total"] - r["wrong"], confusions=r["confusions"])),
              open(os.path.join(HERE, "build_library.json"), "w"), indent=1)
