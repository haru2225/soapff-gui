"""Library references of the bulk sheet silicates (../bulk/make_bulk_sheet.py), typed from the bond graph (mincharge
ff_types) with the plain ClayFF charges (Cygan et al. 2004, Table 1), as build_clayff_refs.py does for the first five.
Writes refs_clayff/<name>_bulk.data.
    python build_clayff_refs_sheet.py"""
import os, sys
from collections import Counter
import numpy as np
import ase.io
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "mincharge"))
from mincharge.ff import ff_types, write_lammps                 # noqa: E402
from soapff_cui.clayff import Q                                  # noqa: E402
OUT = os.path.join(HERE, "refs_clayff")


def write(at, types, name, title):                               # as build_clayff_refs.write
    q = np.array([Q[t] for t in types])
    f = os.path.join(OUT, name + ".data")
    write_lammps(at, q, f, f.replace(".data", ".settings"), types=types, title=title, wrap=True)
    os.remove(f.replace(".data", ".settings"))
    print(name, len(at), "net charge %+.4f" % q.sum(), sorted(set(types)))

NAMES = ("muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "kaolinite", "lizerdite", "talc", "pyrophyllite_2")
for name in NAMES:
    at = ase.io.read(os.path.join(HERE, "..", "bulk", f"{name}_bulk.xyz"))
    types = ff_types(at)
    write(at, types, f"{name}_bulk", f"{name} bulk crystal, ClayFF")
    print("   ", dict(Counter(types)))
