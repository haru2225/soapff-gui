"""Library references with the plain ClayFF + flexible SPC labels (soapff_cui/clayff.py, Cygan et al. 2004 Tables 1, 2):
  - the Na-montmorillonite initial structure of JAEA (2W, 2,502 atoms; water placed artificially, overlapping the clay),
    relabelled: every atom takes the ClayFF charge of its type, water SPC (-0.82, +0.41)
  - the bulk crystals of quartz, albite, brucite and pyrophyllite (supercells >= 10 A, ../bulk/make_bulk.py), typed from
    the bond graph (mincharge ff_types), ClayFF charges.
Writes refs_clayff/<name>.data (LAMMPS data, atom_style full, type names in the Masses comments).
    python build_clayff_refs.py"""
import os, sys
import numpy as np
import ase.io
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "mincharge"))
from mincharge.ff import ff_types, write_lammps
from soapff_cui.clayff import Q
from soapff_cui.core import read_lammps_data, TYPE_ALIAS

JAEA = os.path.join(HERE, "..", "jaea", "initial.xyz")
JAEA_DATA = "/Users/harutokono/clay/test8/montmorillonite_sandwich_2W.lammps"
OUT = os.path.join(HERE, "refs_clayff"); os.makedirs(OUT, exist_ok=True)


def write(at, types, name, title):
    q = np.array([Q[t] for t in types])
    f = os.path.join(OUT, name + ".data")
    write_lammps(at, q, f, f.replace(".data", ".settings"), types=types, title=title, wrap=True)
    os.remove(f.replace(".data", ".settings"))                 # the force field is written by soapff-cui make
    print(name, len(at), "net charge %+.4f" % q.sum(), sorted(set(types)))


at, q0, names = read_lammps_data(JAEA_DATA)
write(at, [TYPE_ALIAS.get(n, n) for n in names], "montmorillonite_2W", "Na-montmorillonite 2W, JAEA initial structure, ClayFF + SPC")
for name in ("quartz", "albite", "brucite", "pyrophyllite"):
    at = ase.io.read(os.path.join(HERE, "..", "bulk", f"{name}_bulk.xyz"))
    write(at, ff_types(at), f"{name}_bulk", f"{name} bulk crystal, ClayFF")
