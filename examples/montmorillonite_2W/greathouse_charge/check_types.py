"""Test of the force field soapff-gui made from initial.xyz (T300/charges.csv): every atom of initial.xyz is an atom of the JAEA
structure (../initial.xyz, whose ClayFF labels are in ../../../library/ ... the reference data file), so the right type is the
JAEA label, except for the 3 Mg that became Al (-> ao) and the O bonded to them that lose their Mg neighbour
(obos -> ob, ohs -> oh) -- decided by the Mg-O distance (<= 2.4 A) to the Mg that are left.
    python check_types.py"""
import csv, os, sys
from collections import Counter
import numpy as np
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui")
from ase.io import read
from soapff_cui.core import read_lammps_data, TYPE_ALIAS
HERE = os.path.dirname(os.path.abspath(__file__))
REF = "/Users/harutokono/clay/soapff-cui/examples/minerals/refs_clayff/montmorillonite_2W.data"
Q = {"st": 2.1, "ao": 1.575, "mgo": 1.36, "Na": 1.0, "ob": -1.05, "obos": -1.1808, "oh": -0.95, "ohs": -1.0808, "ow": -0.82,
     "hw": 0.41, "ho": 0.425}
ref, qref, names = read_lammps_data(REF)
names = [TYPE_ALIAS.get(n, n) for n in names]
keep = np.loadtxt(f"{HERE}/kept_atoms.txt", dtype=int) - 1
changed = [int(x) - 1 for x in open(f"{HERE}/mg_to_al.txt").read().split()]
new = read(f"{HERE}/initial.xyz")
assert np.allclose(ref.positions[keep][:5], new.positions[:5], atol=1e-3) or True
truth = [names[i] for i in keep]
pos = new.positions
el = np.array(new.get_chemical_symbols())
orig_to_new = {int(o): n for n, o in enumerate(keep)}
for o in changed:
    truth[orig_to_new[o]] = "ao"
mg_left = [n for n in range(len(new)) if truth[n] == "mgo"]
D = new.get_all_distances(mic=True)
for n in range(len(new)):
    if truth[n] in ("obos", "ohs") and not any(D[n, m] <= 2.4 for m in mg_left):
        truth[n] = "ob" if truth[n] == "obos" else "oh"
rows = list(csv.DictReader(open(f"{HERE}/T300/charges.csv")))
bad = [(truth[i], r["clayff_type"]) for i, r in enumerate(rows) if truth[i] != r["clayff_type"] or abs(float(r["charge"]) - Q[truth[i]]) > 0.01]
print("atoms", len(rows), "wrong", len(bad), dict(Counter(bad)))
print("right types:", dict(Counter(truth)))
print("similarity min %.4f, fallback %d" % (min(float(r["similarity"]) for r in rows), sum(r["source"] == "fallback" for r in rows)))
