"""The solution library of the edge library (library_clayff_edge_water.npz) on structures it has not seen: the edge particle (1 layer, 2,292 atoms) and the
particle twice along c (4,584 atoms), each with a packmol water film 2-12 A above it (1 g/cm3), 4 or 8 Na+ and Cl-, with other packmol seeds than the
reference of the library (../montmorillonite_edge/particle_water.extxyz, seed 1).  Not relaxed, no MD frame.  Truth: the particle types of the dry
particle (the edge library on the structure without water), ow/hw for the water, the element for the added ions.
  (1) soapff-gui as it runs (runner.prepare: the mineral first, then the water and the ions with the solution library): wrong types, fallback
  (2) the SOAP of the whole structure against the water-free edge library: atoms of the particle below 0.95 and wrong
  (3) the SOAP of the whole structure against the solution library: wrong, below 0.95, minimum similarity per kind
    python test_edge_water.py     (with the Python of soapff-gui)  ->  test_edge_water.json"""
import csv
import json
import os
import tempfile

import ase.io
import numpy as np

from soapff_cui import core, water as W
from soapff_gui import runner
from water_common import nearest

HERE = os.path.dirname(os.path.abspath(__file__))
ME = os.path.join(HERE, "..", "montmorillonite_edge")
LIBM = core.load_library(runner.LIBRARY_EDGE)
LIBS = core.load_library(runner.LIBRARY_EDGE_WATER)
quiet = lambda *a: None
THR = 0.95
res = {}
for name, f, seed, nion in [("particle_seed2", "clay_montmorillonite.extxyz", 2, 4), ("particle_twice_seed3", "clay_montmorillonite_2.extxyz", 3, 8)]:
    dry = ase.io.read(os.path.join(ME, f))
    d0 = tempfile.mkdtemp()
    core.make(dry, d0, runner.LIBRARY_EDGE, md=False, stop_on_fallback=True, log=quiet)
    tdry = np.array([r["clayff_type"] for r in csv.DictReader(open(os.path.join(d0, "charges.csv")))])
    fr = W.frame_of(dry)
    cons = [("inside", "box", [0.0, 0.0, fr["z_atoms"][1] + 2.0, fr["ax"], fr["by"], fr["z_atoms"][1] + 12.0])]
    items = [dict(molecule="water", number=W.estimate(dry, cons, density=1.0)["n_water"], constraints=cons),
             dict(molecule="Na", number=nion, constraints=cons), dict(molecule="Cl", number=nion, constraints=cons)]
    s, rp = W.pack(dry, items, seed=seed, log=quiet)
    n0 = len(dry)
    sym = np.array(s.get_chemical_symbols())
    added = np.asarray(s.arrays["added"], bool)
    truth = np.array(list(tdry) + ["x"] * (len(s) - n0), dtype=object)
    mol = np.nonzero(added)[0]
    for i in mol:
        truth[i] = ("ow" if sym[i] == "O" else "hw") if sym[i] in ("O", "H") else sym[i]
    assert not (np.asarray(s.arrays["added"], bool)[:n0]).any()
    # (1) soapff-gui
    d1 = tempfile.mkdtemp()
    rep = runner.prepare(s, d1, runner.LIBRARY_EDGE, log=quiet, md=False)
    t1 = np.array([r["clayff_type"] for r in csv.DictReader(open(os.path.join(d1, "charges.csv")))], dtype=object)
    # (2) whole structure against the water-free edge library: the particle atoms
    t2, s2 = nearest(s, LIBM)
    pm = np.arange(len(s)) < n0
    # (3) whole structure against the solution library
    t3, s3 = nearest(s, LIBS)
    water_ = np.isin(truth, ["ow", "hw"]); ion = np.isin(truth, ["Na", "Cl"]) & added
    out = dict(atoms=len(s), particle_atoms=n0, water=int((truth == "ow").sum()), ions_added=int(ion.sum()),
               prepare=dict(wrong=int((t1 != truth).sum()), fallback=rep["fallback"], type_check_differ=rep["type_check"]["differ"], charge_sum=round(rep["charge_sum"], 4)),
               whole_vs_water_free=dict(low=int((s2[pm] < THR).sum()), wrong=int((t2[pm] != truth[pm]).sum()), min=float(s2[pm].min())),
               whole_vs_solution=dict(wrong=int((t3 != truth).sum()), low=int((s3 < THR).sum()), min_water=float(s3[water_].min()),
                                      min_ion=float(s3[ion].min()), min_particle=float(s3[pm].min())))
    res[name] = out
    print(name, out, flush=True)
json.dump(res, open(os.path.join(HERE, "test_edge_water.json"), "w"), indent=1)
