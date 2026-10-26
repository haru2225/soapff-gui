"""Crystals in supercells that differ from the library's (paper: Table "larger crystals" and Fig. 1, the workflow): brucite
5x5x4 (500 atoms) and pyrophyllite 3x2x3 (720 atoms), elements, coordinates and box only.  soapff-gui's force field (the
bundled library) is compared with the truth: the ClayFF type from the bond graph (soapff_cui._clayff.ff.ff_types) and its
ClayFF Table 1 charge.  "wrong": an atom whose type differs or whose charge differs by more than 0.01 e.

    python run.py            (with the Python of soapff-gui)  ->  out_<name>/ (structure.extxyz and the GUI's files: system.data,
                                                                system.settings, in.lammps, charges.csv, report.json), results.json"""
import csv
import json
import os
import shutil

from ase.io import read

from soapff_cui import core
from soapff_cui._clayff.ff import ff_types
from soapff_cui.clayff import Q

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")


def compare(structure, out):
    """soapff-gui's types and charges for `structure` (written to `out`) against the bond-graph ClayFF types"""
    at = read(structure)
    truth = ff_types(at)
    rep = core.make(structure, out, LIB, md=False, stop_on_fallback=True, log=lambda *a: None)
    shutil.copy(structure, os.path.join(out, "structure.extxyz"))
    rows = list(csv.DictReader(open(os.path.join(out, "charges.csv"))))
    wrong = [i + 1 for i, (r, t) in enumerate(zip(rows, truth))
             if r["clayff_type"] != t or abs(float(r["charge"]) - Q[t]) > 0.01]
    return dict(atoms=len(rows), wrong=len(wrong), wrong_atoms=wrong[:20], fallback=rep["fallback"],
                min_similarity=round(min(float(r["similarity"]) for r in rows), 3), type_check_differ=rep["type_check"]["differ"],
                consistency=rep["consistency"], charge_sum=round(rep["charge_sum"], 4))


if __name__ == "__main__":
    res = {}
    for name in ("brucite_5x5x4", "pyrophyllite_3x2x3"):
        res[name] = compare(os.path.join(HERE, name + ".extxyz"), os.path.join(HERE, "out_" + name))
        print(name, res[name], flush=True)
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
