"""The 8 sheet silicates of the paper (Table "sheet silicates"): muscovite 2M1 and 1M, phlogopite 1M, nanpingite, kaolinite,
lizardite, talc, pyrophyllite (second source).
  (1) <mineral>.extxyz, a supercell one unit cell longer along a than the library reference (elements, coordinates, box);
  (2) the last frame of a 300 K MD with soapff-gui's force field (md_<mineral>/: minimisation, 10 ps NVT relaxation, 10 ps
      NVT, 40 ps NPT with the full cell and tilt at 1 atm, dump every 1 ps; traj.lammpstrj, log.lammps), the same atoms.
Both are typed with the bundled library and compared with the truth: the bond-graph ClayFF type of the supercell (same atom
order) and its Table 1 charge.  The lattice parameters averaged over the last 30 ps of the NPT, per unit cell, are compared
with the input structure; from log.lammps: the temperature range of the production and the mean temperature of the NVT; the
layer spacing d001 of the two muscovites (the unit-cell height normal to the a-b face per layer: 2 layers in 2M1, 1 in 1M).

    python run.py            (with the Python of soapff-gui)  ->  out_<mineral>/, out_<mineral>_md/ (structure.extxyz and the
                                                                GUI's files for it), results.json"""
import csv
import json
import os

import numpy as np
from ase.io import read

from soapff_cui import core
from soapff_cui._clayff.ff import ff_types
from soapff_cui.clayff import Q

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")
LAYERS = {"muscovite_2M1": 2, "muscovite_1M": 1}                # layers per unit cell along c
MINERALS = ["muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "kaolinite", "lizerdite", "talc", "pyrophyllite_2"]


def compare(at, truth, out):
    rep = core.make(at, out, LIB, md=False, stop_on_fallback=True, log=lambda *a: None)
    at.write(os.path.join(out, "structure.extxyz"), format="extxyz")
    rows = list(csv.DictReader(open(os.path.join(out, "charges.csv"))))
    wrong = [i + 1 for i, (r, t) in enumerate(zip(rows, truth)) if r["clayff_type"] != t or abs(float(r["charge"]) - Q[t]) > 0.01]
    return dict(atoms=len(rows), wrong=len(wrong), fallback=rep["fallback"],
                min_similarity=round(min(float(r["similarity"]) for r in rows), 3), types=sorted(set(truth)))


def thermo(log):
    """the thermo rows of log.lammps: [(run index, step, temp, press)]"""
    rows, run, head = [], -1, None
    for l in open(log):
        p = l.split()
        if p[:2] == ["Step", "Temp"]:
            run, head = run + 1, p
            continue
        if head and len(p) == len(head) and p[0].isdigit():
            rows.append((run, int(p[0]), float(p[1]), float(p[2])))
        elif head and p[:1] == ["Loop"]:
            head = None
    return rows


def volume(p):
    c = np.cos(np.radians(p[3:]))
    return np.prod(p[:3]) * np.sqrt(1 - (c ** 2).sum() + 2 * np.prod(c))


res = {}
for m in MINERALS:
    test = read(os.path.join(HERE, m + ".extxyz"))
    truth = ff_types(test)
    r = dict(supercell=compare(test, truth, os.path.join(HERE, f"out_{m}")))
    frames = core.read_dump_frames(os.path.join(HERE, f"md_{m}", "traj.lammpstrj"), test.get_chemical_symbols())
    r["md_last_frame"] = compare(frames[-1], truth, os.path.join(HERE, f"out_{m}_md"))
    sc = json.load(open(os.path.join(HERE, f"md_{m}", "supercell.json")))
    rep_ = np.array(sc["repeat"], float)
    cp = np.array([f.cell.cellpar() for k, f in enumerate(frames) if k >= 20])      # 1 ps per frame: the last 30 ps of the NPT
    mean, ref = cp.mean(0), np.array(sc["cell0"])
    mean[:3] /= rep_; ref[:3] /= rep_
    r["lattice"] = dict(md=mean.round(3).tolist(), input=ref.round(3).tolist(),
                        abc_dev_percent=((mean[:3] - ref[:3]) / ref[:3] * 100).round(2).tolist(),
                        angle_dev_deg=(mean[3:] - ref[3:]).round(2).tolist(),
                        volume_dev_percent=round(float((volume(mean) / volume(ref) - 1) * 100), 2))
    rows = thermo(os.path.join(HERE, f"md_{m}", "log.lammps"))
    runs = sorted({x[0] for x in rows})                     # minimisations, relaxation NVT, NVT, NPT
    T = [x[2] for x in rows if x[0] in runs[-2:]]
    r["md"] = dict(temperature_range=[round(min(T)), round(max(T))],
                   nvt_mean_temperature=round(float(np.mean([x[2] for x in rows if x[0] == runs[-2]])), 1))
    if m in LAYERS:                                         # d001: the cell height normal to the a-b face per layer
        r["md"]["d001_input_md"] = [round(float(volume(c) / (c[0] * c[1] * np.sin(np.radians(c[5]))) / LAYERS[m]), 2)
                                    for c in (ref, mean)]
    res[m] = r
    print(m, r["supercell"]["atoms"], "wrong", r["supercell"]["wrong"], "| MD frame wrong", r["md_last_frame"]["wrong"],
          "min sim", r["md_last_frame"]["min_similarity"], "fallback", r["md_last_frame"]["fallback"],
          "| V", r["lattice"]["volume_dev_percent"], "%", flush=True)
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
