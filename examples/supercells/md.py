"""MD of the two supercells (paper: "MD with LAMMPS"): soapff-gui's force field and its MD input (minimisation, 10 ps NVT
relaxation, 10 ps NVT, 40 ps NPT with the full cell and tilt at 1 atm; 300 K, 0.5 fs, frames every 1 ps), run with LAMMPS.
Reported from log.lammps and the frames: the temperature range of the production, the mean temperature and pressure of the
NVT (the experimental box), and the volume averaged over the last 20 ps of the NPT against the input box.

    python md.py [mpirun lmp np]   (with the Python of soapff-gui; default: mpirun lmp 2)  ->  md_<name>/, md_results.json"""
import json
import os
import subprocess
import sys

import numpy as np
from ase.io import read

from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")
MPIRUN, LMP, NP = (sys.argv[1:4] + ["mpirun", "lmp", "2"][len(sys.argv[1:4]):])


def thermo(log):
    """the thermo rows of the runs after the minimisation: [(run index, step, temp, press, vol)]"""
    rows, run, head = [], -1, None
    for l in open(log):
        p = l.split()
        if p[:2] == ["Step", "Temp"]:
            run, head = run + 1, p
            continue
        if head and len(p) == len(head) and p[0].isdigit():
            rows.append((run, int(p[0]), float(p[1]), float(p[2]), float(p[head.index("Volume")])))
        elif head and p[:1] == ["Loop"]:
            head = None
    return rows


res = {}
for name in ("brucite_5x5x4", "pyrophyllite_3x2x3"):
    d = os.path.join(HERE, "md_" + name)
    at = read(os.path.join(HERE, name + ".extxyz"))
    core.make(os.path.join(HERE, name + ".extxyz"), d, LIB, nvt_ps=10, npt_ps=40, baro="tri", relax_ps=10, dump_ps=1.0,
              log=lambda *a: None)
    subprocess.run(f"{MPIRUN} -np {NP} {LMP} -in in.lammps -log log.lammps -screen none", shell=True, cwd=d, check=True)
    rows = thermo(os.path.join(d, "log.lammps"))
    runs = sorted({r[0] for r in rows})                      # minimisations (2), relaxation NVT, NVT, NPT
    relax, nvt, npt = runs[-3], runs[-2], runs[-1]
    T = [r[2] for r in rows if r[0] in (nvt, npt)]
    nvt_rows = [r for r in rows if r[0] == nvt]
    frames = core.read_dump_frames(os.path.join(d, "traj.lammpstrj"), at.get_chemical_symbols())   # 1 ps: NVT 0-10, NPT 10-50
    v_npt = np.mean([f.get_volume() for k, f in enumerate(frames) if k >= 30])
    res[name] = dict(temperature_range=[round(min(T)), round(max(T))], nvt_mean_temperature=round(float(np.mean([r[2] for r in nvt_rows])), 1),
                     nvt_mean_pressure_GPa=round(float(np.mean([r[3] for r in nvt_rows])) * 1.01325e-4, 2),
                     npt_last20ps_volume_change_percent=round(float((v_npt / at.get_volume() - 1) * 100), 1))
    print(name, res[name], flush=True)
json.dump(res, open(os.path.join(HERE, "md_results.json"), "w"), indent=1)
