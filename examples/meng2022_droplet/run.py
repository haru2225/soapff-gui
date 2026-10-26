"""The initial structure of Meng et al., Appl. Clay Sci. 224, 106513 (2022): a NaCl droplet on montmorillonite (0 0 1)
(their Droplet-1: 9171 water, 112 Na+ and 112 Cl-, 3.83 wt%, diameter ~80 A, ~5 A above the surface, box 209.2 x 217.4 x 201 A)
made with soapff-gui (paper: test (6) and Fig. "Meng").
  paper/  as in the paper: the (0 0 1) slab of dry Na montmorillonite (pymatgen, the vacuum above it) with c set to 201 A,
          repeated 20 x 12 (206.9 x 215.5 A); water, Na+ and Cl- inside one sphere (packmol: inside sphere) holding 9171
          molecules at 1 g/cm3 (radius 40.3 A), its bottom 5 A above the surface; the force field (bundled library).
  gui/    the same with the GUI's steps only: Surface with the vacuum above the slab, 183 A thick (pymatgen rounds the cell up
          to whole layers: c = 204.4 A), Supercell 20 x 12 x 1, Add water with the same sphere.
The structures are large (105,977 atoms); packmol's xyz files are written but not kept in the repository.

    python run.py [paper|gui]   (with the Python of soapff-gui; packmol must be installed; ~5 minutes each)
                                ->  <variant>/ (structure.extxyz; the GUI's system.data, system.settings, in.lammps, charges.csv,
                                    report.json; packmol/), results.json"""
import json
import os
import shutil
import sys
import time

import ase.io
import numpy as np

from soapff_cui import core, surfaces, water as W

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "..", "library", "refs_bulk")
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")
N_WATER, N_SALT, GAP = 9171, 112, 5.0


def build(variant):
    t0 = time.time()
    d = os.path.join(HERE, variant)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    at, q, types = core.read_lammps_data(os.path.join(REFS, "montmorillonite_dry_Na.data"))
    if variant == "paper":
        sl, typ, inf = surfaces.make_slab(at, list(types), (0, 0, 1), min_thickness=16.0, vacuum=16.0, vacuum_side="above")
        sl.set_cell(np.vstack([sl.cell[0], sl.cell[1], sl.cell[2] / sl.cell.lengths()[2] * 201.0]), scale_atoms=False)
    else:
        sl, typ, inf = surfaces.make_slab(at, list(types), (0, 0, 1), min_thickness=16.0, vacuum=183.0, vacuum_side="above")
    slab = sl.repeat((20, 12, 1))
    fr = W.frame_of(slab)
    z1 = fr["z_atoms"][1]
    R = (N_WATER * W.M_WATER / W.AVOGADRO * 1e24 * 3 / (4 * np.pi)) ** (1 / 3)      # 9171 molecules at 1 g/cm3
    drop = [("inside", "sphere", [fr["ax"] / 2, fr["by"] / 2, z1 + GAP + R, R])]
    items = [dict(molecule="water", number=N_WATER, constraints=drop), dict(molecule="Na", number=N_SALT, constraints=drop),
             dict(molecule="Cl", number=N_SALT, constraints=drop)]
    new, rep = W.pack(slab, items, tolerance=2.0, seed=1, workdir=os.path.join(d, "packmol"), log=lambda *a: None)
    f = os.path.join(d, "structure.extxyz")
    ase.io.write(f, new, format="extxyz")
    t_pack = time.time() - t0
    r = core.make(f, d, LIB, md=False, stop_on_fallback=True, log=lambda *a: None)
    O = new.positions[len(slab):len(slab) + 3 * N_WATER:3]
    return dict(slab_atoms=len(slab), atoms=len(new), cell=[round(x, 2) for x in new.cell.cellpar()[:3]],
                sphere_radius=round(float(R), 2), packmol_converged=rep["converged"],
                placed={x["molecule"]: x["n"] for x in rep["items"]}, nearest_atom=min(x["min_distance"] for x in rep["items"]),
                droplet_bottom_above_surface=round(float(O[:, 2].min() - z1), 2), fallback=r["fallback"],
                type_check_differ=r["type_check"]["differ"], charge_sum=round(r["charge_sum"], 4),
                seconds_packmol=round(t_pack), seconds_total=round(time.time() - t0))


if __name__ == "__main__":
    f = os.path.join(HERE, "results.json")
    res = json.load(open(f)) if os.path.isfile(f) else {}
    for v in sys.argv[1:] or ["paper", "gui"]:
        res[v] = build(v)
        print(v, res[v], flush=True)
        json.dump(res, open(f, "w"), indent=1)
