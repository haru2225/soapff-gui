"""Initial structures on a basal (0 0 1) surface made with packmol's constraints (paper: test (5); Fig. "shapes" shows the
last one, the water film with Na+ and Cl-).  The (0 0 1) slab of dry Na montmorillonite or muscovite 2M1 as the GUI's Surface
cuts it (pymatgen, >= 16 A of mineral, 40 A of vacuum above it), repeated 3 x 2 in the a-b plane; the water (and ions) placed
by soapff_cui.water.pack (the GUI's "Add water", tolerance 2 A, seed 1) in a region given by constraints in packmol's frame
(x' = x - y cot(gamma), z normal to the slab; z1 = the top of the slab):
  sphere      a droplet: inside sphere (radius 10 A, centre 13 A above z1)
  hemisphere  a droplet on the surface: inside sphere (radius 13 A, centre 2 A above z1) + over plane z = z1 + 2 (muscovite)
  cylinder    a water column from the surface to the bottom of the slab above (periodic image): inside cylinder, radius 8 A
  cavity      a water film (2-22 A above z1) with a spherical cavity: inside box + outside sphere (radius 7 A)
  xygauss     a water film with a mound: inside box + below xygauss
  ions        a water film (2-22 A above z1) with 8 Na+ and 8 Cl-: inside box (muscovite)
The number of water molecules: 1 g/cm3 in the free volume of the region (water.estimate, the GUI's "Number from density").
Then the force field of soapff-gui (bundled library) and a LAMMPS minimisation (Ewald 1e-4, as in ../interfaces_001_water).

    python run.py [lmp]      (with the Python of soapff-gui; packmol must be installed)
                              ->  <case>/ (structure.extxyz; the GUI's system.data, system.settings, in.lammps, charges.csv,
                                  report.json; packmol/; the minimisation in.min, log.min), results.json"""
import json
import os
import shutil
import subprocess
import sys

import ase.io

from soapff_cui import core, surfaces, water as W

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "..", "library", "refs_bulk")
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")
LMP = sys.argv[1] if len(sys.argv) > 1 else "lmp_serial"
EWALD = "kspace_style ewald 1.0e-4  # this test: Ewald 1e-4 instead of the PPPM of system.settings\n"


def slab_of(name, reps):
    at, q, types = core.read_lammps_data(os.path.join(REFS, name + ".data"))
    sl, typ, inf = surfaces.make_slab(at, list(types), (0, 0, 1), min_thickness=16.0, vacuum=40.0, vacuum_side="above")
    return sl.repeat(reps)


def case(slab, kind):
    """the items of water.pack for one kind of initial structure"""
    fr = W.frame_of(slab)
    ax, by, cz = fr["ax"], fr["by"], fr["cz"]
    z0, z1 = fr["z_atoms"]
    cx, cy = ax / 2, by / 2
    film = ("inside", "box", [0.0, 0.0, z1 + 2.0, ax, by, z1 + 22.0])
    if kind == "sphere":
        cons = [("inside", "sphere", [cx, cy, z1 + 13.0, 10.0])]
    elif kind == "hemisphere":
        cons = [("inside", "sphere", [cx, cy, z1 + 2.0, 13.0]), ("over", "plane", [0.0, 0.0, 1.0, z1 + 2.0])]
    elif kind == "cylinder":
        cons = [("inside", "cylinder", [cx, cy, z1 + 2.0, 0.0, 0.0, 1.0, 8.0, cz + z0 - z1 - 4.0])]
    elif kind == "cavity":
        cons = [film, ("outside", "sphere", [cx, cy, z1 + 12.0, 7.0])]
    elif kind == "xygauss":
        cons = [("inside", "box", [0.0, 0.0, z1 + 2.0, ax, by, z1 + 26.0]),
                ("below", "xygauss", [cx, cy, 6.0, 6.0, z1 + 10.0, 14.0])]
    else:
        cons = [film]
    items = [dict(molecule="water", number=W.estimate(slab, cons, density=1.0)["n_water"], constraints=cons)]
    if kind == "ions":
        items += [dict(molecule="Na", number=8, constraints=cons), dict(molecule="Cl", number=8, constraints=cons)]
    return items


def minimize(d, lmp=LMP):
    """in.min: the GUI's in.lammps with Ewald 1e-4 after the include (system.settings is kept as written)"""
    s = open(os.path.join(d, "in.lammps")).read().replace("include system.settings\n", "include system.settings\n" + EWALD, 1)
    open(os.path.join(d, "in.min"), "w").write(s)
    p = subprocess.run([lmp, "-in", "in.min", "-log", "log.min"], cwd=d, capture_output=True, text=True)
    L = open(os.path.join(d, "log.min")).read().splitlines()
    k = [i for i, l in enumerate(L) if "Energy initial, next-to-last, final" in l]
    return dict(returncode=p.returncode, stop=[l.split("=")[1].strip() for l in L if "Stopping criterion" in l],
                energy_initial_final=[float(L[k[0] + 1].split()[0]), float(L[k[-1] + 1].split()[2])] if k else None)


if __name__ == "__main__":
    mmt, musc = slab_of("montmorillonite_dry_Na", (3, 2, 1)), slab_of("muscovite_2M1_bulk", (3, 2, 1))
    res = {}
    for kind, slab in (("sphere", mmt), ("hemisphere", musc), ("cylinder", mmt), ("cavity", mmt), ("xygauss", mmt),
                       ("ions", musc)):
        d = os.path.join(HERE, kind)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        items = case(slab, kind)
        s, rep = W.pack(slab, items, tolerance=2.0, seed=1, workdir=os.path.join(d, "packmol"), log=lambda *a: None)
        ase.io.write(os.path.join(d, "structure.extxyz"), s, format="extxyz")
        r = core.make(os.path.join(d, "structure.extxyz"), d, LIB, md=False, stop_on_fallback=True, log=lambda *a: None)
        res[kind] = dict(atoms=len(s), packmol_converged=rep["converged"], placed={x["molecule"]: x["n"] for x in rep["items"]},
                         nearest_atom=min(x["min_distance"] for x in rep["items"]), fallback=r["fallback"],
                         type_check_differ=r["type_check"]["differ"], charge_sum=round(r["charge_sum"], 4),
                         constraints=items[0]["constraints"], minimization=minimize(d))
        print(kind, {k: v for k, v in res[kind].items() if k != "constraints"}, flush=True)
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
