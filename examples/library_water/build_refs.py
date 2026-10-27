"""Reference systems with water and aqueous ions for the library with water in the SOAP (library_clayff_water.npz): structures
made with soapff-gui (packmol, soapff_cui.water) and its force field (the water-free library; water and ions typed by rule,
flexible SPC water, ClayFF Table 1 ions), then a 300 K MD with LAMMPS whose frames go into the library (build_library.py).

  water/           liquid water, 512 molecules in a 24.83 A cube
  NaCl/ KCl/ CsCl/ 8 cations + 8 Cl- in 496 water molecules (about 0.9 mol/L)
  CaCl2/ BaCl2/    6 cations + 12 Cl- in 494 water molecules (about 0.7 mol/L)
  if_<mineral>/    the library's basal (0 0 1) slab (../library/refs_basal/) with the vacuum above it filled with water (1 g/cm3,
                   packmol seed 2; the test structures of ../interfaces_001_water/ use seed 1)
  hydrate_2W/      the 2-layer Na-montmorillonite hydrate of ../montmorillonite_2W/ (hydrate_2W_start.data, the
                   system.data of its diffusion/start/ run); the library's refs_hydrate/montmorillonite_2W.data has water overlapping
                   the clay and cannot start an MD)

MD (in.lammps of soapff-gui): minimisation, 5 ps NVT relaxation, then the liquids and the hydrate 5 ps NVT + 20 ps NPT at 1 atm
(liquids isotropic, hydrate full cell), the interfaces 25 ps NVT; 0.5 fs, frames every 1 ps.

    python build_refs.py [lmp] [jobs] [system ...]     (with the Python of soapff-gui; default lmp, 3 jobs at a time;
                                                        lmp "none": only the inputs, run elsewhere, e.g. run_remote.py)"""
import concurrent.futures
import os
import subprocess
import sys
import tempfile

import numpy as np

import ase.io
from ase import Atoms

from soapff_cui import core, water

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")
REFS = os.path.join(HERE, "..", "library")
LMP = sys.argv[1] if len(sys.argv) > 1 else "lmp"
JOBS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
L = 24.83                                                    # cube edge (A): 512 water molecules at 1 g/cm3
SOLUTIONS = {"water": (512, []), "NaCl": (496, [("Na", 8), ("Cl", 8)]), "KCl": (496, [("K", 8), ("Cl", 8)]),
             "CsCl": (496, [("Cs", 8), ("Cl", 8)]), "CaCl2": (494, [("Ca", 6), ("Cl", 12)]),
             "BaCl2": (494, [("Ba", 6), ("Cl", 12)])}
MINERALS = ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
            "pyrophyllite_2", "kaolinite", "lizerdite", "talc", "brucite"]
quiet = lambda *a: None


def liquid(nw, ions, seed=2, tolerance=2.0):
    """water (+ ions) filling the periodic cube: packmol directly, periodic in x, y and z (water.pack needs a structure to place
    them around); the water geometry of soapff_cui.water (O-H 1.0 A, H-O-H 109.47 deg)"""
    h = np.radians(water.HOH / 2)
    mol = {"water": (["O", "H", "H"], [[0, 0, 0], [water.OH * np.sin(h), 0, water.OH * np.cos(h)],
                                       [-water.OH * np.sin(h), 0, water.OH * np.cos(h)]])}
    with tempfile.TemporaryDirectory() as d:
        inp = [f"tolerance {tolerance}", "filetype xyz", "output packed.xyz", f"seed {seed}", f"pbc 0. 0. 0. {L} {L} {L}"]
        for m, n in [("water", nw)] + list(ions):
            sym, pos = mol.get(m, ([m], [[0, 0, 0]]))
            water._xyz(os.path.join(d, m + ".xyz"), sym, np.array(pos, float))
            inp += ["", f"structure {m}.xyz", f"  number {n}", f"  inside box 0. 0. 0. {L} {L} {L}", "end structure"]
        open(os.path.join(d, "pack.inp"), "w").write("\n".join(inp) + "\n")
        subprocess.run("packmol < pack.inp > packmol.log", shell=True, cwd=d, check=True)
        if "Success!" not in open(os.path.join(d, "packmol.log")).read():
            raise RuntimeError(f"packmol did not converge for {ions}")
        sym, pos = water._read_xyz(os.path.join(d, "packed.xyz"))
    return Atoms(sym, positions=pos, cell=[L, L, L], pbc=True)


def structure(name):
    """the initial structure of a reference system and its MD settings for core.make"""
    if name in SOLUTIONS:
        return liquid(*SOLUTIONS[name]), dict(relax_ps=5, nvt_ps=5, npt_ps=20, baro="iso")
    if name == "hydrate_2W":
        return core.read_structure(os.path.join(HERE, "hydrate_2W_start.data")), \
            dict(relax_ps=5, nvt_ps=5, npt_ps=20, baro="tri")
    at, q, types = core.read_lammps_data(os.path.join(REFS, "refs_basal", name[3:] + "_001.data"))
    gap = water.gaps(at, min_gap=8.0, density=1.0)[0]                 # the vacuum above the slab
    s, rep = water.add_water(at, [(gap, gap["n_water"])], seed=2, log=quiet)
    return s, dict(relax_ps=5, nvt_ps=25, npt_ps=0)


def run(name):
    d = os.path.join(HERE, "refs", name)
    os.makedirs(d, exist_ok=True)
    s, md = structure(name)
    ase.io.write(os.path.join(d, "structure.extxyz"), s, format="extxyz")
    r = core.make(os.path.join(d, "structure.extxyz"), d, LIB, dump_ps=1.0, stop_on_fallback=True, log=quiet, **md)
    if LMP == "none":
        return name, r["atoms"], None
    p = subprocess.run(f"{LMP} -in in.lammps -log log.lammps -screen none", shell=True, cwd=d)
    return name, r["atoms"], p.returncode


if __name__ == "__main__":
    names = sys.argv[3:] or list(SOLUTIONS) + ["hydrate_2W"] + ["if_" + m for m in MINERALS]
    with concurrent.futures.ProcessPoolExecutor(JOBS) as ex:
        for name, n, rc in ex.map(run, names):
            print(name, n, "atoms, lammps rc", rc, flush=True)
