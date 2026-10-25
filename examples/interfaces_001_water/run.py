"""Water on the basal (0 0 1) slabs (paper: test (4) and Fig. "interfaces"): the library slab of each of the 10 layered
minerals of the figure (../library/refs_basal/; the second pyrophyllite source is left out, as in the paper) with its
widest empty slice filled with water by packmol at 1 g/cm3 (soapff_cui.water, the GUI's "Add water": tolerance 2 A, seed 1),
the force field of soapff-gui (bundled library, flexible SPC water) and a LAMMPS minimisation.  The minimisation uses Ewald
(1e-4) instead of PPPM, as in the paper (PPPM stopped at the start of the minimisation on the LAMMPS used there).

    python run.py [lmp]      (with the Python of soapff-gui; lmp: the LAMMPS executable, default lmp_serial)
                              ->  <mineral>/ (structure.extxyz, system.data, system.settings, in.lammps, log.min), results.json"""
import json
import os
import shutil
import subprocess
import sys

import ase.io

from soapff_cui import core, water

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "..", "library", "refs_basal")
LIB = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz")
LMP = sys.argv[1] if len(sys.argv) > 1 else "lmp_serial"
MINERALS = ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
            "kaolinite", "lizerdite", "talc", "brucite"]


def minimize(d, lmp=LMP):
    """the minimisation of in.lammps (make md=False) with Ewald 1e-4; returns the stopping criteria and the energies"""
    s = open(os.path.join(d, "system.settings")).read().replace("kspace_style pppm 1.0e-5", "kspace_style ewald 1.0e-4")
    open(os.path.join(d, "system.settings"), "w").write(s)
    p = subprocess.run([lmp, "-in", "in.lammps", "-log", "log.min"], cwd=d, capture_output=True, text=True)
    L = open(os.path.join(d, "log.min")).read().splitlines()
    k = [i for i, l in enumerate(L) if "Energy initial, next-to-last, final" in l]
    return dict(returncode=p.returncode, stop=[l.split("=")[1].strip() for l in L if "Stopping criterion" in l],
                energy_initial_final=[float(L[k[0] + 1].split()[0]), float(L[k[-1] + 1].split()[2])] if k else None)


if __name__ == "__main__":
    res = {}
    for m in MINERALS:
        d = os.path.join(HERE, m)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        at, q, types = core.read_lammps_data(os.path.join(REFS, m + "_001.data"))
        gap = water.gaps(at, min_gap=8.0, density=1.0)[0]                       # the vacuum above the slab
        s, rep = water.add_water(at, [(gap, gap["n_water"])], log=lambda *a: None)
        ase.io.write(os.path.join(d, "structure.extxyz"), s, format="extxyz")
        r = core.make(os.path.join(d, "structure.extxyz"), d, LIB, md=False, stop_on_fallback=True, log=lambda *a: None)
        res[m] = dict(slab_atoms=len(at), water=rep[0]["n"], atoms=r["atoms"], packmol_converged=rep[0]["converged"],
                      nearest_atom=rep[0]["min_distance"], fallback=r["fallback"], type_check_differ=r["type_check"]["differ"],
                      charge_sum=round(r["charge_sum"], 4), minimization=minimize(d))
        print(m, res[m], flush=True)
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
