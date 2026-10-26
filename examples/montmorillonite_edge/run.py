"""The optional library "Montmorillonite edge (Underwood & Bourg 2020)" on the montmorillonite particle with edges on every
side (paper: test (7)).  The particle files are in ../library/refs_edge/ (clay_montmorillonite.data, .settings and
clay_montmorillonite_2.data, the particle twice along c).
  1. the particle (2,292 atoms) and the particle twice (4,584 atoms) as elements, coordinates and box only -> soapff-gui's force
     field with the edge library (src/soapff_gui/data/library_clayff_edge.npz); types and charges against the files, the
     pair, bond and angle coefficients against the settings, the numbers of bonds and angles;
  2. with the default library (ClayFF) the particle stops: its edge atoms are not in that library;
  3. a LAMMPS single point (Ewald 1e-5, cutoff 10 A, special_bonds 0 0 0) of the written system.data / system.settings against
     the original files (the original data file needs a blank line after "Angle Coeffs" to be read by LAMMPS: a copy is fixed);
  4. minimisations (Ewald 1e-4): the particle, the particle twice, and the particle with a water film (1 g/cm3, 2-12 A above it)
     and 4 Na+ + 4 Cl- (packmol): the water and Cl- take ClayFF (flexible SPC), the particle and Na+ the edge model.

    python run.py [lmp]      (with the Python of soapff-gui; lmp default lmp_serial)
                              ->  out_*/ (structure.extxyz; the GUI's system.data, system.settings, in.lammps, charges.csv,
                                  report.json; the minimisation in.min, log.min), single_point/, results.json"""
import csv
import json
import os
import shutil
import subprocess
import sys

import ase.io
from ase import Atoms

from soapff_cui import core, water as W

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "..", "library", "refs_edge")
DATA = os.path.join(HERE, "..", "..", "src", "soapff_gui", "data")
EDGE, BASIC = os.path.join(DATA, "library_clayff_edge.npz"), os.path.join(DATA, "library_clayff_bulk.npz")
LMP = sys.argv[1] if len(sys.argv) > 1 else "lmp_serial"
EWALD = "kspace_style ewald 1.0e-4  # this test: Ewald 1e-4 instead of the PPPM of system.settings\n"
quiet = lambda *a: None


def bare(name):
    """the elements, coordinates and box of a particle file -> <name>.extxyz; returns (path, charges, types)"""
    at, q, t = core.read_lammps_data(os.path.join(REFS, name + ".data"))
    f = os.path.join(HERE, name + ".extxyz")
    ase.io.write(f, Atoms(at.get_chemical_symbols(), positions=at.positions, cell=at.cell, pbc=True), format="extxyz")
    return f, q, t


def lammps(d, inp, log):
    p = subprocess.run([LMP, "-in", inp, "-log", log], cwd=d, capture_output=True, text=True)
    return p.returncode, open(os.path.join(d, log)).read().splitlines()


def minimize(d):
    """in.min: the GUI's in.lammps with Ewald 1e-4 after the include (system.settings is kept as written)"""
    s = open(os.path.join(d, "in.lammps")).read().replace("include system.settings\n", "include system.settings\n" + EWALD, 1)
    open(os.path.join(d, "in.min"), "w").write(s)
    rc, L = lammps(d, "in.min", "log.min")
    k = [i for i, l in enumerate(L) if "Energy initial, next-to-last, final" in l]
    return dict(returncode=rc, stop=[l.split("=")[1].strip() for l in L if "Stopping criterion" in l],
                energy_initial_final=[float(L[k[0] + 1].split()[0]), float(L[k[-1] + 1].split()[2])] if k else None)


res = {}
ff = core.load_library(EDGE)["forcefield"]
# 1. the particle and the particle twice with the edge library
for name in ("clay_montmorillonite", "clay_montmorillonite_2"):
    f, q, t = bare(name)
    d = os.path.join(HERE, "out_" + name)
    shutil.rmtree(d, ignore_errors=True)
    rep = core.make(f, d, EDGE, md=False, stop_on_fallback=True, log=quiet)
    shutil.copy(f, os.path.join(d, "structure.extxyz"))
    rows = list(csv.DictReader(open(os.path.join(d, "charges.csv"))))
    D = open(os.path.join(d, "system.data")).read().splitlines()
    S = open(os.path.join(d, "system.settings")).read().splitlines()
    tid = {l.split("#")[1].strip(): int(l.split()[0]) for l in core._section(D, "Masses", raw=True)}
    inv = {v: k for k, v in tid.items()}
    lj_differ = [l for l in S if l.startswith("pair_coeff") and [float(x) for x in l.split()[3:5]] != ff["lj"][inv[int(l.split()[1])]]]
    res[name] = dict(atoms=len(rows), wrong_type=sum(r["clayff_type"] != x for r, x in zip(rows, t)),
                     wrong_charge=int(sum(abs(float(r["charge"]) - x) > 1e-6 for r, x in zip(rows, q))), fallback=rep["fallback"],
                     min_similarity=round(min(float(r["similarity"]) for r in rows), 4), type_check_differ=rep["type_check"]["differ"],
                     pair_coeff_differ=lj_differ, bonds=next(l for l in D if l.endswith(" bonds")),
                     angles=next(l for l in D if l.endswith(" angles")),
                     bonded_coeffs=[l for l in S if l.startswith(("bond_coeff", "angle_coeff"))])
    res[name]["minimization"] = minimize(d)
    print(name, res[name], flush=True)
# 2. the default library
try:
    core.make(os.path.join(HERE, "clay_montmorillonite.extxyz"), os.path.join(HERE, "out_default_library"), BASIC, md=False,
              stop_on_fallback=True, log=quiet)
    res["default_library"] = "no fallback"
except core.FallbackError as e:
    res["default_library"] = str(e)[:120]
shutil.rmtree(os.path.join(HERE, "out_default_library"), ignore_errors=True)
print("default library:", res["default_library"], flush=True)
# 3. single point: the original files against the written ones
sp = os.path.join(HERE, "single_point")
shutil.rmtree(sp, ignore_errors=True)
os.makedirs(sp)
s = open(os.path.join(REFS, "clay_montmorillonite.data")).read().replace("Angle Coeffs # harmonic\n1 ", "Angle Coeffs # harmonic\n\n1 ", 1)
open(os.path.join(sp, "original.data"), "w").write(s)
shutil.copy(os.path.join(REFS, "clay_montmorillonite.settings"), os.path.join(sp, "original.settings"))
shutil.copy(os.path.join(HERE, "out_clay_montmorillonite", "system.data"), os.path.join(sp, "system.data"))
s = open(os.path.join(HERE, "out_clay_montmorillonite", "system.settings")).read().replace("kspace_style pppm 1.0e-5", "kspace_style ewald 1.0e-5")
open(os.path.join(sp, "written.settings"), "w").write(s)
tail = "thermo_style custom pe ebond eangle evdwl ecoul elong\nthermo_modify format float %.6f\nrun 0\n"
open(os.path.join(sp, "in.original"), "w").write(
    "units real\natom_style full\nbond_style harmonic\nangle_style harmonic\npair_style lj/cut/coul/long 10.0\n"
    "special_bonds lj/coul 0.0 0.0 0.0\nread_data original.data\ninclude original.settings\n" + tail)
open(os.path.join(sp, "in.written"), "w").write("include written.settings\n" + tail)
res["single_point"] = {}
for x in ("original", "written"):
    rc, L = lammps(sp, "in." + x, "log." + x)
    k = next(i for i, l in enumerate(L) if l.split()[:1] == ["PotEng"])
    res["single_point"][x] = dict(zip(("pe", "ebond", "eangle", "evdwl", "ecoul", "elong"), map(float, L[k + 1].split())))
print("single point", res["single_point"], flush=True)
# 4. the particle with a water film and NaCl
at = ase.io.read(os.path.join(HERE, "clay_montmorillonite.extxyz"))
fr = W.frame_of(at)
cons = [("inside", "box", [0.0, 0.0, fr["z_atoms"][1] + 2.0, fr["ax"], fr["by"], fr["z_atoms"][1] + 12.0])]
items = [dict(molecule="water", number=W.estimate(at, cons, density=1.0)["n_water"], constraints=cons),
         dict(molecule="Na", number=4, constraints=cons), dict(molecule="Cl", number=4, constraints=cons)]
s, rp = W.pack(at, items, seed=1, log=quiet)
f = os.path.join(HERE, "particle_water.extxyz")
ase.io.write(f, s, format="extxyz")
d = os.path.join(HERE, "out_particle_water")
shutil.rmtree(d, ignore_errors=True)
rep = core.make(f, d, EDGE, md=False, stop_on_fallback=True, log=quiet)
shutil.copy(f, os.path.join(d, "structure.extxyz"))
res["particle_water"] = dict(atoms=rep["atoms"], placed={x["molecule"]: x["n"] for x in rp["items"]}, roles=rep["roles"],
                             fallback=rep["fallback"], type_check_differ=rep["type_check"]["differ"],
                             charge_sum=round(rep["charge_sum"], 4), minimization=minimize(d))
print("particle + water", res["particle_water"], flush=True)
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
