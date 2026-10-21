"""Na-montmorillonite of Cygan et al. (2004), Na3(Si31Al)(Al14Mg2)O80(OH)16, built like theirs: pyrophyllite (the
dioctahedral 2:1 layer of the soapff library, COD 9000809 with H added) as 2 x 2 x 1 unit cells (4 cells, one TOT layer per
cell), one tetrahedral Si -> Al, then the two octahedral Al -> Mg chosen by the static ClayFF energy (every pair of the 16
octahedral sites, with 3 Na+ in the interlayer, each relaxed with the box free; 4 random Na placements per pair, the lowest
energy kept).  Pairs that would make an O bonded to both the tetrahedral Al and a Mg (ClayFF obss, not in the soapff
library) are reported but not chosen.  Dry (no water); d001 starts at 10 A.
    python build_structure.py [nproc]  -> dry.extxyz (elements, positions, cell only), screening.csv
"""
import csv
import itertools
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import numpy as np
from ase import Atoms
from ase.io import write

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../src"))
from soapff_cui.core import read_lammps_data, write_system, TYPE_ALIAS  # noqa: E402
from soapff_cui._clayff.ff import ff_types, CLAYFF_Q  # noqa: E402

LMP = os.environ.get("LMP", "/Users/harutokono/ScoreMD/.pixi/envs/default/bin/lmp")
REF = "/Users/harutokono/clay/soapff-cui/examples/minerals/refs_clayff/pyrophyllite_bulk.data"
D001 = 10.0

ref, _, names = read_lammps_data(REF)
f = ref.get_scaled_positions()
keep = (f[:, 2] < 0.25) | (f[:, 2] >= 0.75)              # the TOT layer centred at z = 0 (octahedral sheet), 2 x 2 x 1
cell = np.array(ref.cell)
cell[2] *= 0.5
lay = Atoms(ref.get_chemical_symbols(), positions=ref.positions, cell=ref.cell, pbc=True)[keep]
fz = lay.get_scaled_positions(wrap=False)[:, 2]
fz = np.where(fz >= 0.75, fz - 1.0, fz)                  # layer centred on the octahedral plane
pos = lay.positions + np.outer(fz - lay.get_scaled_positions(wrap=False)[:, 2], ref.cell[2])
c = cell[2] / np.linalg.norm(cell[2])
cell[2] = c * D001 / (c @ np.cross(cell[0], cell[1]) / np.linalg.norm(np.cross(cell[0], cell[1])))   # d001 = 10 A
pos += 0.5 * cell[2]                                      # octahedral plane in the middle of the cell
layer = Atoms(lay.get_chemical_symbols(), positions=pos, cell=cell, pbc=True)
el = np.array(layer.get_chemical_symbols())
assert len(layer) == 160 and (el == "Al").sum() == 16 and (el == "Si").sum() == 32, len(layer)
si = np.nonzero(el == "Si")[0]
al_t = si[np.argmax(layer.positions[si, 2])]               # one Si of the upper tetrahedral sheet -> Al
octa = np.nonzero(el == "Al")[0]


def structure(mg, na_frac):
    sym = el.copy()
    sym[al_t] = "Al"
    sym[list(mg)] = "Mg"
    zmid = (0.5 + 0.5) % 1.0                              # interlayer midplane: frac z = 0 (= 1)
    na = [np.array([fx, fy, zmid]) @ cell for fx, fy in na_frac]
    return Atoms(list(sym) + ["Na"] * 3, positions=np.vstack([layer.positions, na]), cell=cell, pbc=True)


def energy(args, keep=None):
    """relaxed ClayFF energy of one (Mg pair, Na placement); keep: a path to write the relaxed structure to."""
    mg, k = args
    rng = np.random.default_rng(1000 + 17 * k + 131 * mg[0] + mg[1])
    na_frac = rng.random((3, 2))
    at = structure(mg, na_frac)
    t = ff_types(at)
    if "obss" in t:
        return mg, k, None, True
    q = np.array([CLAYFF_Q[x] for x in t])
    with tempfile.TemporaryDirectory() as d:
        write_system(at, q, t, d)
        open(f"{d}/in.min", "w").write("""include system.settings
neighbor 2.0 bin
thermo_style custom step pe press
min_style cg
min_modify dmax 0.05
minimize 1.0e-6 1.0e-6 2000 20000
fix br all box/relax tri 0.0 vmax 0.001
minimize 1.0e-8 1.0e-8 2000 20000
unfix br
minimize 1.0e-8 1.0e-8 2000 20000
print "PE $(pe:%.6f) LZ $(lz:%.4f)"
write_data relaxed.data
""")
        r = subprocess.run([LMP, "-in", "in.min", "-log", "none", "-echo", "none"], cwd=d, capture_output=True, text=True)
        line = next((l for l in r.stdout.splitlines() if l.startswith("PE ")), None)
        if line is None:
            return mg, k, None, False
        w = line.split()
        if keep:
            a_, _, _ = read_lammps_data(f"{d}/relaxed.data")
            write(keep, a_, format="extxyz")
        return mg, k, (float(w[1]), float(w[3]), na_frac.tolist()), False


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    jobs = [(mg, k) for mg in itertools.combinations(octa.tolist(), 2) for k in range(4)]
    with Pool(nproc) as p:
        res = p.map(energy, jobs)
    best = {}
    obss = set()
    for mg, k, r, is_obss in res:
        if is_obss:
            obss.add(mg)
        elif r and (mg not in best or r[0] < best[mg][0]):
            best[mg] = r
    rows = sorted(best.items(), key=lambda kv: kv[1][0])
    e0 = rows[0][1][0]
    with open(os.path.join(HERE, "screening.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mg1", "mg2", "mg_mg_distance_A", "pe_kcal_mol", "pe_minus_best", "d001_A"])
        for mg, (pe, lz, _) in rows:
            dmm = layer.get_distance(mg[0], mg[1], mic=True)
            w.writerow([mg[0] + 1, mg[1] + 1, f"{dmm:.3f}", f"{pe:.4f}", f"{pe - e0:.4f}", f"{lz:.4f}"])
    (mg, (pe, lz, na_frac)) = rows[0]
    kbest = next(k for m_, k, r, _ in res if m_ == mg and r and r[0] == pe)
    energy((mg, kbest), keep=os.path.join(HERE, "dry.extxyz"))      # the relaxed structure (box and atoms) of the best
    print(f"pairs {len(best)} (+ {len(obss)} with obss, skipped); best Mg {mg[0] + 1} {mg[1] + 1} "
          f"(Mg-Mg {layer.get_distance(mg[0], mg[1], mic=True):.2f} A), next {rows[1][1][0] - e0:+.3f}, "
          f"worst {rows[-1][1][0] - e0:+.3f} kcal/mol; relaxed d001 {lz:.3f} A")
