#!/usr/bin/env python3
"""LAMMPS data file (gouge.data, 01_equil.data, 03_consolidate_<s>.data, ...) -> extended XYZ for OVITO (no numpy needed).
  python3 analysis/data2xyz.py RUN_N54_w0_nro/gouge.data            -> RUN_N54_w0_nro/gouge.xyz
  python3 analysis/data2xyz.py RUN_N54_w0_nro/gouge.data --wrap     the atoms into the box (default: as in the file; gouge.data keeps each
                                                                    particle whole, so a particle can stick out of the box)
  python3 analysis/data2xyz.py RUN_N54_w0_nro/gouge.data --clay     only the clay (no water)
Columns: element, x y z, mol (one number per particle; water: its own), lmp_type, particle (1, 2, ... = clay particle; 0 = water).
OVITO reads the box from Lattice=.  In OVITO: Color coding by "particle" colours each clay particle; Expression selection "particle == 0"
+ Delete selected hides the water.
"""
import argparse, math, os, re, sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('data')
ap.add_argument('-o', '--out', default=None, help='default: the data file with .xyz')
ap.add_argument('--wrap', action='store_true', help='wrap the atoms into the box')
ap.add_argument('--clay', action='store_true', help='leave out the water')
a = ap.parse_args()

ELEM = (('H', 1.008), ('Li', 6.94), ('C', 12.011), ('N', 14.007), ('O', 15.999), ('F', 18.998), ('Na', 22.990), ('Mg', 24.305),
        ('Al', 26.982), ('Si', 28.086), ('K', 39.098), ('Ca', 40.078), ('Fe', 55.845), ('Cl', 35.45), ('Cs', 132.91))
def element(m):
    return min(ELEM, key=lambda e: abs(e[1] - m))[0]

lo = [0.0] * 3; hi = [0.0] * 3; tilt = [0.0] * 3
masses = {}; atoms = []; sec = None; nwater_types = set()
with open(a.data) as f:
    for line in f:
        s = line.split('#')[0].strip()
        m = re.match(r'^([-\d.eE+]+)\s+([-\d.eE+]+)\s+([xyz])lo\s+[xyz]hi', s)
        if m:
            k = 'xyz'.index(m.group(3)); lo[k], hi[k] = float(m.group(1)), float(m.group(2)); continue
        m = re.match(r'^([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+xy\s+xz\s+yz', s)
        if m: tilt = [float(m.group(i)) for i in (1, 2, 3)]; continue
        head = line.strip().split('#')[0].strip()
        if head in ('Masses', 'Atoms', 'Velocities', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs'):
            sec = head; continue
        if not s: continue
        w = s.split()
        if sec == 'Masses' and len(w) >= 2:
            masses[int(w[0])] = float(w[1])
            c = line.split('#')[1].strip().lower() if '#' in line else ''
            if c in ('hw', 'ow'): nwater_types.add(int(w[0]))
        elif sec == 'Atoms' and len(w) >= 7:
            atoms.append((int(w[0]), int(w[1]), int(w[2]), float(w[4]), float(w[5]), float(w[6]),
                          tuple(int(x) for x in w[7:10]) if len(w) >= 10 else (0, 0, 0)))
if not atoms: sys.exit(f'{a.data}: no Atoms section')
if not nwater_types:                       # data files written by LAMMPS lose the comments: water = the two highest types
    nwater_types = set(sorted(masses)[-2:])
L = [hi[k] - lo[k] for k in range(3)]
xy, xz, yz = tilt
atoms.sort()
# the clay particles: their mol numbers (the water molecules have the higher ones)
pmol = sorted({at[1] for at in atoms if at[2] not in nwater_types})
pidx = {m: i + 1 for i, m in enumerate(pmol)}
out = a.out or os.path.splitext(a.data)[0] + '.xyz'
rows = []
for (i, mol, t, x, y, z, img) in atoms:
    if a.clay and t in nwater_types: continue
    # image flags of a data file written by LAMMPS: unwrap with them; --wrap: into the box
    x += img[0] * L[0] + img[1] * xy + img[2] * xz; y += img[1] * L[1] + img[2] * yz; z += img[2] * L[2]
    if a.wrap:
        k = math.floor((z - lo[2]) / L[2]); z -= k * L[2]; y -= k * yz; x -= k * xz
        k = math.floor((y - lo[1]) / L[1]); y -= k * L[1]; x -= k * xy
        k = math.floor((x - lo[0]) / L[0]); x -= k * L[0]
    rows.append(f'{element(masses[t])} {x:.4f} {y:.4f} {z:.4f} {mol} {t} {pidx.get(mol, 0)}\n')
with open(out, 'w') as f:
    f.write(f'{len(rows)}\n')
    f.write(f'Lattice="{L[0]:.5f} 0 0 {xy:.5f} {L[1]:.5f} 0 {xz:.5f} {yz:.5f} {L[2]:.5f}" Origin="{lo[0]:.5f} {lo[1]:.5f} {lo[2]:.5f}" '
            f'Properties=species:S:1:pos:R:3:mol:I:1:lmp_type:I:1:particle:I:1 pbc="T T T"\n')
    f.writelines(rows)
print(f'{out}: {len(rows)} atoms ({len(pmol)} clay particles), box {L[0]:.1f} x {L[1]:.1f} x {L[2]:.1f} A'
      + (' (wrapped)' if a.wrap else ' (particles whole)'))
