"""Water placed with packmol in the empty regions of a structure (vacuum above a surface, an interlayer).

The structure is turned to the standard orientation (a along x, b in the xy plane), so that the layers and surfaces of a
cell with c out of the a-b plane lie along z.  gaps(): the empty slices along z between the atoms (periodic), with their
area and the number of water molecules at a density.  add_water(): packmol fills the chosen gaps with rigid water
molecules (the geometry of the water model: O-H 1.0 A, H-O-H 109.47 deg), keeping `tolerance` from every atom.  packmol
works in a space where the in-plane cell is rectangular (x' = x - (b_x/b_y) y, exact for the lattice), with its periodic
box in x and y; each molecule is moved back rigidly.  The water is typed by rule in soapff_cui.core.make."""
import os
import shutil
import subprocess
import tempfile

import numpy as np
from ase import Atoms
from ase.geometry import cellpar_to_cell, get_distances

AVOGADRO = 6.02214076e23
M_WATER = 18.01528
OH, HOH = 1.0, 109.47


def standard(at):
    """the same structure with the cell in the standard orientation (a along x, b in the xy plane)"""
    s = at.copy()
    s.set_cell(cellpar_to_cell(at.cell.cellpar()), scale_atoms=True)
    s.pbc = True
    s.wrap()
    return s


def gaps(at, min_gap=2.5, density=1.0):
    """empty slices along z of the standard-oriented structure: list of dict(lo, hi, thickness, area, n_water) (A, A^2;
    n_water at `density` g/cm3 in area x thickness, the distance between the atom planes that bound it), thickest first"""
    s = standard(at)
    cz = s.cell[2, 2]
    z = np.sort(s.get_scaled_positions()[:, 2] * cz)
    area = abs(s.cell[0, 0] * s.cell[1, 1])
    zz = np.append(z, z[0] + cz)
    out = []
    for k in range(len(z)):
        t = zz[k + 1] - zz[k]
        if t >= min_gap:
            out.append(dict(lo=float(zz[k]), hi=float(zz[k + 1]), thickness=float(t), area=float(area),
                            n_water=n_water(area * t, density)))
    return sorted(out, key=lambda g: -g["thickness"])


def n_water(volume_A3, density=1.0):
    """number of water molecules in volume_A3 at density g/cm3"""
    return int(round(density * volume_A3 * 1e-24 * AVOGADRO / M_WATER))


def find_packmol(path=None):
    for p in (path, os.environ.get("PACKMOL"), shutil.which("packmol"),
              os.path.join(os.path.expanduser("~"), ".soapff-gui", "packmol", "bin", "packmol"),
              os.path.join(os.path.expanduser("~"), ".soapff-gui", "lammps", "bin", "packmol")):
        if p and os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    raise FileNotFoundError("packmol not found: install it (install.sh, conda-forge 'packmol', or brew install packmol) "
                            "or set PACKMOL")


def _xyz(path, sym, pos):
    with open(path, "w") as f:
        f.write(f"{len(sym)}\nsoapff\n")
        for e, p in zip(sym, pos):
            f.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


def _read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    rows = [l.split() for l in L[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def add_water(at, chosen, tolerance=2.0, seed=1, packmol=None, log=print):
    """at + water: chosen = [(gap dict from gaps(), number of molecules), ...].  Returns the new Atoms (standard
    orientation, the water appended) and per gap the number of molecules placed and the smallest water-other distance."""
    exe = find_packmol(packmol)
    s = standard(at)
    ax, bx, by = s.cell[0, 0], s.cell[1, 0], s.cell[1, 1]
    c = np.array(s.cell[2])
    shear = bx / by
    h = np.radians(HOH / 2)
    template = np.array([[0, 0, 0], [OH * np.sin(h), OH * np.cos(h), 0], [-OH * np.sin(h), OH * np.cos(h), 0]])
    report = []
    for gap, n in chosen:
        if n <= 0:
            continue
        # move the gap to the middle of the cell along c (a rigid translation), so that it does not wrap
        mid = 0.5 * (gap["lo"] + gap["hi"])
        f = 0.5 - mid / c[2]
        s.positions += f * c
        s.wrap()
        lo, hi = gap["lo"] + f * c[2], gap["hi"] + f * c[2]
        # the rectangular space of packmol: x' = x - shear * y; the atoms, with their images along +-c
        P = s.positions.copy()
        P[:, 0] -= shear * P[:, 1]
        P[:, 0] %= ax
        P[:, 1] %= by
        co = c.copy(); co[0] -= shear * c[1]
        sym = s.get_chemical_symbols()
        fixed = np.vstack([P, P + co, P - co])
        tmp = tempfile.mkdtemp(prefix="soapff_packmol_")
        try:
            _xyz(os.path.join(tmp, "solid.xyz"), sym * 3, fixed)
            _xyz(os.path.join(tmp, "water.xyz"), ["O", "H", "H"], template)
            zmin, zmax = fixed[:, 2].min() - 10, fixed[:, 2].max() + 10
            m = min(0.5, 0.25 * (hi - lo))
            inp = (f"tolerance {tolerance}\nfiletype xyz\noutput packed.xyz\nseed {seed}\n"
                   f"pbc 0.0 0.0 {zmin:.4f} {ax:.6f} {by:.6f} {zmax:.4f}\n"
                   f"structure solid.xyz\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
                   f"structure water.xyz\n  number {n}\n  inside box 0.0 0.0 {lo + m:.4f} {ax:.6f} {by:.6f} {hi - m:.4f}\nend structure\n")
            open(os.path.join(tmp, "pack.inp"), "w").write(inp)
            r = subprocess.run([exe], stdin=open(os.path.join(tmp, "pack.inp")), cwd=tmp, capture_output=True, text=True)
            out = os.path.join(tmp, "packed.xyz")
            if not os.path.isfile(out):
                out = os.path.join(tmp, "packed.xyz_FORCED")
            if not os.path.isfile(out):
                raise RuntimeError("packmol wrote no structure:\n" + r.stdout[-2000:])
            converged = "Success!" in r.stdout
            _, W = _read_xyz(out)
            W = W[len(fixed):].reshape(-1, 3, 3)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        # back to the cell: each molecule moved rigidly by the shear of its O
        W[:, :, 0] += shear * W[:, :1, 1]
        water = Atoms(["O", "H", "H"] * len(W), positions=W.reshape(-1, 3), cell=s.cell, pbc=True)
        dmin = float(get_distances(water.positions, s.positions, cell=s.cell, pbc=True)[1].min()) if len(water) else 9.9
        s += water
        s.positions -= f * c
        s.wrap()
        report.append(dict(lo=gap["lo"], hi=gap["hi"], n=len(W), converged=converged, min_distance=round(dmin, 2)))
        log(f"water: {len(W)} molecules in z {gap['lo']:.2f}-{gap['hi']:.2f} A (packmol {'converged' if converged else 'not converged'}"
            f", nearest other atom {dmin:.2f} A)")
    return s, report
