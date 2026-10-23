"""Water placed with packmol in the empty regions of a structure (vacuum above a surface, an interlayer).

The structure is turned to the standard orientation (a along x, b in the xy plane), so that the layers and surfaces of a
cell with c out of the a-b plane lie along z.  gaps(): the empty slices along z between the atoms (periodic), with their
area and the number of water molecules at a density.  region(): any slice of z, optionally limited in a and b (fractions of
the cell), e.g. only the upper part of a vacuum or a pore of a slab.  add_water(): packmol fills the chosen regions with rigid
water molecules (the geometry of the water model: O-H 1.0 A, H-O-H 109.47 deg) and, if asked, with ions of ClayFF Table 1
(Na, K, Cs, Ca, Ba, Cl; ION_CHARGE), keeping `tolerance` from every atom and from each other.  packmol
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
ION_CHARGE = {"Na": 1, "K": 1, "Cs": 1, "Ca": 2, "Ba": 2, "Cl": -1}        # the ions of ClayFF Table 1 (formal charges)


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


def region(at, lo, hi, xr=(0.0, 1.0), yr=(0.0, 1.0), density=1.0):
    """a box to fill: z from lo to hi (A, standard orientation), the fractions xr of a and yr of b; the same dict as the
    items of gaps() (area of the part of the in-plane cell, thickness hi - lo, n_water at `density`)"""
    s = standard(at)
    area = abs(s.cell[0, 0] * s.cell[1, 1]) * (xr[1] - xr[0]) * (yr[1] - yr[0])
    if not (hi > lo and 0.0 <= xr[0] < xr[1] <= 1.0 and 0.0 <= yr[0] < yr[1] <= 1.0):
        raise ValueError("a region needs lo < hi and 0 <= from < to <= 1 for the a and b fractions")
    if hi - lo >= s.cell[2, 2]:
        raise ValueError("the region is thicker than the cell along z")
    return dict(lo=float(lo), hi=float(hi), thickness=float(hi - lo), area=float(area), xr=tuple(xr), yr=tuple(yr),
                n_water=n_water(area * (hi - lo), density))


def net_charge(ions):
    """formal charge of [(element, n), ...] (ClayFF Table 1 ions)"""
    return sum(ION_CHARGE[e] * n for e, n in ions)


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


def add_water(at, chosen, tolerance=2.0, seed=1, packmol=None, log=print, ions=()):
    """at + water (+ ions): chosen = [(region dict from gaps() or region(), number of molecules), ...]; ions = [(element, n), ...]
    with the elements of ION_CHARGE, shared among the regions in proportion to their molecules (an equal share when no
    region has water).  A region may carry xr, yr (fractions of a and b; default the whole cell).  Returns the new Atoms
    (standard orientation, the water and then the ions appended in each region) and per region the numbers placed and
    the smallest distance to the other atoms."""
    ions = [(e, int(n)) for e, n in ions if int(n) > 0]
    for e, _ in ions:
        if e not in ION_CHARGE:
            raise ValueError(f"ion {e}: not an ion of ClayFF Table 1 ({', '.join(ION_CHARGE)})")
    chosen = [(g, int(n)) for g, n in chosen if int(n) > 0 or ions]
    if not chosen:
        return standard(at), []
    exe = find_packmol(packmol)
    s = standard(at)
    ax, bx, by = s.cell[0, 0], s.cell[1, 0], s.cell[1, 1]
    c = np.array(s.cell[2])
    shear = bx / by
    h = np.radians(HOH / 2)
    template = np.array([[0, 0, 0], [OH * np.sin(h), OH * np.cos(h), 0], [-OH * np.sin(h), OH * np.cos(h), 0]])
    # the share of each ion in each region (the last region takes the remainder)
    weights = np.array([max(n, 0) for _, n in chosen], float)
    weights = weights / weights.sum() if weights.sum() > 0 else np.full(len(chosen), 1.0 / len(chosen))
    share = []
    for e, n in ions:
        k = [int(round(n * w)) for w in weights]
        k[-1] += n - sum(k)
        if k[-1] < 0:
            k[-1] = 0
            k[int(np.argmax(k))] += n - sum(k)
        share.append(k)
    report = []
    for ir, (gap, n) in enumerate(chosen):
        here_ions = [(e, share[k][ir]) for k, (e, _) in enumerate(ions) if share[k][ir] > 0]
        if n <= 0 and not here_ions:
            continue
        # move the region to the middle of the cell along c (a rigid translation), so that it does not wrap
        mid = 0.5 * (gap["lo"] + gap["hi"])
        f = 0.5 - mid / c[2]
        s.positions += f * c
        s.wrap()
        lo, hi = gap["lo"] + f * c[2], gap["hi"] + f * c[2]
        xr, yr = gap.get("xr", (0.0, 1.0)), gap.get("yr", (0.0, 1.0))
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
            box = (f"inside box {xr[0] * ax:.4f} {yr[0] * by:.4f} {lo + m:.4f} {xr[1] * ax:.6f} {yr[1] * by:.6f} {hi - m:.4f}")
            inp = (f"tolerance {tolerance}\nfiletype xyz\noutput packed.xyz\nseed {seed}\n"
                   f"pbc 0.0 0.0 {zmin:.4f} {ax:.6f} {by:.6f} {zmax:.4f}\n"
                   f"structure solid.xyz\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n")
            if n > 0:
                inp += f"structure water.xyz\n  number {n}\n  {box}\nend structure\n"
            for e, k in here_ions:
                _xyz(os.path.join(tmp, f"ion_{e}.xyz"), [e], [[0.0, 0.0, 0.0]])
                inp += f"structure ion_{e}.xyz\n  number {k}\n  {box}\nend structure\n"
            open(os.path.join(tmp, "pack.inp"), "w").write(inp)
            r = subprocess.run([exe], stdin=open(os.path.join(tmp, "pack.inp")), cwd=tmp, capture_output=True, text=True)
            out = os.path.join(tmp, "packed.xyz")
            if not os.path.isfile(out):
                out = os.path.join(tmp, "packed.xyz_FORCED")
            if not os.path.isfile(out):
                raise RuntimeError("packmol wrote no structure:\n" + r.stdout[-2000:])
            converged = "Success!" in r.stdout
            sy, X = _read_xyz(out)
            X = X[len(fixed):]
            W = X[:3 * n].reshape(-1, 3, 3)
            I = X[3 * n:]
            Isym = [e for e, k in here_ions for _ in range(k)]
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        # back to the cell: each molecule moved rigidly by the shear of its O; each ion by its own shear
        W[:, :, 0] += shear * W[:, :1, 1]
        I[:, 0] += shear * I[:, 1]
        water = Atoms(["O", "H", "H"] * len(W), positions=W.reshape(-1, 3), cell=s.cell, pbc=True)
        ion = Atoms(Isym, positions=I, cell=s.cell, pbc=True)
        new = water + ion
        dmin = float(get_distances(new.positions, s.positions, cell=s.cell, pbc=True)[1].min()) if len(new) else 9.9
        s += new
        s.positions -= f * c
        s.wrap()
        report.append(dict(lo=gap["lo"], hi=gap["hi"], n=len(W), ions={e: k for e, k in here_ions}, converged=converged,
                           min_distance=round(dmin, 2)))
        log(f"water: {len(W)} molecules" + "".join(f", {k} {e}" for e, k in here_ions)
            + f" in z {gap['lo']:.2f}-{gap['hi']:.2f} A (packmol {'converged' if converged else 'not converged'}"
            f", nearest other atom {dmin:.2f} A)")
    return s, report
