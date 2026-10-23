"""Water and ions placed with packmol: in the empty regions of a structure (vacuum above a surface, an interlayer), or in any
region built from the geometric constraints of packmol (inside/outside box, cube, sphere, ellipsoid, cylinder; over/below
plane, xygauss), as in PACKMOL-GUI (VMD plugin), restricted to the molecules ClayFF types: the water (rigid, the geometry of
the water model: O-H 1.0 A, H-O-H 109.47 deg) and the ions of ClayFF Table 1 (Na, K, Cs, Ca, Ba, Cl; ION_CHARGE).

The structure is turned to the standard orientation (a along x, b in the xy plane), so that the layers and surfaces of a
cell with c out of the a-b plane lie along z.  packmol works in a space where the in-plane cell is rectangular
(x' = x - (b_x/b_y) y = x - y cot(gamma), exact for the lattice), with its periodic box in x' and y: the coordinates of the
constraints are given in that frame (x' = x when gamma = 90 deg), z as in the standard orientation (it may run past the top
of the cell: the regions are moved together along c, rigidly, so that none crosses the cell boundary).  Each molecule is
moved back rigidly.
  gaps(): the empty slices along z between the atoms (periodic), with their area and the number of water molecules at a
density.  region(): any slice of z, optionally limited in a and b (fractions of the cell).  estimate(): the volume of a
region of constraints (Monte Carlo), the part of it outside the van der Waals spheres of the atoms, and the molecules that
fill it at a density.  pack(): packmol places items (molecule, number, constraints), keeping `tolerance` from every atom
and from each other; packmol_input() is the input it runs.  add_water(): the regions of gaps()/region() filled with water
and ions (each region an inside-box constraint).  The water and ions are typed by rule in soapff_cui.core.make."""
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
MOLECULES = ("water",) + tuple(ION_CHARGE)                              # what can be placed (ClayFF types them)
# packmol constraints: the parameters of each shape (packmol's order) and the locations it takes
SHAPES = {
    "box": (("xmin", "ymin", "zmin", "xmax", "ymax", "zmax"), ("inside", "outside")),
    "cube": (("xmin", "ymin", "zmin", "side"), ("inside", "outside")),
    "sphere": (("x", "y", "z", "radius"), ("inside", "outside")),
    "ellipsoid": (("x", "y", "z", "a", "b", "c", "d"), ("inside", "outside")),      # ((x-x0)/a)^2 + ... <= d
    "cylinder": (("x", "y", "z", "dx", "dy", "dz", "radius", "length"), ("inside", "outside")),   # from (x y z) along (dx dy dz)
    "plane": (("a", "b", "c", "d"), ("over", "below")),                             # a x + b y + c z >= d (over)
    "xygauss": (("x0", "y0", "sx", "sy", "z0", "height"), ("over", "below")),       # z >= z0 + height exp(-...) (over)
}


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




def _frame(at):
    """the standard-oriented copy and the packmol frame: (s, ax, by, shear, c, co); x' = x - shear y, co = c in that frame"""
    s = standard(at)
    ax, bx, by = s.cell[0, 0], s.cell[1, 0], s.cell[1, 1]
    shear = bx / by
    c = np.array(s.cell[2])
    co = c.copy()
    co[0] -= shear * c[1]
    return s, ax, by, shear, c, co


def frame_of(at):
    """the structure in the frame of the constraints: dict(positions (x' in [0, ax), y in [0, by), z), symbols, ax, by,
    shear (x' = x - shear y), cz (the height of the cell), z_atoms (lowest and highest z))"""
    s, ax, by, shear, c, co = _frame(at)
    P = s.positions.copy()
    P[:, 0] -= shear * P[:, 1]
    P[:, 0] %= ax
    return dict(positions=P, symbols=s.get_chemical_symbols(), ax=float(ax), by=float(by), shear=float(shear),
                cz=float(c[2]), z_atoms=(float(P[:, 2].min()), float(P[:, 2].max())))


def check_constraint(con):
    """(location, shape, parameters) -> the same with float parameters; ValueError if packmol would not take it"""
    loc, shape, par = con
    if shape not in SHAPES:
        raise ValueError(f"shape {shape}: one of {', '.join(SHAPES)}")
    names, locs = SHAPES[shape]
    if loc not in locs:
        raise ValueError(f"a {shape} takes {' or '.join(locs)}, not {loc}")
    par = [float(x) for x in par]
    if len(par) != len(names):
        raise ValueError(f"a {shape} needs {len(names)} numbers ({' '.join(names)})")
    if shape == "box" and not all(par[k] < par[k + 3] for k in range(3)):
        raise ValueError("box: xmin < xmax, ymin < ymax, zmin < zmax")
    if shape in ("cube", "sphere") and par[3] <= 0:
        raise ValueError(f"{shape}: the {names[3]} must be positive")
    if shape == "ellipsoid" and min(par[3:]) <= 0:
        raise ValueError("ellipsoid: a, b, c and d must be positive")
    if shape == "cylinder" and (np.linalg.norm(par[3:6]) == 0 or par[6] <= 0 or par[7] <= 0):
        raise ValueError("cylinder: a direction, a positive radius and a positive length")
    if shape == "plane" and np.linalg.norm(par[:3]) == 0:
        raise ValueError("plane: a, b, c (the normal) cannot all be 0")
    if shape == "xygauss" and (par[2] <= 0 or par[3] <= 0):
        raise ValueError("xygauss: sx and sy must be positive")
    return (loc, shape, par)


def _zspan(con):
    """z range (frame) of a constraint that bounds the region (inside a box, cube, sphere, ellipsoid, cylinder); else None"""
    loc, shape, p = con
    if loc != "inside":
        return None
    if shape == "box":
        return p[2], p[5]
    if shape == "cube":
        return p[2], p[2] + p[3]
    if shape == "sphere":
        return p[2] - p[3], p[2] + p[3]
    if shape == "ellipsoid":
        r = p[5] * np.sqrt(p[6])
        return p[2] - r, p[2] + r
    if shape == "cylinder":
        d = np.array(p[3:6]) / np.linalg.norm(p[3:6])
        e = p[2] + p[7] * d[2]
        rz = p[6] * np.sqrt(max(0.0, 1.0 - d[2] ** 2))
        return min(p[2], e) - rz, max(p[2], e) + rz
    return None


def _shift(con, t):
    """the constraint moved by the vector t (frame)"""
    loc, shape, p = con
    p = list(p)
    if shape == "box":
        for k in range(3):
            p[k] += t[k]
            p[k + 3] += t[k]
    elif shape in ("cube", "sphere", "ellipsoid", "cylinder"):
        for k in range(3):
            p[k] += t[k]
    elif shape == "plane":
        p[3] += p[0] * t[0] + p[1] * t[1] + p[2] * t[2]
    elif shape == "xygauss":
        p[0] += t[0]
        p[1] += t[1]
        p[4] += t[2]
    return (loc, shape, p)


def _test(con, X):
    """which points X (n x 3, frame) satisfy the constraint (packmol's definitions)"""
    loc, shape, p = con
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    if shape == "box":
        ok = (x >= p[0]) & (x <= p[3]) & (y >= p[1]) & (y <= p[4]) & (z >= p[2]) & (z <= p[5])
    elif shape == "cube":
        ok = (x >= p[0]) & (x <= p[0] + p[3]) & (y >= p[1]) & (y <= p[1] + p[3]) & (z >= p[2]) & (z <= p[2] + p[3])
    elif shape == "sphere":
        ok = ((X - p[:3]) ** 2).sum(1) <= p[3] ** 2
    elif shape == "ellipsoid":
        ok = ((x - p[0]) / p[3]) ** 2 + ((y - p[1]) / p[4]) ** 2 + ((z - p[2]) / p[5]) ** 2 <= p[6]
    elif shape == "cylinder":
        d = np.array(p[3:6]) / np.linalg.norm(p[3:6])
        v = X - np.array(p[:3])
        a = v @ d
        ok = (a >= 0) & (a <= p[7]) & ((v ** 2).sum(1) - a ** 2 <= p[6] ** 2)
    elif shape == "plane":
        ok = p[0] * x + p[1] * y + p[2] * z - p[3] >= 0
    else:                                                    # xygauss
        ok = z >= p[4] + p[5] * np.exp(-((x - p[0]) ** 2 / (2 * p[2] ** 2) + (y - p[1]) ** 2 / (2 * p[3] ** 2)))
    return ~ok if loc in ("outside", "below") else ok


def _item_span(item):
    zs = [z for z in (_zspan(c) for c in item["constraints"]) if z]
    if not zs:
        raise ValueError(f"{item['molecule']}: give an inside constraint (box, cube, sphere, ellipsoid or cylinder) that bounds "
                         "the region")
    lo, hi = max(z[0] for z in zs), min(z[1] for z in zs)
    if hi <= lo:
        raise ValueError(f"{item['molecule']}: the inside constraints do not overlap along z")
    return lo, hi


def _arrange(items, cz):
    """per item the number of cells k (along c) its constraints are moved down by, and the z range of all the regions
    then; ValueError when they do not fit in one cell along z"""
    spans = [_item_span(it) for it in items]
    best = None
    for z0, _ in spans:
        ks = [int(np.floor((lo - z0) / cz + 1e-9)) for lo, _ in spans]
        los = [lo - k * cz for (lo, _), k in zip(spans, ks)]
        his = [hi - k * cz for (_, hi), k in zip(spans, ks)]
        span = max(his) - min(los)
        if best is None or span < best[0]:
            best = (span, ks, min(los), max(his))
    span, ks, zlo, zhi = best
    if span >= cz:
        raise ValueError(f"the regions span {span:.1f} A along z, more than the cell ({cz:.1f} A)")
    return ks, zlo, zhi


def _vdw(el):
    from ase.data import atomic_numbers, vdw_radii
    r = vdw_radii[atomic_numbers[el]]
    return 1.5 if not np.isfinite(r) else float(r)


def estimate(at, constraints, density=1.0, n_points=200000, seed=0):
    """Monte Carlo over one region (the constraints of an item, frame coordinates): its volume, the part of it outside the
    van der Waals spheres of the atoms (with their periodic images), the water molecules that fill that part at `density`
    g/cm3, and up to 3000 points of it (frame coordinates) for drawing"""
    from scipy.spatial import cKDTree
    cons = [check_constraint(c) for c in constraints]
    s, ax, by, shear, c, co = _frame(at)
    lo, hi = _item_span(dict(molecule="region", constraints=cons))
    rng = np.random.default_rng(seed)
    X = rng.random((n_points, 3)) * [ax, by, hi - lo] + [0.0, 0.0, lo]
    ok = np.ones(len(X), bool)
    for con in cons:
        ok &= _test(con, X)
    vol = float(ok.mean() * ax * by * (hi - lo))
    P = X[ok].copy()
    P[:, 0] += shear * P[:, 1]                               # frame -> Cartesian
    free = np.ones(len(P), bool)
    if len(P):
        cz = c[2]
        ks = range(int(np.floor(lo / cz)) - 1, int(np.ceil(hi / cz)) + 2)
        shifts = np.array([i * s.cell[0] + j * s.cell[1] + k * c for i in (-1, 0, 1) for j in (-1, 0, 1) for k in ks])
        sym = np.array(s.get_chemical_symbols())
        for e in sorted(set(sym)):
            A = (s.positions[sym == e][None, :, :] + shifts[:, None, :]).reshape(-1, 3)
            d, _ = cKDTree(A).query(P, distance_upper_bound=_vdw(e))
            free &= ~np.isfinite(d)
    vfree = vol * float(free.mean()) if len(P) else 0.0
    pts = X[ok][free]
    if len(pts) > 3000:
        pts = pts[rng.choice(len(pts), 3000, replace=False)]
    return dict(volume=vol, free_volume=vfree, n_water=n_water(vfree, density), z=(lo, hi), points=pts)


def _molecule_file(m):
    return "water.xyz" if m == "water" else f"ion_{m}.xyz"


def _prepare(at, items, tolerance=2.0, seed=1):
    items = [dict(molecule=it["molecule"], number=int(it["number"]), constraints=[check_constraint(c) for c in it["constraints"]])
             for it in items]
    for it in items:
        if it["molecule"] not in MOLECULES:
            raise ValueError(f"{it['molecule']}: only the water and the ions of ClayFF Table 1 ({', '.join(ION_CHARGE)}) can be placed")
        if it["number"] < 0:
            raise ValueError(f"{it['molecule']}: a negative number")
    items = [it for it in items if it["number"] > 0]
    if not items:
        raise ValueError("nothing to place (every number is 0)")
    s, ax, by, shear, c, co = _frame(at)
    ks, zlo, zhi = _arrange(items, c[2])
    f = (0.5 * c[2] - 0.5 * (zlo + zhi)) / c[2]              # the regions to the middle of the cell along c
    s.positions += f * c
    s.wrap()
    P = s.positions.copy()
    P[:, 0] -= shear * P[:, 1]
    fixed = np.vstack([P, P + co, P - co])                   # the structure with its images along c, in the packmol box
    fixed[:, 0] %= ax
    fixed[:, 1] %= by
    zmin, zmax = fixed[:, 2].min() - 10.0, fixed[:, 2].max() + 10.0
    L = [f"# soapff-cui (soapff_cui.water.pack): the structure ({len(s)} atoms) fixed, with its images along c",
         f"# frame: x' = x - y cot(gamma) (shear {shear:.6f}); structure and regions moved by {f:+.4f} c along c "
         "(the coordinates below are moved)",
         f"tolerance {tolerance:g}", "filetype xyz", "output packed.xyz", f"seed {int(seed)}",
         f"pbc 0.0 0.0 {zmin:.4f} {ax:.6f} {by:.6f} {zmax:.4f}", "",
         "structure solid.xyz", "  number 1", "  fixed 0. 0. 0. 0. 0. 0.", "end structure"]
    for it, k in zip(items, ks):
        t = (f - k) * co
        L += ["", f"structure {_molecule_file(it['molecule'])}", f"  number {it['number']}"]
        for con in it["constraints"]:
            loc, shape, p = _shift(con, t)
            L.append(f"  {loc} {shape} " + " ".join(f"{x:.4f}" for x in p))
        L.append("end structure")
    return dict(s=s, f=f, c=c, shear=shear, fixed=fixed, sym=s.get_chemical_symbols() * 3, items=items, input="\n".join(L) + "\n")


def packmol_input(at, items, tolerance=2.0, seed=1):
    """the packmol input pack() runs for these items (text)"""
    return _prepare(at, items, tolerance, seed)["input"]


def pack(at, items, tolerance=2.0, seed=1, packmol=None, workdir=None, log=print):
    """at + molecules placed by packmol.  items: [dict(molecule="water" or an ion of ION_CHARGE, number=N,
    constraints=[(location, shape, parameters), ...]), ...] (frame coordinates, see the module).  Every item needs an inside
    constraint that bounds it (box, cube, sphere, ellipsoid, cylinder).  workdir: keep packmol's files there (pack.inp, the
    xyz files, packed.xyz, packmol.log).  Returns the new Atoms (standard orientation, the molecules appended in the order of
    the items) and a report: converged, per item the number placed and the nearest atom of the structure, the input."""
    exe = find_packmol(packmol)
    pr = _prepare(at, items, tolerance, seed)
    s, f, c, shear, fixed = pr["s"], pr["f"], pr["c"], pr["shear"], pr["fixed"]
    h = np.radians(HOH / 2)
    template = np.array([[0, 0, 0], [OH * np.sin(h), OH * np.cos(h), 0], [-OH * np.sin(h), OH * np.cos(h), 0]])
    tmp = workdir or tempfile.mkdtemp(prefix="soapff_packmol_")
    os.makedirs(tmp, exist_ok=True)
    try:
        _xyz(os.path.join(tmp, "solid.xyz"), pr["sym"], fixed)
        for m in {it["molecule"] for it in pr["items"]}:
            if m == "water":
                _xyz(os.path.join(tmp, "water.xyz"), ["O", "H", "H"], template)
            else:
                _xyz(os.path.join(tmp, _molecule_file(m)), [m], [[0.0, 0.0, 0.0]])
        open(os.path.join(tmp, "pack.inp"), "w").write(pr["input"])
        for old in ("packed.xyz", "packed.xyz_FORCED"):
            if os.path.isfile(os.path.join(tmp, old)):
                os.remove(os.path.join(tmp, old))
        r = subprocess.run([exe], stdin=open(os.path.join(tmp, "pack.inp")), cwd=tmp, capture_output=True, text=True)
        if workdir:
            open(os.path.join(tmp, "packmol.log"), "w").write(r.stdout + r.stderr)
        out = os.path.join(tmp, "packed.xyz")
        if not os.path.isfile(out):
            out = os.path.join(tmp, "packed.xyz_FORCED")
        if not os.path.isfile(out):
            raise RuntimeError("packmol wrote no structure:\n" + (r.stdout + r.stderr)[-2000:])
        _, X = _read_xyz(out)
    finally:
        if not workdir:
            shutil.rmtree(tmp, ignore_errors=True)
    converged = "Success!" in r.stdout
    X = X[len(fixed):]
    n0 = len(s)
    rep, k = [], 0
    for it in pr["items"]:
        na = 3 if it["molecule"] == "water" else 1
        M = X[k:k + na * it["number"]].reshape(-1, na, 3)
        k += na * it["number"]
        M[:, :, 0] += shear * M[:, :1, 1]                    # back to Cartesian: each molecule rigidly, by the shear of its first atom
        sym = ["O", "H", "H"] if na == 3 else [it["molecule"]]
        new = Atoms(sym * len(M), positions=M.reshape(-1, 3), cell=s.cell, pbc=True)
        dmin = float(get_distances(new.positions, s.positions[:n0], cell=s.cell, pbc=True)[1].min()) if len(new) else 9.9
        rep.append(dict(molecule=it["molecule"], n=len(M), min_distance=round(dmin, 2)))
        s += new
    s.positions -= f * c
    s.wrap()
    log(f"packmol {'converged' if converged else 'did NOT converge (packed.xyz_FORCED: overlaps may remain)'}: "
        + ", ".join(f"{r_['n']} {r_['molecule']} (nearest atom {r_['min_distance']:.2f} A)" for r_ in rep))
    return s, dict(converged=converged, items=rep, input=pr["input"], log=r.stdout[-4000:])


def add_water(at, chosen, tolerance=2.0, seed=1, packmol=None, log=print, ions=()):
    """at + water (+ ions): chosen = [(region dict from gaps() or region(), number of molecules), ...]; ions = [(element, n), ...]
    with the elements of ION_CHARGE, shared among the regions in proportion to their molecules (an equal share when no
    region has water).  A region may carry xr, yr (fractions of a and b; default the whole cell); it is an inside-box
    constraint of pack() (0.5 A in from its z ends).  All regions are packed in one packmol run (each region by itself when
    they do not fit in one cell along z).  Returns the new Atoms (standard orientation) and per region the numbers placed,
    whether packmol converged and the smallest distance to the other atoms."""
    ions = [(e, int(n)) for e, n in ions if int(n) > 0]
    for e, _ in ions:
        if e not in ION_CHARGE:
            raise ValueError(f"ion {e}: not an ion of ClayFF Table 1 ({', '.join(ION_CHARGE)})")
    chosen = [(g, int(n)) for g, n in chosen if int(n) > 0 or ions]
    if not chosen:
        return standard(at), []
    s, ax, by, shear, c, co = _frame(at)
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
    groups = []                                              # per region: its items
    for ir, (gap, n) in enumerate(chosen):
        lo, hi = gap["lo"], gap["hi"]
        xr, yr = gap.get("xr", (0.0, 1.0)), gap.get("yr", (0.0, 1.0))
        m = min(0.5, 0.25 * (hi - lo))
        box = ("inside", "box", [xr[0] * ax, yr[0] * by, lo + m, xr[1] * ax, yr[1] * by, hi - m])
        its = ([dict(molecule="water", number=n, constraints=[box])] if n > 0 else []) + \
              [dict(molecule=e, number=share[k][ir], constraints=[box]) for k, (e, _) in enumerate(ions) if share[k][ir] > 0]
        groups.append((gap, its))
    try:
        _arrange([it for _, its in groups for it in its], c[2])
        runs = [[g for g in groups if g[1]]]
    except ValueError:
        runs = [[g] for g in groups if g[1]]                 # the regions do not fit in one cell along z: one run each
    report, cur = [], at
    for run in runs:
        new, rep = pack(cur, [it for _, its in run for it in its], tolerance=tolerance, seed=seed, packmol=packmol,
                        log=lambda *a: None)
        k = 0
        for gap, its in run:
            r_ = rep["items"][k:k + len(its)]
            k += len(its)
            nw = sum(x["n"] for x in r_ if x["molecule"] == "water")
            ion_n = {x["molecule"]: x["n"] for x in r_ if x["molecule"] != "water"}
            dmin = min(x["min_distance"] for x in r_)
            report.append(dict(lo=gap["lo"], hi=gap["hi"], n=nw, ions=ion_n, converged=rep["converged"], min_distance=dmin))
            log(f"water: {nw} molecules" + "".join(f", {v} {e}" for e, v in ion_n.items())
                + f" in z {gap['lo']:.2f}-{gap['hi']:.2f} A (packmol {'converged' if rep['converged'] else 'not converged'}"
                f", nearest other atom {dmin:.2f} A)")
        cur = new
    return cur, report
