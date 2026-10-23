"""Surface (slab in vacuum) references for the SOAP library: the basal plane of a layered mineral.

A face (h k l) is given in the unit cell of the bulk reference (the reference itself may be a supercell of it).  Only a face
that can be cut without cutting any cation-O bond is accepted: the basal plane between the layers of a layered mineral.  A
face that cuts bonds wherever it is placed (an edge, any face of quartz) is refused; nothing is added to the atoms.  The slab:
  1. orientation: the crystal in a cell with two lattice vectors in the plane and a lattice vector out of it (oriented);
  2. termination: a window of a whole number of periods (>= min_thickness) along the normal; every framework cation in the
     window keeps all its O and their H, small fragments are dropped.  Of the window positions that cut no cation-O bond,
     the one with the smallest net formal charge, then the middle of such a range (a cleavage in the middle of the gap), is
     taken;
  3. types: every atom keeps the ClayFF type of the bulk reference (its charge: ClayFF Table 1, clayff.type_charges);
  4. ions (K, Na, Cs ...): half of the ions of the cleaved interlayer on each face.  Water of the reference is not used
     (the library leaves water out of the SOAP).
The slab is periodic in the plane with `vacuum` A of vacuum along the normal (or cut inside the cell: make_slab).  Its net
formal charge is reported."""
import itertools
import math

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list


# cation-O bond lengths (A, upper limits) by ClayFF type, and the formal charges of the ions
BOND_MAX = {"st": 2.0, "at": 2.05, "ao": 2.3, "mgo": 2.45, "mgh": 2.45, "cao": 2.75, "cah": 2.75, "feo": 2.35, "lio": 2.5}
ION_Z = {"Na": 1, "K": 1, "Cs": 1, "Li": 1, "Rb": 1, "Ca": 2, "Ba": 2, "Sr": 2, "Cl": -1}
CATION_Z = {"st": 4, "at": 3, "ao": 3, "mgo": 2, "mgh": 2, "cao": 2, "cah": 2, "feo": 3, "lio": 1}   # formal valences
MAX_INDEX = 1                     # only the faces |h|, |k|, |l| <= 1 of the unit cell (those of the library) can be cut
GAP_MIN = 6.0                     # the empty part of a cell cut inside the cell: wider than the SOAP cut-off (5 A)
AUTO_THICKNESS = 10.0             # the slab of the automatic cut inside the cell is at least this thick (A)


def unit_repeats(at, types, tol=0.25):
    """repeats (n1, n2, n3) of the unit cell in the supercell `at`: the largest n along each axis for which a translation
    by a_i / n maps the O sublattice onto itself (substitutions are ignored)"""
    pos = at.positions[np.array(at.get_chemical_symbols()) == "O"]
    cell = np.array(at.cell)
    inv = np.linalg.inv(cell)
    reps = []
    for i in range(3):
        best = 1
        for n in range(8, 1, -1):
            f = ((pos + cell[i] / n)[:, None, :] - pos[None, :, :]) @ inv
            if (np.linalg.norm((f - np.round(f)) @ cell, axis=2).min(1) < tol).all():
                best = n
                break
        reps.append(best)
    return tuple(reps)


def faces(max_index=1):
    """the Miller indices (h k l) with |h|, |k|, |l| <= max_index, one of each +/- pair, no common divisor"""
    out = []
    for h, k, l in itertools.product(range(-max_index, max_index + 1), repeat=3):
        if (h, k, l) == (0, 0, 0) or math.gcd(math.gcd(abs(h), abs(k)), abs(l)) != 1:
            continue
        if next(x for x in (h, k, l) if x) > 0:
            out.append((h, k, l))
    return sorted(out, key=lambda m: (sum(map(abs, m)), [-abs(x) for x in m], m))


def oriented(at, miller, tol=1e-10):
    """the same crystal in a cell whose first two vectors span the plane (h k l) (cell coordinates of `at`) and whose third
    is a lattice vector out of the plane (unimodular change of basis, as in ase.build.surface, but without replacing the
    third vector by its projection on the normal, which is not a lattice vector)"""
    from ase.build.general_surface import ext_gcd
    h, k, l = miller
    h0, k0, l0 = h == 0, k == 0, l == 0
    if (h0 and k0) or (h0 and l0) or (k0 and l0):
        if not h0:
            C = [(0, 1, 0), (0, 0, 1), (1, 0, 0)]
        elif not k0:
            C = [(0, 0, 1), (1, 0, 0), (0, 1, 0)]
        else:
            C = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    else:
        p, q = ext_gcd(k, l)
        a1, a2, a3 = np.array(at.cell)
        k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3), l * a2 - k * a3)
        k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3), l * a2 - k * a3)
        if abs(k2) > tol:
            i = -int(round(k1 / k2))
            p, q = p + i * l, q - i * k
        a, b = ext_gcd(p * k + q * l, h)
        C = [(p * k + q * l, -p * h, -q * h), tuple(np.array((0, l, -k)) // abs(math.gcd(l, k))), (b, a * p, a * q)]
    C = np.array(C)
    if round(abs(np.linalg.det(C))) != 1:
        raise ValueError(f"basis for {miller} is not unimodular")
    if np.linalg.det(C @ np.array(at.cell)) < 0:
        C[2] = -C[2]
    o = at.copy()
    o.set_cell(C @ np.array(at.cell), scale_atoms=False)
    o.pbc = True
    o.wrap()
    return o


def _bonds(at, types):
    """cation-O and H-O bonds of a periodic structure: list of (i_cation_or_H, j_O)"""
    t = np.array(types)
    el = np.array(at.get_chemical_symbols())
    I, J, D = neighbor_list("ijd", at, max(BOND_MAX.values()))
    out = []
    for i, j, d in zip(I, J, D):
        if el[j] != "O":
            continue
        if t[i] in BOND_MAX and d < BOND_MAX[t[i]]:
            out.append((i, j))
        elif el[i] == "H" and t[i] == "ho" and d < 1.25:
            out.append((i, j))
    return out


def make_slab(bulk, types, miller_unit, min_thickness=16.0, vacuum=16.0, step=0.1, seed=0, keep_cell=False, repeat=(1, 1, 1),
              periods=None):
    """slab of the bulk reference (Atoms, ClayFF type per atom; no water) for the face miller_unit of its unit cell.
    keep_cell: the vacuum is added along the lattice vector that leaves the plane (for a face of the cell, e.g. (100): a is
    lengthened, b and c stay as they are), so the cell keeps the shape and the directions of the input cell; the vacuum
    along the normal is still `vacuum`.  Otherwise the third vector is the normal (a rectangular stacking direction).
    repeat: (na, nb, nc), the bulk cell is repeated before the cut (the size of the cell).
    periods: (an int, or "auto" = all unit cells but one) cut inside the cell, without changing it: the slab is `periods` unit cells thick (unit cells along the direction out
    of the plane, as counted in the cell) and the rest of the cell is empty; `vacuum` is then only the extra length added to
    the cell (0: the cell is exactly the input cell, in the same directions).  The empty part must be at least 2 A wide.
    Only a face that can be cut without cutting any cation-O bond (the basal plane between the layers of a layered mineral) is
    accepted, at a position that cuts none; a face that cuts bonds everywhere (an edge, every face of quartz) is refused.
    Returns (Atoms, types, info)."""
    types = list(types)
    if max(abs(int(x)) for x in miller_unit) > MAX_INDEX:
        raise ValueError(f"face {tuple(miller_unit)}: only faces with |h|, |k|, |l| <= {MAX_INDEX} of the unit cell can be cut (the faces "
                         "of the library; faces with larger indices are not validated)")
    if tuple(repeat) != (1, 1, 1):
        bulk = bulk.repeat(tuple(int(x) for x in repeat))
        types = types * int(np.prod(repeat))
    reps = unit_repeats(bulk, types)
    m = [x * n for x, n in zip(miller_unit, reps)]
    g = math.gcd(math.gcd(abs(m[0]), abs(m[1])), abs(m[2]))
    m = tuple(x // g for x in m)
    b = bulk.copy()
    b.set_array("tid", np.arange(len(b)))
    ori = oriented(b, m)
    a1, a2, a3 = np.array(ori.cell)
    nrm = np.cross(a1, a2); area = np.linalg.norm(nrm); nrm /= area
    d = abs(np.dot(a3, nrm))
    nper = max(1, math.ceil(min_thickness / d))
    n_u = 1                                              # unit cells along the vector out of the plane
    cell_n = 1                                           # periods of the output cell along it (more than 1: the cell is lengthened)
    u = d                                                # the length of a unit cell along it
    if periods is not None:
        for i in range(3):
            if np.allclose(a3, np.array(b.cell[i])) or np.allclose(a3, -np.array(b.cell[i])):
                n_u = reps[i]
        u = d / n_u
        if periods == "auto":                            # all unit cells but one stay, and at least AUTO_THICKNESS of slab
            periods = max(n_u - 1, math.ceil(AUTO_THICKNESS / u - 1e-9), 1)
        if periods < 1:
            raise ValueError("keep at least 1 unit cell")
        nper = max(1, math.ceil(periods * u / d - 1e-9))
        cell_n = max(1, math.ceil((periods * u + GAP_MIN) / d - 1e-9))   # the empty part of the cell is at least GAP_MIN wide
    big = ori.repeat((1, 1, nper + 3))
    tb = [types[k] for k in big.get_array("tid")]
    el = np.array(big.get_chemical_symbols())
    z = big.positions @ nrm
    z -= z.min()
    is_cat = np.array([x in BOND_MAX for x in tb])
    is_ion = np.array([x in ION_Z for x in el]) & ~is_cat
    zcat = np.array([CATION_Z.get(x, 0) for x in tb])
    cat_of_o, h_of_o = {}, {}
    for i, j in _bonds(big, tb):
        (cat_of_o if is_cat[i] else h_of_o).setdefault(j, []).append(i)
    o_list = np.array(sorted(cat_of_o))
    o_cats = [np.array(cat_of_o[o]) for o in o_list]
    o_nh = np.array([len(h_of_o.get(o, [])) for o in o_list])
    T = nper * d if periods is None else periods * u

    def evaluate(z0):
        """(cation-O bonds cut, |net formal charge|) of the window at z0: its cations with all their O and H, and its ions"""
        inside = is_cat & (z >= z0) & (z < z0 + T)
        cut, formal = 0, int(zcat[inside].sum())
        for cs, nh in zip(o_cats, o_nh):
            k = inside[cs]
            if k.any():
                cut += int((~k).sum())
                formal += int(nh) - 2
        formal += sum(ION_Z[el[i]] for i in np.nonzero(is_ion & (z >= z0) & (z < z0 + T))[0])
        return cut, abs(formal)
    z0s = np.arange(d, 2 * d, step)
    score = [evaluate(z0) for z0 in z0s]
    if min(s_[0] for s_ in score) > 0:
        raise ValueError(f"face {tuple(miller_unit)} cuts cation-O bonds wherever it is placed: only faces that cut no bond (the "
                         "basal plane between the layers of a layered mineral, e.g. (0 0 1)) are supported; edges are not")
    score = [s_[1] if s_[0] == 0 else math.inf for s_ in score]          # only the positions that cut no bond
    best = min(score)
    runs, cur = [], []
    for k, s_ in enumerate(score):
        if s_ == best:
            cur.append(k)
        elif cur:
            runs.append(cur); cur = []
    if cur:
        runs.append(cur)
    if len(runs) > 1 and score[0] == best and score[-1] == best:          # the run wraps around the period
        runs[0] = runs[-1] + runs[0]; runs.pop()
    run = max(runs, key=len)
    z0 = z0s[run[len(run) // 2] % len(z0s)]
    inside = is_cat & (z >= z0) & (z < z0 + T)
    # drop small fragments (cations linked through shared O)
    cats = np.nonzero(inside)[0]
    parent = {c: c for c in cats}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for cs in o_cats:
        k = [c for c in cs if inside[c]]
        for c in k[1:]:
            parent[find(c)] = find(k[0])
    comp = {}
    for c in cats:
        comp.setdefault(find(c), []).append(c)
    biggest = max(len(v) for v in comp.values())
    dropped = 0
    for v in comp.values():
        if len(v) < max(4, 0.1 * biggest):
            inside[v] = False; dropped += len(v)
    keep = inside.copy()
    o_keep = [o for o, cs in zip(o_list, o_cats) if inside[cs].any()]
    keep[o_keep] = True
    for o in o_keep:
        keep[h_of_o.get(o, [])] = True
    # ions: the cleaved interlayer (the gap of the framework around z0) gives half of its ions to each face
    ion_idx = np.nonzero(is_ion)[0]
    zf = z[~is_ion]
    gap = (zf[zf < z0].max(), zf[zf >= z0].min())
    band = [i for i in ion_idx if gap[0] < z[i] < gap[1]]
    keep[[i for i in ion_idx if gap[1] <= z[i] <= gap[0] + T]] = True
    rng = np.random.default_rng(seed)
    stay = set(np.array(band)[rng.permutation(len(band))[: len(band) // 2]].tolist()) if band else set()
    moved = {}
    for i in band:
        keep[i] = True
        if i not in stay:
            moved[i] = (T / d) * a3                                    # the same interlayer site on the other face
    idx = list(np.nonzero(keep)[0])
    pos = np.array([big.positions[i] + moved.get(i, 0.0) for i in idx])
    zc = pos @ nrm
    if periods is not None:                                   # inside the cell: the slab in the middle, the rest of the cell empty
        extent = zc.max() - zc.min()
        while cell_n * d + vacuum - extent < GAP_MIN:         # the faces must not see each other (SOAP cut-off): lengthen the cell
            cell_n += 1
        height = cell_n * d + vacuum
        shift = nrm * (zc.min() - (height - extent) / 2)
        keep_cell = True
    else:
        shift = nrm * (zc.min() - vacuum / 2)
        height = (zc.max() - zc.min()) + vacuum
    third = nrm * height
    if keep_cell:                                              # along the lattice vector out of the plane, perpendicular height as above
        third = a3 / abs(np.dot(a3, nrm)) * height
    slab = Atoms([el[i] for i in idx], positions=pos - shift, cell=np.array([a1, a2, third]), pbc=True)
    typ = [tb[i] for i in idx]
    slab.wrap()
    sym = slab.get_chemical_symbols()
    I, J, D = neighbor_list("ijd", slab, 1.5)
    bonded_oh = {(min(i, j), max(i, j)) for i, j, dd in zip(I, J, D) if {sym[i], sym[j]} == {"O", "H"} and dd < 1.1}
    close = [dd for i, j, dd in zip(I, J, D) if i < j and (i, j) not in bonded_oh]
    formal = (sum(CATION_Z.get(x, 0) for x in typ) + typ.count("ho")
              + sum(ION_Z[s] for s, x in zip(sym, typ) if x == s and s in ION_Z) - 2 * sym.count("O"))
    info = dict(miller=tuple(miller_unit), miller_cell=m, repeats=reps, period=round(d, 3), periods=nper, area=round(area, 2),
                cut_bonds=0, basal=True, atoms=len(slab), fragments_dropped=dropped,
                ions=int(sum(1 for s, x in zip(sym, typ) if x == s and s in ION_Z)), formal_charge=int(round(formal)),
                close_contacts=len(close), min_contact=round(min(close, default=9.9), 2),
                keep_cell=bool(keep_cell), cell_periods=periods, cell_doubled=(cell_n > 1), cell_length_periods=cell_n)
    return slab, typ, info


def hkl_name(m):
    return "".join(f"m{-x}" if x < 0 else str(x) for x in m)


def write_surfaces(data, out_dir, max_index=MAX_INDEX, min_thickness=16.0, vacuum=16.0, keep_cell=False, repeat=(1, 1, 1), periods=None, log=print):
    """slabs of every face |h|,|k|,|l| <= max_index of the bulk reference `data` (LAMMPS data file with type names; its
    water is left out) that cuts no cation-O bond -> out_dir/<name>_<hkl>.data (ClayFF charges; references for build-library)
    and .xyz; the faces that cut bonds are skipped; returns {tag: info}"""
    import os
    import ase.io
    from .clayff import type_charges
    from . import core
    os.makedirs(out_dir, exist_ok=True)
    at, q, types = core.read_lammps_data(data)
    keep = ~np.isin(types, ["ow", "hw"])
    at, types = at[keep], [t for t, k in zip(types, keep) if k]
    name = os.path.splitext(os.path.basename(data))[0].replace("_bulk", "")
    summary = {}
    for m in faces(max_index):
        tag = f"{name}_{hkl_name(m)}"
        try:
            slab, typ, info = make_slab(at, types, m, min_thickness=min_thickness, vacuum=vacuum, keep_cell=keep_cell, repeat=repeat,
                                          periods=periods)
        except ValueError as e:
            if "cuts cation-O bonds" not in str(e):
                raise
            log(f"{tag}: skipped (cuts cation-O bonds)")
            continue
        tq = type_charges("clayff2004", types=typ)
        qs = np.array([tq[t] for t in typ])
        info["partial_charge"] = round(float(qs.sum()), 4)
        tmp = os.path.join(out_dir, "." + tag)
        os.makedirs(tmp, exist_ok=True)
        core.write_system(slab, qs, typ, tmp, boundary="p p p")
        os.replace(os.path.join(tmp, "system.data"), os.path.join(out_dir, tag + ".data"))
        for f in os.listdir(tmp):
            os.remove(os.path.join(tmp, f))
        os.rmdir(tmp)
        ase.io.write(os.path.join(out_dir, tag + ".xyz"), slab, format="extxyz")
        summary[tag] = info
        log(f"{tag}: {info['atoms']} atoms (basal cleavage, no bond cut), net formal charge {info['formal_charge']:+d}")
    return summary
