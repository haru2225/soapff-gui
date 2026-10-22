"""Surface (slab in vacuum) references for the SOAP library: every low-index face of a bulk reference.

A face (h k l) is given in the unit cell of the bulk reference (the reference itself may be a supercell of it).  The slab:
  1. orientation: the crystal in a cell with two lattice vectors in the plane and a lattice vector out of it (oriented);
  2. termination: a window of a whole number of periods (>= min_thickness) along the normal; every framework cation in the
     window keeps all its O (no undercoordinated Si or tetrahedral Al), small fragments are dropped.  Of all window
     positions, the one with the smallest net formal charge after the healing, then the fewest OH removed to neutralise
     it, then the fewest cut cation-O bonds, then the middle of such a range (a cleavage in the middle of a gap) is taken;
  3. healing (White & Zelazny at near-neutral pH, as in the edge procedure of Lammers et al. used by Underwood & Bourg),
     by the bond valence BVc = sum Z/CN of the cations an O keeps:  BVc > 1: as it is (e.g. Si-O-Al(oct),
     Al(tet)-O-Al(oct)); 0.5 < BVc <= 1: one H in all (Si-OH, Al(tet)-OH, Al2-OH, Mg2-OH), an added H points towards the
     lost cations and is turned away from the other atoms when needed; BVc <= 0.5 (an O left on a single octahedral Al
     or Mg, whatever H it had): a metal with one such O has an OH2 there, a metal with several has one OH (e.g. Al-OH)
     and OH2 on the others (=Al(OH)(OH2) of White & Zelazny); every OH2 is removed and its charge given to the metal
     (Underwood & Bourg): <metal>_e<n>.  When the slab is still negative, OH of such metals become OH2 (removed) as well,
     one per missing charge, alternating between the faces (a protonation towards the point of zero charge);
  4. types: an O that lost cations is oe_<codes> (clayff.py), its charge by the bond-valence rule; other atoms keep the
     type of the bulk reference;
  5. ions (K, Na, Cs ...): a basal cleavage (no bond cut) gets half of the ions of the cleaved interlayer on each face;
     otherwise the ions inside the window.  Water of the reference is not used (the library leaves water out of the SOAP).
The slab is periodic in the plane with `vacuum` A of vacuum along the normal (more than twice the SOAP cut-off, so that
the two faces do not see each other).  Its net formal charge is reported: a slab that these rules cannot neutralise stays
charged (its protonation, i.e. pH, is not changed further)."""
import itertools
import math

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list

from .clayff import BV_ZCN, o_type_name, register_types

# cation-O bond lengths (A, upper limits) by ClayFF type, and the formal charges of the ions
BOND_MAX = {"st": 2.0, "at": 2.05, "ao": 2.3, "mgo": 2.45, "mgh": 2.45, "cao": 2.75, "cah": 2.75, "feo": 2.35, "lio": 2.5}
ION_Z = {"Na": 1, "K": 1, "Cs": 1, "Li": 1, "Rb": 1, "Ca": 2, "Ba": 2, "Sr": 2, "Cl": -1}
OH_LEN = 1.0
H_CLEAR = 1.8                     # an added H is turned until it is at least this far from every other atom (if possible)


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
    """cation-O and H-O bonds of a periodic structure: list of (i_cation_or_H, j_O, vector i->j)"""
    t = np.array(types)
    el = np.array(at.get_chemical_symbols())
    I, J, D, V = neighbor_list("ijdD", at, max(BOND_MAX.values()))
    out = []
    for i, j, d, v in zip(I, J, D, V):
        if el[j] != "O":
            continue
        if t[i] in BOND_MAX and d < BOND_MAX[t[i]]:
            out.append((i, j, v))
        elif el[i] == "H" and t[i] == "ho" and d < 1.25:
            out.append((i, j, v))
    return out


def _place_h(slab, o_index, u0):
    """position of an H on the O slab[o_index]: along u0 (towards the lost cations) or, when that is crowded, the direction
    on cones around u0 (up to 80 deg) that is farthest from the other atoms"""
    o = slab.positions[o_index]
    u0 = u0 / np.linalg.norm(u0)
    a = np.cross(u0, [1.0, 0, 0] if abs(u0[0]) < 0.9 else [0, 1.0, 0]); a /= np.linalg.norm(a)
    b = np.cross(u0, a)
    others = [k for k in range(len(slab)) if k != o_index]
    best, best_d = None, -1.0
    for th in (0, 20, 40, 60, 80):
        for ph in (range(0, 360, 30) if th else (0,)):
            t, p = np.radians(th), np.radians(ph)
            u = np.cos(t) * u0 + np.sin(t) * (np.cos(p) * a + np.sin(p) * b)
            x = o + OH_LEN * u
            probe = slab.copy(); probe.append("H"); probe.positions[-1] = x
            dmin = probe.get_distances(len(probe) - 1, others, mic=True).min() if others else 9.9
            if dmin > best_d:
                best, best_d = x, dmin
        if best_d >= H_CLEAR:
            break
    return best, best_d


def make_slab(bulk, types, miller_unit, min_thickness=16.0, vacuum=16.0, step=0.1, seed=0):
    """slab of the bulk reference (Atoms, ClayFF type per atom; no water) for the face miller_unit of its unit cell.
    Returns (Atoms, types, info)."""
    types = list(types)
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
    big = ori.repeat((1, 1, nper + 3))
    tb = [types[k] for k in big.get_array("tid")]
    el = np.array(big.get_chemical_symbols())
    z = big.positions @ nrm
    z -= z.min()
    is_cat = np.array([x in BOND_MAX for x in tb])
    is_ion = np.array([x in ION_Z for x in el]) & ~is_cat
    zcat = np.array([BV_ZCN[x][0] if x in BV_ZCN and x != "ho" else 0 for x in tb])
    bvc = np.array([BV_ZCN[x][0] / BV_ZCN[x][1] if x in BV_ZCN and x != "ho" else 0.0 for x in tb])
    cat_of_o, h_of_o = {}, {}
    for i, j, v in _bonds(big, tb):
        (cat_of_o if is_cat[i] else h_of_o).setdefault(j, []).append((i, v))
    o_list = np.array(sorted(cat_of_o))
    o_cats = [np.array([i for i, _ in cat_of_o[o]]) for o in o_list]
    o_nh = np.array([len(h_of_o.get(o, [])) for o in o_list])
    T = nper * d

    def classify(inside):
        """per O touched by the window: kept as it is, "terminal" (BV of the kept cations <= 0.5: on a single octahedral
        metal) or healed to its number of H; returns (cut bonds, formal charge without the terminal O, metal -> terminal O)"""
        cut, formal, terminal = 0, int(zcat[inside].sum()), {}
        for o, cs, nh in zip(o_list, o_cats, o_nh):
            k = inside[cs]
            if not k.any():
                continue
            if k.all():
                formal += nh - 2
                continue
            cut += int((~k).sum())
            bvcat = bvc[cs[k]].sum()
            if bvcat <= 0.5 + 1e-9 and k.sum() == 1:
                terminal.setdefault(int(cs[k][0]), []).append(o)
                continue
            formal += max(nh, 1 if bvcat <= 1.0 + 1e-9 else 0) - 2
        return cut, formal, terminal

    def evaluate(z0):
        """(net formal charge left after the neutralisation, OH removed for it, cut bonds) of the window at z0"""
        inside = is_cat & (z >= z0) & (z < z0 + T)
        cut, formal, terminal = classify(inside)
        n_oh = sum(1 for os_ in terminal.values() if len(os_) >= 2)   # one OH (formal -1) per metal with several terminal O
        formal += -n_oh + sum(ION_Z[el[i]] for i in np.nonzero(is_ion & (z >= z0) & (z < z0 + T))[0])
        used = min(n_oh, max(0, -formal))                              # OH turned into OH2 (removed): +1 each
        return abs(formal + used), used, cut
    z0s = np.arange(d, 2 * d, step)
    score = [evaluate(z0) for z0 in z0s]
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
    cut_best = score[run[len(run) // 2] % len(z0s)][2]
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
        for h, _ in h_of_o.get(o, []):
            keep[h] = True
    # ions
    basal = cut_best == 0
    ion_idx = np.nonzero(is_ion)[0]
    moved = {}
    if basal:
        zf = z[~is_ion]                       # the framework of the whole crystal: the interlayer gap around z0
        gap = (zf[zf < z0].max(), zf[zf >= z0].min())
        band = [i for i in ion_idx if gap[0] < z[i] < gap[1]]
        keep[[i for i in ion_idx if gap[1] <= z[i] <= gap[0] + T]] = True
        rng = np.random.default_rng(seed)
        stay = set(np.array(band)[rng.permutation(len(band))[: len(band) // 2]].tolist()) if band else set()
        for i in band:
            keep[i] = True
            if i not in stay:
                moved[i] = nper * a3                                   # the same interlayer site on the other face
    else:
        keep[[i for i in ion_idx if z0 <= z[i] < z0 + T]] = True
    # healing
    new_t = {i: tb[i] for i in np.nonzero(keep)[0]}
    add_h, remove, metal_n = [], set(), {}
    n_oe = 0
    zmid = z[keep].mean()
    _, formal0, terminal = classify(inside)
    formal0 += sum(ION_Z[el[i]] for i in np.nonzero(keep & is_ion)[0])
    # one OH (the most exposed O) on a metal with several terminal O; when the slab is negative, as many of these OH as
    # needed become OH2 (removed) too, alternating between the two faces
    as_oh = {max(os_, key=lambda o: abs(z[o] - zmid)) for os_ in terminal.values() if len(os_) >= 2}
    formal0 -= len(as_oh)
    neutralised = 0
    if formal0 < 0 and as_oh:
        top = sorted((o for o in as_oh if z[o] >= zmid), key=lambda o: -z[o])
        bot = sorted((o for o in as_oh if z[o] < zmid), key=lambda o: z[o])
        order = [o for pair in itertools.zip_longest(top, bot) for o in pair if o is not None]
        for o in order[:-formal0]:
            as_oh.discard(o); neutralised += 1
    term_o = {o for os_ in terminal.values() for o in os_}
    for o, cs in zip(o_list, o_cats):
        if not keep[o]:
            continue
        kept = [c for c in cs if inside[c]]
        lost = [(c, v) for c, v in cat_of_o[o] if not inside[c]]
        if not lost:
            continue
        hs = [h for h, _ in h_of_o.get(o, []) if keep[h]]
        if o in term_o and o not in as_oh:                         # OH2: removed, its charge goes to the metal
            remove.add(o); remove.update(hs)
            for c in kept:
                metal_n[c] = metal_n.get(c, 0) + 1
            continue
        bvcat = sum(bvc[c] for c in kept)
        n_add = max(0, (1 if (o in as_oh or bvcat <= 1.0 + 1e-9) else 0) - len(hs))
        if n_add:
            u = -sum(v / np.linalg.norm(v) for _, v in lost)               # v: cation -> O; the H points to the lost cations
            add_h.append((o, u if np.linalg.norm(u) > 1e-6 else -lost[0][1]))
        new_t[o] = o_type_name([tb[c] for c in kept] + ["ho"] * (len(hs) + n_add))
        n_oe += 1
    for c, n in metal_n.items():
        new_t[c] = f"{tb[c]}_e{n}"
    idx = [i for i in sorted(new_t) if i not in remove]
    pos = np.array([big.positions[i] + moved.get(i, 0.0) for i in idx])
    zc = pos @ nrm
    shift = nrm * (zc.min() - vacuum / 2)
    height = (zc.max() - zc.min()) + vacuum
    slab = Atoms([el[i] for i in idx], positions=pos - shift, cell=np.array([a1, a2, nrm * height]), pbc=True)
    typ = [new_t[i] for i in idx]
    where = {i: k for k, i in enumerate(idx)}
    h_min = 9.9
    for o, u in add_h:
        x, dmin = _place_h(slab, where[o], u)
        slab.append("H"); slab.positions[-1] = x; typ.append("ho")
        h_min = min(h_min, dmin)
    slab.wrap()
    register_types(typ)
    sym = slab.get_chemical_symbols()
    I, J, D = neighbor_list("ijd", slab, 1.5)
    bonded_oh = {(min(i, j), max(i, j)) for i, j, dd in zip(I, J, D) if {sym[i], sym[j]} == {"O", "H"} and dd < 1.1}
    close = [dd for i, j, dd in zip(I, J, D) if i < j and (i, j) not in bonded_oh]
    formal = (sum(BV_ZCN[x.split("_e")[0]][0] for x in typ if x.split("_e")[0] in BV_ZCN)
              + sum(ION_Z[s] for s, x in zip(sym, typ) if x == s and s in ION_Z) - 2 * sym.count("O"))
    info = dict(miller=tuple(miller_unit), miller_cell=m, repeats=reps, period=round(d, 3), periods=nper, area=round(area, 2),
                cut_bonds=int(cut_best), basal=bool(basal), atoms=len(slab), oe=n_oe, h_added=len(add_h),
                oh2_removed=len([o for o in remove if el[o] == "O"]), oh_neutralised=neutralised, fragments_dropped=dropped,
                ions=int(sum(1 for s, x in zip(sym, typ) if x == s and s in ION_Z)), formal_charge=int(round(formal)),
                close_contacts=len(close), min_contact=round(min(close, default=9.9), 2), h_clearance=round(h_min, 2))
    return slab, typ, info


def hkl_name(m):
    return "".join(f"m{-x}" if x < 0 else str(x) for x in m)


def write_surfaces(data, out_dir, max_index=1, min_thickness=16.0, vacuum=16.0, log=print):
    """slabs of every face |h|,|k|,|l| <= max_index of the bulk reference `data` (LAMMPS data file with type names; its
    water is left out) -> out_dir/<name>_<hkl>.data (ClayFF charges; references for build-library) and .xyz; returns
    {tag: info}"""
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
        slab, typ, info = make_slab(at, types, m, min_thickness=min_thickness, vacuum=vacuum)
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
        log(f"{tag}: {info['atoms']} atoms, {info['cut_bonds']} cation-O bonds cut{' (basal cleavage)' if info['basal'] else ''}, "
            f"{info['h_added']} H added, {info['oh2_removed']} OH2 removed, net formal charge {info['formal_charge']:+d}")
    return summary
