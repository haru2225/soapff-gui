"""Library building, charge assignment and the LAMMPS writer of soapff-cui."""
import collections
import csv
import itertools
import json
import os

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list

from ._clayff.chem import MASS
from ._clayff.ff import ELEMENT, ff_types, lammps_frame

from .clayff import (ANGLE_HOH, ANGLE_MOH, BOND_OH, CHARGE_SCHEMES, LJ, WATER_LABELS,
                     WATER_MODELS, WATER_Q,
                     type_charges)
from .clayff import Q as CLAYFF_Q

SOAP_DEFAULT = dict(r_cut=5.0, n_max=4, l_max=3, sigma=0.4)
THRESHOLD = 0.95
FRAMEWORK = ("Si", "Al", "Mg", "Fe", "Ti")
TYPE_ALIAS = {"OW": "ow", "HW": "hw", "OH": "oh", "HO": "ho"}   # library type names -> ClayFF type names (clayff.py)
TIMESTEP_FS = 0.5                                                  # flexible SPC water: 0.5 fs (as in ClayFF studies)
IONS = ("Li", "Na", "K", "Rb", "Cs", "Ca", "Sr", "Ba", "Cl", "Br", "F", "I")
FORMAL_ION = {"Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Ca": 2, "Sr": 2, "Ba": 2, "Cl": -1, "Br": -1, "F": -1, "I": -1}


# ------------------------------------------------------------------------------------------------- reading structures
def read_lammps_data(path):
    """atom_style full data file -> (Atoms, charges, type names or None).  Elements from the masses."""
    L = open(path).read().splitlines()
    box = {}
    tilt = (0.0, 0.0, 0.0)
    for l in L[:60]:
        p = l.split()
        if len(p) >= 4 and p[2] in ("xlo", "ylo", "zlo"):
            box[p[2][0]] = (float(p[0]), float(p[1]))
        if len(p) >= 6 and p[3] == "xy":
            tilt = tuple(map(float, p[:3]))
    i = L.index(next(l for l in L if l.strip().startswith("Masses"))) + 2
    mass, name = {}, {}
    while L[i].strip():
        p = L[i].split()
        mass[int(p[0])] = float(p[1])
        name[int(p[0])] = L[i].split("#")[1].strip() if "#" in L[i] else None
        i += 1
    el_of = {}
    for t, m in mass.items():
        if name[t] in ELEMENT:            # a known type name fixes the element
            el_of[t] = ELEMENT[name[t]]
            continue
        el_of[t] = min(MASS, key=lambda e: abs(MASS[e] - m)) if min(abs(MASS[e] - m) for e in MASS) < 0.1 else _mass_to_el(m)
    i = [k for k, l in enumerate(L) if l.strip().startswith("Atoms")][0] + 2
    rows = []
    while i < len(L) and L[i].strip():
        p = L[i].split()
        rows.append((int(p[0]), int(p[2]), float(p[3]), float(p[4]), float(p[5]), float(p[6])))
        i += 1
    rows.sort()
    lx = box["x"][1] - box["x"][0]; ly = box["y"][1] - box["y"][0]; lz = box["z"][1] - box["z"][0]
    cell = [[lx, 0, 0], [tilt[0], ly, 0], [tilt[1], tilt[2], lz]]
    lo = np.array([box["x"][0], box["y"][0], box["z"][0]])
    at = Atoms([el_of[r[1]] for r in rows], positions=np.array([r[3:6] for r in rows]) - lo, cell=cell, pbc=True)
    names = [name[r[1]] for r in rows] if all(name.values()) else None
    return at, np.array([r[2] for r in rows]), names


def _mass_to_el(m):
    from ase.data import atomic_masses, chemical_symbols
    k = int(np.argmin(np.abs(atomic_masses[1:] - m))) + 1
    return chemical_symbols[k]


def read_structure(path):
    """any structure with a cell -> Atoms (charges/types from a LAMMPS data file are ignored here)."""
    if path.endswith((".data", ".lammps", ".lmp")) or "data" in os.path.basename(path).split(".")[-1]:
        at, _, _ = read_lammps_data(path)
        return at
    import ase.io
    at = ase.io.read(path)
    if not at.cell.any():
        raise ValueError(f"{path}: no cell; give a periodic structure (extxyz with Lattice, cif, POSCAR or LAMMPS data)")
    at.pbc = True
    return at


def read_dump_frames(path, symbols):
    """frames of a LAMMPS dump (id + x/y/z, xu/yu/zu or xs/ys/zs; orthogonal or triclinic box) -> list of Atoms
    (positions wrapped into the cell, atoms sorted by id)."""
    out = []
    with open(path) as f:
        while True:
            l = f.readline()
            if not l:
                return out
            if l.startswith("ITEM: TIMESTEP"):
                f.readline(); f.readline(); n = int(f.readline())
                head = f.readline()
                b = [list(map(float, f.readline().split())) for _ in range(3)]
                if "xy" in head:
                    xy, xz, yz = b[0][2], b[1][2], b[2][2]
                    xlo = b[0][0] - min(0.0, xy, xz, xy + xz); xhi = b[0][1] - max(0.0, xy, xz, xy + xz)
                    ylo = b[1][0] - min(0.0, yz); yhi = b[1][1] - max(0.0, yz)
                else:
                    xy = xz = yz = 0.0
                    xlo, xhi, ylo, yhi = b[0][0], b[0][1], b[1][0], b[1][1]
                cell = np.array([[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, b[2][1] - b[2][0]]])
                cols = f.readline().split()[2:]
                keys = ("xs", "ys", "zs") if "xs" in cols else ("xu", "yu", "zu") if "xu" in cols else ("x", "y", "z")
                # only the id and position columns (a dump may also hold strings, e.g. the type names of soapff-cui)
                a = np.loadtxt([f.readline() for _ in range(n)], usecols=[cols.index("id")] + [cols.index(c) for c in keys], ndmin=2)
                a = a[np.argsort(a[:, 0])]
                xyz = a[:, 1:4] @ cell if keys[0] == "xs" else a[:, 1:4] - np.array([xlo, ylo, b[2][0]])
                at = Atoms(symbols, positions=xyz, cell=cell, pbc=True)
                at.wrap()
                out.append(at)


# ------------------------------------------------------------------------------------------------- roles and checks
def roles(at):
    """water (an O with two H at 0.9-1.1 A making an H-O-H angle of 95-125 deg, the shape of a (flexible) water molecule, and those
    H; also when it overlaps the mineral in a built structure, or sits on a metal as an OH2, which the edge model of
    Underwood & Bourg treats as water), ion (isolated ions), mineral (everything else).  Also returns warnings for atoms
    closer than physically possible."""
    el = np.array(at.get_chemical_symbols())
    n = len(at)
    I, J, D = neighbor_list("ijd", at, 2.4)
    nbr = [[] for _ in range(n)]
    for i, j, d in zip(I, J, D):
        nbr[i].append((j, d))
    role = np.array(["mineral"] * n, dtype="U12")
    warn = []
    for h in np.nonzero(el == "H")[0]:
        os_ = [(j, d) for j, d in nbr[h] if el[j] == "O" and d < 1.2]
        if len(os_) > 1:
            warn.append(f"H{h + 1} is within 1.2 A of {len(os_)} O atoms (" + ", ".join(f"O{j + 1} {d:.2f} A" for j, d in os_)
                        + "): overlapping atoms, minimise the structure first")
    for o in np.nonzero(el == "O")[0]:
        hs = [j for j, d in nbr[o] if el[j] == "H" and 0.9 <= d <= 1.1 and role[j] != "water"]
        if len(hs) < 2:
            continue
        pair = None
        for a, b in itertools.combinations(hs, 2):
            if 95.0 <= at.get_angle(a, o, b, mic=True) <= 125.0:
                pair = (a, b)
                break
        if pair and len(hs) == 2:
            role[o] = "water"; role[list(pair)] = "water"
    for i in np.nonzero(np.isin(el, IONS))[0]:
        role[i] = "ion"
    # overlapping atoms: closer than 0.75 A (an O-H bond 0.6-1.2 A excepted), or two atoms that are not H closer than 1.2 A
    close = [(i, j, d) for i, j, d in zip(I, J, D) if i < j and (
        (d < 0.75 and not ({el[i], el[j]} == {"O", "H"} and d > 0.6)) or (d < 1.2 and "H" not in (el[i], el[j])))]
    for i, j, d in close[:20]:
        warn.append(f"atoms {i + 1} ({el[i]}) and {j + 1} ({el[j]}) are {d:.2f} A apart: overlapping atoms in the structure")
    return role, warn


def fallback_charges(at, role):
    """ClayFF type charges (water flexible SPC, ions formal) - used only for atoms the library cannot describe."""
    at = at.copy()
    at.set_array("role", np.array(["water" if r == "water" else "ion" if r == "ion" else "bulk" for r in role], dtype="U12"))
    types = ff_types(at)
    q = np.array([CLAYFF_Q.get(t, 0.0) for t in types])
    qO, qH = WATER_Q
    el = at.get_chemical_symbols()
    for i, r in enumerate(role):
        if r == "water":
            q[i] = qO if el[i] == "O" else qH
        elif r == "ion":
            q[i] = FORMAL_ION.get(el[i], q[i])
    return q, types


# ------------------------------------------------------------------------------------------------- SOAP library
def soap_vectors(at, species, r_cut, n_max, l_max, sigma, weighting=None, species_weights=None):
    """normalised SOAP vectors (dscribe).  weighting: dscribe radial weighting of the neighbour density, e.g.
    {"function": "pow", "r0": 2.0, "m": 4, "c": 1, "d": 1}: w(r) = c / (d + (r/r0)^m) (Willatt et al. 2018), which
    makes the first coordination shell (the bonded cations) count more than the outer environment."""
    from dscribe.descriptors import SOAP
    kw = dict(weighting=weighting) if weighting else {}
    soap = SOAP(species=list(species), r_cut=r_cut, n_max=n_max, l_max=l_max, sigma=sigma, periodic=True, sparse=False, **kw)
    X = soap.create(at)
    if species_weights:                          # scale the components of each species pair (Z1, Z2) by w(Z1) w(Z2)
        X = np.array(X, dtype=float)
        sp = list(species)
        for i, a in enumerate(sp):
            for b in sp[i:]:
                f = species_weights.get(a, 1.0) * species_weights.get(b, 1.0)
                if f != 1.0:
                    X[:, soap.get_location((a, b))] *= f
    return (X / np.linalg.norm(X, axis=1, keepdims=True)).astype(np.float32)


def _prune_new(x, el, tn, q, kept, threshold):
    """entries of one new structure to keep: those that no kept entry (dict (element, type, charge) -> vectors) of the same
    element, type and charge, and no entry kept before them in this structure, matches with similarity >= threshold"""
    sel = np.zeros(len(x), bool)
    qr = np.round(q, 6)
    for key in set(zip(el, tn, qr)):
        idx = np.nonzero((el == key[0]) & (tn == key[1]) & (qr == key[2]))[0]
        cand = idx
        if key in kept:
            cand = idx[(x[idx] @ kept[key].T).max(1) < threshold]
        chosen = []
        for i in cand:
            if not chosen or (x[chosen] @ x[i]).max() < threshold:
                chosen.append(i)
        sel[chosen] = True
    return sel


def prune_entries(X, q, types, elements, threshold, log=print, protect=0):
    """indices of the entries to keep: an entry is dropped when an already kept entry of the same element, charge and
    type is at least `threshold` similar (it adds no new environment); near entries with another charge or type are
    always kept, so the one-to-one information is not lost.  Greedy, in the order of the library; the first `protect`
    entries are always kept."""
    keep = []
    for e in sorted(set(elements)):
        idx = np.nonzero(elements == e)[0]
        kept = []
        for k0 in range(0, len(idx), 512):
            blk = idx[k0:k0 + 512]
            for i in blk:
                if kept and i >= protect:
                    K = np.array(kept)
                    sim = X[K] @ X[i]
                    same = (np.abs(q[K] - q[i]) < 1e-6) & (types[K] == types[i])
                    if (sim[same] >= threshold).any():
                        continue
                kept.append(i)
        keep += kept
        log(f"  prune {e}: {len(idx)} -> {len(kept)} entries")
    return np.sort(np.array(keep))


def build_library(ref_data, frame_files=(), out="library.npz", every=1, soap=None, prune=None, water_in_soap=False,
                  protect=0, log=print):
    """reference LAMMPS data file(s) (charges and type names in the Masses comments) + optional MD dump frames of any
    of them -> library of (element, SOAP vector, charge, type) for every atom of every structure.  With several reference
    systems the SOAP species are the union of their elements, so that all vectors are comparable.
    frame_files: "N:dump" (frames of the N-th reference, 1-based; same atoms in the same order) or "dump" (= 1:dump).
    Thermal frames add the environments of a structure in motion (an initial structure alone misses some of them).
    water_in_soap=False: the water molecules (types ow/hw) are left out of every reference before the SOAP (no water entries,
    and the water is not part of the environment of the other atoms); make() then types water by rule and leaves it out of the
    SOAP too.  The ions stay (their positions belong to the mineral).
    protect: with prune, the first `protect` references are never pruned (only later ones, e.g. slabs, lose duplicates)."""
    soap = dict(SOAP_DEFAULT, **(soap or {}))
    for k in ("weighting", "species_weights"):
        if not soap.get(k):
            soap.pop(k, None)
    refs = [ref_data] if isinstance(ref_data, str) else list(ref_data)
    loaded = [read_lammps_data(r) for r in refs]
    species = sorted(set(e for at0, _, _ in loaded for e in at0.get_chemical_symbols()))
    structs = []
    masks = [np.ones(len(at0), bool) for at0, _, _ in loaded]
    if not water_in_soap:                                    # water left out of the references
        dry = []
        for k, (at0, q, names) in enumerate(loaded):
            masks[k] = ~np.isin(names, list(WATER_TYPES)) if names else roles(at0)[0] != "water"
            keep = masks[k]
            dry.append((at0[keep], q[keep], [n for n, x in zip(names, keep) if x] if names else None))
        loaded = dry
    for r, (at0, q, names) in zip(refs, loaded):
        names = names or [""] * len(at0)
        structs.append((os.path.basename(r), at0, q, names))
    for spec in frame_files:
        head, sep, tail = spec.partition(":")
        k_ref, fpath = (int(head) - 1, tail) if sep and head.isdigit() else (0, spec)
        if not 0 <= k_ref < len(refs):
            raise ValueError(f"--frames {spec}: there are {len(refs)} reference files")
        at0, q, names = loaded[k_ref][0], loaded[k_ref][1], structs[k_ref][3]
        full = read_lammps_data(refs[k_ref])[0]
        frames = read_dump_frames(fpath, full.get_chemical_symbols())[::every]
        if frames and len(frames[0]) != len(full):
            raise ValueError(f"--frames {spec}: {len(frames[0])} atoms, reference {refs[k_ref]} has {len(full)}")
        frames = [a[masks[k_ref]] for a in frames]
        for k, a in enumerate(frames):
            structs.append((f"{os.path.basename(refs[k_ref])}|{os.path.basename(fpath)}#{k}", a, q, names))
    X, Q, T, E, SRC = [], [], [], [], []
    kept = {}                    # (element, type, charge) -> SOAP vectors kept so far (pruning, structure by structure)
    n_all = 0
    for k, (name, a, q, names) in enumerate(structs):
        x = soap_vectors(a, species, **soap)
        el = np.array(a.get_chemical_symbols()); tn = np.array(names); q = np.asarray(q, float)
        sel = np.ones(len(a), bool)
        n_all += len(a)
        if prune and k >= protect:
            sel = _prune_new(x, el, tn, q, kept, prune)
        for key in set(zip(el[sel], tn[sel], np.round(q[sel], 6))):
            m = sel & (el == key[0]) & (tn == key[1]) & (np.round(q, 6) == key[2])
            kept[key] = np.vstack([kept[key], x[m]]) if key in kept else x[m]
        X.append(x[sel]); Q.append(q[sel]); T += list(tn[sel]); E += list(el[sel]); SRC += [name] * int(sel.sum())
        log(f"  structure {k + 1}/{len(structs)} ({name}): {len(a)} atoms" + (f", {int(sel.sum())} kept" if prune else ""))
    X, Q, T, E, SRC = np.vstack(X), np.concatenate(Q), np.array(T), np.array(E), np.array(SRC)
    if prune:
        log(f"pruned (similarity >= {prune} to a kept entry of the same element, charge and type; the first {protect} "
            f"references kept whole): {n_all} -> {len(Q)} entries")
    np.savez_compressed(out, X=X, q=Q, type=T, element=E, origin=SRC,
                        species=np.array(species), source=np.array([os.path.abspath(r) for r in refs]),
                        n_structures=len(structs), weighting=json.dumps(soap.get("weighting") or {}), water_in_soap=water_in_soap,
                        species_weights=json.dumps(soap.get("species_weights") or {}),
                        **{k: v for k, v in soap.items() if k not in ("weighting", "species_weights")})
    log(f"library: {len(Q)} entries from {len(structs)} structures, species {species} -> {out}")
    return out


def check_library(path, log=print):
    """developer check: is the library one-to-one?  For every entry, the most similar OTHER entry of the same element
    (excluding identical copies from the same structure is not needed: every entry is its own atom) must carry the same
    charge.  Returns counts per reference system and the confusions."""
    d = np.load(path, allow_pickle=False)
    X, q, el = d["X"], d["q"], d["element"]
    org = d["origin"] if "origin" in d.files else np.array(["library"] * len(q))
    wrong = np.zeros(len(q), bool); nn = np.zeros(len(q), int)
    for e in sorted(set(el)):
        m = np.nonzero(el == e)[0]
        for k in range(0, len(m), 4000):
            b = m[k:k + 4000]
            S = X[b].astype(np.float64) @ X[m].astype(np.float64).T
            S[np.arange(len(b)), np.arange(k, k + len(b))] = -np.inf
            j = m[S.argmax(1)]
            nn[b] = j; wrong[b] = np.abs(q[j] - q[b]) > 1e-6
    rep = {}
    for o in sorted(set(org)):
        mo = org == o
        rep[o] = dict(entries=int(mo.sum()), wrong=int((wrong & mo).sum()))
    conf = {}
    for i in np.nonzero(wrong)[0]:
        key = f"{d['type'][i]}({org[i]}) -> {d['type'][nn[i]]}({org[nn[i]]})"
        conf[key] = conf.get(key, 0) + 1
    log(f"library one-to-one check: {len(q) - wrong.sum()} / {len(q)} entries have a nearest other entry with the same charge")
    for o, r in rep.items():
        log(f"  {o}: {r['entries'] - r['wrong']} / {r['entries']}")
    for k, v in conf.items():
        log(f"  wrong: {k} x{v}")
    return dict(per_origin=rep, confusions=conf, total=int(len(q)), wrong=int(wrong.sum()))


def load_library(path):
    d = np.load(path, allow_pickle=False)
    soap = {k: (float(d[k]) if k in ("r_cut", "sigma") else int(d[k])) for k in SOAP_DEFAULT}
    w = json.loads(str(d["weighting"])) if "weighting" in d.files else {}
    if w:
        soap["weighting"] = w
    sw = json.loads(str(d["species_weights"])) if "species_weights" in d.files else {}
    if sw:
        soap["species_weights"] = sw
    return dict(X=d["X"], q=d["q"], type=d["type"], element=d["element"], species=list(d["species"]), soap=soap,
                water_in_soap=bool(d["water_in_soap"]) if "water_in_soap" in d.files else True)


WATER_TYPES = ("ow", "hw", "OW", "HW")
FRAMEWORK_TYPES = ("ao", "at", "mgo", "mgh", "st", "cao", "cah", "feo", "lio")


def hydroxyl_o(t, extra=()):
    """mask of the hydroxyl O types (oh, ohs) among the type array t"""
    return np.array([x in ("oh", "ohs", *extra) for x in t], bool)


def framework_cation(t, si=True):
    """mask of the framework cations; si=False leaves Si out"""
    base = FRAMEWORK_TYPES if si else tuple(x for x in FRAMEWORK_TYPES if x != "st")
    return np.array([x in base for x in t], bool)


def assign(at, lib, threshold=THRESHOLD, block=4000, exclude_types=(), only=None):
    """charge of every atom = charge of the most similar library entry of the same element.
    exclude_types: library entries of these types are not candidates for the atoms in `only` (a boolean mask over the
    atoms; all atoms when None), e.g. the water entries for the atoms that are not water molecules."""
    el = np.array(at.get_chemical_symbols())
    ok_el = np.isin(el, lib["species"])
    idx = np.nonzero(ok_el)[0]
    sub = at[idx]
    F = soap_vectors(sub, lib["species"], **lib["soap"])
    q = np.full(len(at), np.nan); sim = np.zeros(len(at)); match = np.array([""] * len(at), dtype="U12")
    for e in sorted(set(el[idx])):
        m = np.nonzero(el[idx] == e)[0]; r = np.nonzero(lib["element"] == e)[0]
        if not len(r):
            continue
        r_all = r
        r_ok = r[~np.isin(lib["type"][r], list(exclude_types))] if len(exclude_types) else r
        for k in range(0, len(m), block):
            b = m[k:k + block]
            restricted = np.ones(len(b), bool) if only is None else only[idx[b]]
            for sel, rr in ((restricted & (r_ok is not r_all), r_ok), (~(restricted & (r_ok is not r_all)), r_all)):
                if not sel.any() or not len(rr):
                    continue
                bb = b[sel]
                S = F[bb] @ lib["X"][rr].T
                j = S.argmax(1)
                q[idx[bb]] = lib["q"][rr][j]; sim[idx[bb]] = S[np.arange(len(bb)), j]; match[idx[bb]] = lib["type"][rr][j]
    return q, sim, match, ok_el


# ------------------------------------------------------------------------------------------------- boundary conditions
def parse_boundary(boundary):
    """'p p f' -> (['p', 'p', 'f'], [True, True, False]).  Each axis: p, or one or two of f / s / m (lo and hi face),
    as in the LAMMPS boundary command."""
    b = boundary.split()
    if len(b) != 3 or any(not (x == "p" or (1 <= len(x) <= 2 and set(x) <= set("fsm"))) for x in b):
        raise ValueError(f"boundary must be three of p / f / s / m per axis (e.g. 'p p p', 'p p f'), got {boundary!r}")
    return b, [x == "p" for x in b]


def center_nonperiodic(at, periodic):
    """along every non-periodic axis, move the widest empty gap to the box faces so the material sits in the middle
    (a slab given across the periodic boundary would otherwise be cut).  Returns (Atoms, {axis: gap width in A})."""
    s = at.get_scaled_positions(wrap=True)
    gaps = {}
    for k in range(3):
        if periodic[k]:
            continue
        x = np.sort(s[:, k])
        d = np.diff(np.r_[x, x[0] + 1.0])
        j = int(np.argmax(d))
        s[:, k] = np.mod(s[:, k] - (x[j] + d[j] / 2), 1.0)        # gap centre -> face
        h = at.cell.volume / np.linalg.norm(np.cross(*[at.cell[i] for i in range(3) if i != k]))   # box height along k
        gaps["xyz"[k]] = float(d[j] * h)
    at = at.copy()
    at.set_scaled_positions(s)
    return at, gaps


def kspace_setup(kspace, boundary, periodic, tilted):
    """long-range Coulomb for the boundary: (pair_style coul suffix, settings lines, note or None).
    p p p: kspace as given (PPPM).  Only z non-periodic ('f') in an orthogonal box: PPPM/Ewald with the slab correction.
    Anything else (a tilted box, x or y non-periodic, shrink-wrapped faces): MSM, the LAMMPS solver for non-periodic
    boundaries (PPPM/Ewald slab corrections do not support triclinic boxes).
    MSM is run at 1e-4 with the pair cutoff kept (cutoff/adjust no), whatever the kspace accuracy: on a 320-atom
    pyrophyllite slab its energy is as close to PPPM-slab 1e-6 as PPPM-slab 1e-5 is (2.7 vs 4.1 kcal/mol), at ~18x the
    PPPM cost; MSM 1e-5 lets LAMMPS grow the Coulomb cutoff to ~30 A (neighbour list overflow) or, with the cutoff kept,
    costs ~200x."""
    style = kspace.split()[0]
    if all(periodic):
        return "long", [f"kspace_style {kspace}"], None
    if periodic[0] and periodic[1] and boundary[2] == "f" and not tilted and style in ("pppm", "ewald"):
        return "long", [f"kspace_style {kspace}", "kspace_modify slab 3.0"], f"{style} with the slab correction (slab 3.0)"
    why = "tilted (triclinic) box" if tilted else f"boundary {' '.join(boundary)}"
    return "msm", ["kspace_style msm 1.0e-4", "kspace_modify cutoff/adjust no"], \
        f"msm 1.0e-4 ({why}: PPPM/Ewald cannot handle it; ~10-20x slower than PPPM)"


# ------------------------------------------------------------------------------------------------- topology and writer
def topology_from_types(at, types):
    """bonds and angles from the types (not from distances alone, so that overlapping atoms do not create false bonds):
    each water O takes its two nearest water H, each hydroxyl H its nearest hydroxyl O; M-O-H angles for the framework
    cations within 2.5 A (Al) / 2.6 A (Mg) of a hydroxyl O; H-O-H for water."""
    t = np.array(types)
    pos, cell = at.positions, np.array(at.cell)
    inv = np.linalg.inv(cell)

    def d(i, J):
        v = pos[J] - pos[i]
        f = v @ inv
        v = (f - np.round(f)) @ cell
        return np.linalg.norm(v, axis=-1)
    bonds, angles = [], []
    ow, hw = np.nonzero(t == "ow")[0], np.nonzero(t == "hw")[0]
    free = set(hw.tolist())
    pairs = sorted(((float(dd), o, h) for o in ow for h, dd in zip(hw, d(o, hw)) if dd < 1.5))
    nh = {o: [] for o in ow}
    for dd, o, h in pairs:
        if h in free and len(nh[o]) < 2:
            nh[o].append(h); free.discard(h)
    for o, hs in nh.items():
        for h in hs:
            bonds.append((o, h, "OHw"))
        if len(hs) == 2:
            angles.append((hs[0], o, hs[1], ("HOHw",)))
    ohO = np.nonzero(hydroxyl_o(t, extra=("oaq",)))[0]
    fw = np.nonzero(framework_cation(t))[0]
    for h in np.nonzero(t == "ho")[0]:
        if not len(ohO):
            continue
        dd = d(h, ohO); o = ohO[int(np.argmin(dd))]
        bonds.append((o, h, "OH"))
        dm = d(o, fw) if len(fw) else np.array([])
        cat = [m for m, x in zip(fw, dm) if x < {"Si": 2.0, "Al": 2.5, "Mg": 2.6}.get(ELEMENT[t[m]], 2.6)]
        site = "bulk" if len(cat) >= 2 else "surface"
        for m in cat:
            angles.append((m, o, h, ("MOH", ELEMENT[t[m]], site)))
    unbonded = [int(h) for h in free] + [int(o) for o, hs in nh.items() if len(hs) != 2]
    return bonds, angles, unbonded


def check_consistency(at, q, types):
    """molecule and charge consistency of the assigned types (a safety check after the SOAP assignment):
      water: every water O has exactly two water H, every water H belongs to a water O, each water molecule is neutral;
      hydroxyl: every hydroxyl H (ho) sits within 1.25 A of a hydroxyl O (oh/ohs) that is bonded to a framework cation,
      and every hydroxyl O carries exactly one hydroxyl H.
    Returns a list of problems (empty when consistent), each naming the atoms (1-based)."""
    t = np.array(types)
    pos, cell = at.positions, np.array(at.cell)
    inv = np.linalg.inv(cell)

    def d(i, J):
        v = pos[J] - pos[i]
        f = v @ inv
        return np.linalg.norm((f - np.round(f)) @ cell, axis=-1)
    bonds, _, unbonded = topology_from_types(at, types)
    out = []
    if unbonded:
        out.append(f"{len(unbonded)} water atoms not in a water molecule (atoms {[i + 1 for i in sorted(unbonded)][:10]})")
    nw_o, nw_h = int((t == "ow").sum()), int((t == "hw").sum())
    if nw_h != 2 * nw_o:
        out.append(f"{nw_h} water H for {nw_o} water O (expected {2 * nw_o})")
    mol = {}
    for o, h, k in bonds:
        if k == "OHw":
            mol.setdefault(o, []).append(h)
    bad = [o + 1 for o, hs in mol.items() if len(hs) == 2 and abs(q[o] + q[hs[0]] + q[hs[1]]) > 1e-6]
    if bad:
        out.append(f"{len(bad)} water molecules not neutral (O atoms {bad[:10]})")
    fw = np.nonzero(framework_cation(t))[0]
    ohO = np.nonzero(hydroxyl_o(t))[0]
    nho = {int(o): 0 for o in ohO}
    far = []
    for o, h, k in bonds:
        if k == "OH":
            nho[int(o)] = nho.get(int(o), 0) + 1
            if d(h, np.array([o]))[0] > 1.25:
                far.append(h + 1)
    if far:
        out.append(f"{len(far)} hydroxyl H farther than 1.25 A from any hydroxyl O (atoms {far[:10]})")
    wrong_n = [o + 1 for o, n in nho.items() if n != 1]
    if wrong_n:
        out.append(f"{len(wrong_n)} hydroxyl O without exactly one hydroxyl H (atoms {wrong_n[:10]})")
    if len(ohO) and len(fw):
        lone = [int(o) + 1 for o in ohO if not (d(o, fw) < 2.6).any()]
        if lone:
            out.append(f"{len(lone)} hydroxyl O not bonded to a framework cation (atoms {lone[:10]})")
    return out


def write_system(at, q, types, out_dir, cutoff=10.0, kspace="pppm 1.0e-5", boundary="p p p", reserve_water=False,
                 charges="clayff2004", water="spc"):
    """reserve_water: the water types (ow, hw), the water O-H bond and H-O-H angle types and room for their bonds are
    written even when the structure has no water, so that water can be added in LAMMPS (swelling curve)."""
    bnd, periodic = parse_boundary(boundary)
    bonds, angles, unbonded = topology_from_types(at, types)
    tnames = sorted(set(types) | ({"ow", "hw"} if reserve_water else set()), key=lambda x: (ELEMENT[x], x))
    tid = {x: k + 1 for k, x in enumerate(tnames)}
    bkeys = sorted(set(b[2] for b in bonds) | ({"OHw"} if reserve_water else set())); bid = {k: n + 1 for n, k in enumerate(bkeys)}
    akeys = sorted(set(a[3] for a in angles) | ({("HOHw",)} if reserve_water else set())); aid = {k: n + 1 for n, k in enumerate(akeys)}
    mol = np.ones(len(at), int); nxt = 2
    for o, h, k in bonds:
        if k == "OHw":
            if mol[o] == 1:
                mol[o] = nxt; nxt += 1
            mol[h] = mol[o]
    for i, x in enumerate(types):
        if x in ("Na", "K", "Ca", "Cs", "Cl", "Ba"):
            mol[i] = nxt; nxt += 1
    p, xl = lammps_frame(at, wrap=True)
    lx, ly, lz, xy, xz, yz = p.get_lammps_prism()
    L = ["soapff-cui: ClayFF types/LJ/bonded terms, charges from SOAP", "", f"{len(at)} atoms", f"{len(tnames)} atom types",
         f"{len(bonds)} bonds", f"{len(bkeys)} bond types", f"{len(angles)} angles", f"{len(akeys)} angle types"]
    if reserve_water:
        L += ["2 extra bond per atom", "1 extra angle per atom", "2 extra special per atom"]
    L += ["",
         f"0.0 {lx:.8f} xlo xhi", f"0.0 {ly:.8f} ylo yhi", f"0.0 {lz:.8f} zlo zhi"]
    if any(abs(v) > 1e-8 for v in (xy, xz, yz)):
        L.append(f"{xy:.8f} {xz:.8f} {yz:.8f} xy xz yz")
    L += ["", "Masses", ""] + [f"{tid[x]} {MASS[ELEMENT[x]]:.5f}  # {x}" for x in tnames]
    # image flags: an H bonded across the cell boundary gets the image that puts it next to its O
    img = np.zeros((len(at), 3), int)
    frac = at.get_scaled_positions(wrap=True)
    for o, h, k in bonds:
        img[h] = -np.round(frac[h] - frac[o]).astype(int)
    cut = [(o, h) for o, h, k in bonds if any(img[h, a] and not periodic[a] for a in range(3))]
    if cut:
        o, h = cut[0]
        raise ValueError(f"{len(cut)} O-H bonds cross a non-periodic boundary (e.g. O{o + 1}-H{h + 1}): the structure has no "
                         f"empty gap along that axis; add vacuum or use a periodic boundary")
    coul, kspace_lines, kspace_note = kspace_setup(kspace, bnd, periodic, any(abs(v) > 1e-8 for v in (xy, xz, yz)))
    L += ["", "Atoms  # full", ""] + [f"{i + 1} {mol[i]} {tid[types[i]]} {q[i]:.8f} {xl[i, 0]:.8f} {xl[i, 1]:.8f} {xl[i, 2]:.8f}"
                                      f" {img[i, 0]} {img[i, 1]} {img[i, 2]}" for i in range(len(at))]
    if bonds:
        L += ["", "Bonds", ""] + [f"{k + 1} {bid[b[2]]} {b[0] + 1} {b[1] + 1}" for k, b in enumerate(bonds)]
    if angles:
        L += ["", "Angles", ""] + [f"{k + 1} {aid[a[3]]} {a[0] + 1} {a[1] + 1} {a[2] + 1}" for k, a in enumerate(angles)]
    open(os.path.join(out_dir, "system.data"), "w").write("\n".join(L) + "\n")
    S = [f"# soapff-cui settings: ClayFF (Cygan et al. 2004, Tables 1 and 2) with {WATER_LABELS[water]} water; types in system.data "
         f"from SOAP, charges: {charges}",
         "units real", "atom_style full", f"boundary {' '.join(bnd)}",
         "bond_style harmonic", "angle_style harmonic", f"pair_style lj/cut/coul/{coul} {cutoff}", "pair_modify mix arithmetic tail no",
         "special_bonds lj/coul 0.0 0.0 0.0", "read_data system.data"]
    for x in tnames:
        e, sg = LJ[x]
        # Table 1 as published: epsilon = D0, sigma = R0 / 2^(1/6) written in full double precision (not rounded)
        S.append(f"pair_coeff {tid[x]} {tid[x]} {e!r} {sg!r}  # {x}")
    for k in bkeys:
        S.append(f"bond_coeff {bid[k]} {BOND_OH[0]} {BOND_OH[1]}  # {'O-H (water)' if k == 'OHw' else 'O-H'}")
    for k in akeys:
        K, th = ANGLE_HOH if k[0] in ("HOH", "HOHw") else ANGLE_MOH
        S.append(f"angle_coeff {aid[k]} {K} {th}  # {'-'.join(map(str, k))}")
    S += kspace_lines
    open(os.path.join(out_dir, "system.settings"), "w").write("\n".join(S) + "\n")
    return dict(type_ids=tid, water_bond=bid.get("OHw"), water_angle=aid.get(("HOHw",)), nbonds=len(bonds),
                nangles=len(angles), unbonded=unbonded, nmol=int(nxt - 1), kspace=kspace_lines, kspace_note=kspace_note,
                tilted=any(abs(v) > 1e-8 for v in (xy, xz, yz)))


# ------------------------------------------------------------------------------------------------- the force field
BARO = {"iso": "iso", "aniso": "aniso", "tri": "tri", "z": "z"}   # NPT coupling: iso/aniso/full cell (tri)/normal only (z)


ION_NAMES = ("Na", "K", "Ca", "Cs", "Ba", "Cl")


def analysis_lines(type_ids, every=100):
    """LAMMPS commands for soapff_cui.analysis: mean squared displacement of water O and of each ion type computed by
    LAMMPS itself during the MD (compute msd, written every `every` steps (0.1 ps) by fix ave/time to msd_<group>.dat)."""
    groups = [("Ow", type_ids["ow"])] if "ow" in type_ids else []
    groups += [(n, type_ids[n]) for n in type_ids if n in ION_NAMES]
    L = ["# mean squared displacement by LAMMPS (read by soapff-cui analyze)"] if groups else []
    for g, t in groups:
        L += [f"group a_{g} type {t}", f"compute msd_{g} a_{g} msd",
              f"fix msd_{g} all ave/time 1 1 {every} c_msd_{g}[*] file msd_{g}.dat"]
    return L


def md_protocol(temp, nvt_ps, npt_ps=0.0, press=1.0, baro="aniso", type_names=(), dump_ps=1.0, extra=(), relax_ps=10.0,
                timestep_fs=TIMESTEP_FS):
    """MD part of in.lammps (after the energy minimisation): relaxation NVT for relax_ps (no trajectory, no MSD; the
    step is reset to 0 afterwards), then the production: NVT for nvt_ps, then (npt_ps > 0) NPT at press (atm) with the
    barostat coupling baro.
    'z' (only the box length normal to the layers/surface changes) is the usual choice for slabs and interfaces;
    'aniso' or 'tri' for bulk crystals ('tri' also relaxes the tilt factors)."""
    if baro not in BARO:
        raise ValueError(f"baro must be one of {list(BARO)}")
    L = ["change_box all triclinic"] if (npt_ps > 0 and baro == "tri") else []   # must precede the dump (LAMMPS rule)
    # the dump carries the ClayFF type name (element column, read by OVITO as the particle type name), the wrapped
    # positions (OVITO shows the atoms in the box) with the image flags ix iy iz (the analysis rebuilds the continuous
    # positions for the mean squared displacement exactly) and the charge
    steps = lambda ps: int(round(ps * 1000.0 / timestep_fs))
    every = max(1, steps(dump_ps))
    if relax_ps > 0:
        # the artificial initial structure (water put on the mineral) relaxes before anything is recorded
        L += ["# relaxation of the initial structure: NVT, not recorded", "fix relax all nvt temp {0} {0} 100.0".format(temp),
              f"run {steps(relax_ps)}", "unfix relax", "reset_timestep 0", "# production"]
    L += ["fix md all nvt temp {0} {0} 100.0".format(temp), f"dump d all custom {every} traj.lammpstrj id type element x y z ix iy iz q",
          "dump_modify d sort id" + (" element " + " ".join(type_names) if type_names else "")]
    L += list(extra)
    if nvt_ps > 0:
        L.append(f"run {steps(nvt_ps)}")
    if npt_ps > 0:
        L += ["unfix md",
              f"fix md all npt temp {temp} {temp} 100.0 {BARO[baro]} {press} {press} 1000.0",
              f"run {steps(npt_ps)}"]
    return "\n".join(L)


class FallbackError(ValueError):
    """atoms whose environment is not in the library (make with stop_on_fallback)."""


# ------------------------------------------------------------------------------------------------- swelling curve
SPC_OH, SPC_HOH = 1.0, 109.47                               # flexible SPC water of ClayFF: O-H (A), H-O-H (deg) at rest


def reserve_water_types(out_dir, log=print):
    """rewrite system.data / system.settings of an output folder with the water types reserved (write_system
    reserve_water), from the structure of report.json and the charges and types of charges.csv (in.lammps is kept).
    Returns the write_system info."""
    rep = json.load(open(os.path.join(out_dir, "report.json")))
    at = read_structure(rep["structure"])
    rows = list(csv.DictReader(open(os.path.join(out_dir, "charges.csv"))))
    if len(rows) != len(at):
        raise ValueError(f"charges.csv has {len(rows)} atoms, the structure {len(at)}")
    q = np.array([float(r["charge"]) for r in rows])
    types = [r["clayff_type"] for r in rows]
    kspace = " ".join(rep.get("kspace", [""])[0].split()[1:]) if rep.get("kspace") else "pppm 1.0e-5"
    info = write_system(at, q, types, out_dir, kspace=kspace or "pppm 1.0e-5", boundary=rep.get("boundary", "p p p"),
                        reserve_water=True, charges=rep.get("charges", "clayff2004"), water=rep.get("water_model", "spc"))
    log(f"system.data / system.settings rewritten with the water types reserved (ow {info['type_ids']['ow']}, "
        f"hw {info['type_ids']['hw']}, O-H bond {info['water_bond']}, H-O-H angle {info['water_angle']})")
    return info


def water_molecule(type_ids, bond_type, angle_type, water="spc"):
    """LAMMPS molecule file of one flexible water (charges of the water model `water`, spc) for create_atoms ... mol."""
    h = np.radians(SPC_HOH / 2)
    xs = [(0.0, 0.0, 0.0), (SPC_OH * np.sin(h), SPC_OH * np.cos(h), 0.0), (-SPC_OH * np.sin(h), SPC_OH * np.cos(h), 0.0)]
    t = [type_ids["ow"], type_ids["hw"], type_ids["hw"]]
    qq = [WATER_MODELS[water][0], WATER_MODELS[water][1], WATER_MODELS[water][1]]
    L = [f"# {WATER_LABELS[water]} water, written by soapff-cui", "", "3 atoms", "2 bonds", "1 angles", "", "Coords", ""]
    L += [f"{k + 1} {x:.6f} {y:.6f} {z:.6f}" for k, (x, y, z) in enumerate(xs)]
    L += ["", "Types", ""] + [f"{k + 1} {v}" for k, v in enumerate(t)]
    L += ["", "Charges", ""] + [f"{k + 1} {v:.4f}" for k, v in enumerate(qq)]
    L += ["", "Bonds", "", f"1 {bond_type} 1 2", f"2 {bond_type} 1 3", "", "Angles", "", f"1 {angle_type} 2 1 3"]
    return "\n".join(L) + "\n"


def swelling_schedule(one_to=23, two_to=99, down=True, step_up2=2, step_down=2):
    """water molecules added to the starting structure, state by state: one at a time up to `one_to`, then `step_up2`
    at a time up to `two_to` (Cygan et al. 2004: 1 ... 23, 25 ... 99); then, if down, removed `step_down` at a time
    back to 1 (deswelling).  -> (up, down) lists of water counts (relative to the starting structure)."""
    up = list(range(1, one_to + 1)) + list(range(one_to + step_up2, two_to + 1, step_up2))
    dn = list(range(up[-1] - step_down, 0, -step_down)) if down else []
    return up, dn


def swelling_lines(type_ids, n_nonwater, n_water0, mol_max, temp, press, timestep_fs, tilted, up, down,
                   eq_ps=30.0, avg_ps=20.0, seed=4928459, overlap=2.0):
    """LAMMPS commands of the swelling curve (Cygan et al. 2004): for each water content of `up` (then `down`), water
    molecules are added to (removed from) the interlayer, the structure is relaxed (minimize), then NPT (x, y, z and the
    in-plane tilt xy for a tilted box, else aniso; the interlayer shear tilts xz, yz fixed) for eq_ps and avg_ps more; every 100 fs (step, waters, phase, cell, energy) is appended
    to swelling.dat (phase 1 = the averaging part; soapff analyze -> d001 = lz against the water content).
    The interlayer is found from the clay: the clay is moved to the middle of the cell in z at each state (displace_atoms) and water is
    added between its top and the box top or its bottom and the box bottom (the interlayer across the periodic boundary).
    When no water fits (create_atoms overlap), the box is lengthened by 0.5 A along z and the insertion tried again.
    n_nonwater: atoms that are not water at the start; n_water0: water molecules at the start; mol_max: the largest
    molecule ID of system.data (added molecules get the IDs above it, removed are the last added)."""
    nonwater = {"ow", "hw"} | set(ION_NAMES)
    clay = " ".join(str(v) for k, v in sorted(type_ids.items(), key=lambda kv: kv[1]) if k not in nonwater)
    s_eq = int(round(eq_ps * 1000.0 / timestep_fs))
    s_avg = int(round(avg_ps * 1000.0 / timestep_fs))
    every = max(1, int(round(100.0 / timestep_fs)))
    # the cell: x, y, z and (tilted box) the in-plane tilt xy; the interlayer shear tilts xz, yz stay fixed: with water in
    # the interlayer the layers slide (Cygan et al. 2004: alpha and gamma drift), the tilts run away and fix npt stops
    # ("tilted box too far in one step"); d001 = lz does not depend on them
    baro = (f"x {{P}} {{P}} 1000.0 y {{P}} {{P}} 1000.0 z {{P}} {{P}} 1000.0 xy 0.0 0.0 1000.0 couple none" if tilted
            else "aniso {P} {P} 1000.0")
    T, P = f"{float(temp):g}", f"{float(press):g}"
    n_pre = int(round(2000.0 / timestep_fs))        # 2 ps
    state = [# the clay back to the middle of the cell in z (once per state; a lengthened box leaves it off centre)
             "displace_atoms all move 0 0 $(0.5*(zlo+zhi)-xcm(sw_clay,z)) units box",
             f"minimize 1.0e-4 1.0e-6 2000 20000",
             # the inserted water starts with large forces: a short damped run (limited displacement per step, Langevin)
             # before the barostat, which otherwise moves the cell too far in one step ("tilted box too far")
             f"fix sw_lim all nve/limit 0.1",
             f"fix sw_lan all langevin {T} {T} 50.0 {seed + 3} zero yes",
             f"run {n_pre}",
             "unfix sw_lan", "unfix sw_lim",
             f"fix sw_npt all npt temp {T} {T} 100.0 {baro.format(P=P)}",
             "variable sw_phase equal 0",
             f'fix sw_out all print {every} "${{sw_step}} ${{sw_nw}} ${{sw_phase}} ${{sw_lx}} ${{sw_ly}} ${{sw_lz}} ${{sw_xy}} '
             f'${{sw_xz}} ${{sw_yz}} ${{sw_pe}} ${{sw_temp}}" append swelling.dat screen no',
             f"run {s_eq}",
             "variable sw_phase equal 1",
             f"run {s_avg}",
             "unfix sw_out", "unfix sw_npt",
             "write_dump all custom swelling_states.lammpstrj id type x y z modify append yes"]
    L = [f"# swelling curve (soapff; Cygan et al. 2004): water added to the interlayer ({len(up)} states), "
         + (f"then removed ({len(down)} states); " if down else "")
         + f"each state NPT {eq_ps:g} ps + {avg_ps:g} ps (recorded); d001 = lz in swelling.dat",
         f"velocity all create {T} {seed} mom yes rot no dist gaussian",
         f"timestep {timestep_fs:g}",
         "shell rm -f swelling_states.lammpstrj",
         "molecule h2o water.mol",
         f"group sw_clay type {clay}",
         f"variable sw_nw equal (atoms-{n_nonwater})/3-{n_water0}",
         "variable sw_step equal step", "variable sw_lx equal lx", "variable sw_ly equal ly", "variable sw_lz equal lz",
         "variable sw_xy equal xy", "variable sw_xz equal xz", "variable sw_yz equal yz", "variable sw_pe equal pe",
         "variable sw_temp equal temp",
         'print "# step added_waters phase lx ly lz xy xz yz pe temp" file swelling.dat screen no',
         "variable sw_k equal 0",
         "region sw_lo block INF INF INF INF EDGE 0.1 units box", "region sw_hi block INF INF INF INF 0.2 EDGE units box",
         "region sw_il union 2 sw_lo sw_hi",
         f"variable sw_add index {' '.join(map(str, up))}",
         "label sw_up",
         "label sw_ins",
         "region sw_il delete", "region sw_lo delete", "region sw_hi delete",
         "variable sw_z1 equal $(bound(sw_clay,zmin)-0.3)", "variable sw_z2 equal $(bound(sw_clay,zmax)+0.3)",
         "region sw_lo block INF INF INF INF EDGE $((v_sw_z1>zlo+0.05)*v_sw_z1+(v_sw_z1<=zlo+0.05)*(zlo+0.05)) units box",
         "region sw_hi block INF INF INF INF $((v_sw_z2<zhi-0.05)*v_sw_z2+(v_sw_z2>=zhi-0.05)*(zhi-0.05)) EDGE units box",
         "region sw_il union 2 sw_lo sw_hi",
         "variable sw_k equal $(v_sw_k+1)",
         "variable sw_before equal $(v_sw_nw)",
         f"create_atoms 0 random 1 $({seed}+7*v_sw_k:%.0f) sw_il mol h2o $({seed}+11*v_sw_k:%.0f) overlap {overlap:g} maxtry 500 units box",
         'if "$(v_sw_nw) == ${sw_before}" then "change_box all z delta 0 0.5 units box"',
         'if "$(v_sw_k) > 5000" then "print \'swelling: water could not be added\'" "quit 1"',
         'if "$(v_sw_nw) < ${sw_add}" then "jump SELF sw_ins"'] + state + [
         "next sw_add", "jump SELF sw_up"]
    if down:
        L += [f"variable sw_rem index {' '.join(map(str, down))}",
              "label sw_dn",
              f"group sw_del molecule <> $({mol_max}+v_sw_rem+1:%.0f) $({mol_max}+v_sw_nw:%.0f)",
              "delete_atoms group sw_del mol yes compress no",
              "group sw_del delete"] + state + ["next sw_rem", "jump SELF sw_dn"]
    return L


def make(structure, out_dir, library, threshold=THRESHOLD, cutoff=10.0, kspace="pppm 1.0e-5", temp=300.0, nvt_ps=50.0,
         npt_ps=0.0, press=1.0, baro="aniso", boundary="p p p", dump_ps=1.0, fallback=True, check=True, water_rule=False,
         relax_ps=10.0, timestep_fs=TIMESTEP_FS, md=True, reserve_water=False, stop_on_fallback=False, charges="clayff2004",
         water="spc", log=print):
    """charges: "clayff2004", the only scheme (clayff.CHARGE_SCHEMES): SOAP gives the type of every atom, ClayFF Table 1 (the
    published values) its charge.  water: "spc", the flexible SPC water of ClayFF (the only model).
    md=False: in.lammps stops after the energy minimisation (force field, neighbour list, thermo output, minimisation and
    the walls of fixed faces); the MD (velocities, time step, fixes, computes, runs) is then written by the user (soapff-gui).
    fallback=True: atoms whose best similarity is below `threshold` get the geometric ClayFF type and charge.
    fallback=False: they keep the type and charge of the most similar library entry (same element) and are only
    reported (charges.csv source "soap-low"); only elements absent from the library fall back.
    stop_on_fallback: any atom that the library cannot describe (fallback) stops make with an error before anything is
    written (soapff-gui): the library lacks that environment and should be extended.
    water_rule: water molecules (an O with exactly two H within 1.2 A and no framework cation within 2.4 A, and those H)
    are typed by that rule (ow/hw, charges of the water model) instead of by SOAP; the SOAP assignment is used for everything else."""
    bnd, periodic = parse_boundary(boundary)
    if charges not in CHARGE_SCHEMES or water not in WATER_MODELS:
        raise ValueError(f"charges must be one of {list(CHARGE_SCHEMES)}, water one of {list(WATER_MODELS)}")
    if npt_ps > 0 and not all(periodic):
        raise ValueError(f"NPT needs a periodic box (boundary {boundary}): the pressure of a box with vacuum or open faces "
                         f"is not defined; use NVT only (npt_ps 0)")
    os.makedirs(out_dir, exist_ok=True)
    at = read_structure(structure) if isinstance(structure, str) else structure
    lib = load_library(library)
    role, warn = roles(at)
    for w in warn:
        log("WARNING " + w)
    water_free = not lib["water_in_soap"]                    # library built without water: water by rule, out of the SOAP
    by_rule = (role == "water") if (water_rule or water_free) else np.zeros(len(at), bool)
    if water_free:
        sub = np.nonzero(~by_rule)[0]
        q = np.full(len(at), np.nan); sim = np.zeros(len(at)); match = np.array([""] * len(at), dtype="U12")
        ok_el = np.zeros(len(at), bool)
        if len(sub):
            q[sub], sim[sub], match[sub], ok_el[sub] = assign(at[sub], lib, threshold)
    else:
        # water by rule: the atoms that are not water molecules (mineral, ions) never take a water type from the library
        q, sim, match, ok_el = assign(at, lib, threshold, exclude_types=WATER_TYPES if water_rule else (), only=~by_rule)
    q_fb, types_geo = fallback_charges(at, role)
    if by_rule.any():                                        # water by rule, not by SOAP
        qO, qH = WATER_Q
        q = np.where(by_rule, np.where(np.array(at.get_chemical_symbols()) == "O", qO, qH), q)
        match = np.where(by_rule, np.where(np.array(at.get_chemical_symbols()) == "O", "ow", "hw"), match)
        sim = np.where(by_rule, 1.0, sim); ok_el = ok_el | by_rule
    unknown = ((~ok_el) | np.isnan(q) | ((sim < threshold) if fallback else False)) & ~by_rule
    # ions the library cannot describe (e.g. Cl or Na in solution): formal charge and the ion type, by rule
    el_arr = np.array(at.get_chemical_symbols())
    ion_rule = unknown & (role == "ion") & np.isin(el_arr, [e for e in FORMAL_ION if e in LJ])
    if ion_rule.any():
        q = np.where(ion_rule, [FORMAL_ION.get(e, 0) for e in el_arr], q)
        match = np.where(ion_rule, el_arr, match); sim = np.where(ion_rule, 1.0, sim)
        unknown = unknown & ~ion_rule; by_rule = by_rule | ion_rule
        log(f"{int(ion_rule.sum())} ions not described by the library: formal charge, ion type (by rule)")
    low = (~unknown) & (sim < threshold)                     # nearest entry used although below the threshold
    if stop_on_fallback and unknown.any():
        el_ = at.get_chemical_symbols()
        idx = np.nonzero(unknown)[0]
        detail = ", ".join(f"{el_[i]}{i + 1} ({'not in the library' if not ok_el[i] else f'similarity {sim[i]:.3f}'})"
                           for i in idx[:10]) + (" ..." if len(idx) > 10 else "")
        ids = {f" {i + 1} (" for i in idx}
        overlap = sorted((w for w in warn if "overlapping atoms" in w), key=lambda w: not any(x in w for x in ids))
        if overlap:
            raise FallbackError(f"{len(idx)} atoms are not in the library (fallback): {detail}. The structure has overlapping "
                                f"atoms ({'; '.join(overlap[:3])}{' ...' if len(overlap) > 3 else ''}): correct the structure")
        raise FallbackError(f"{len(idx)} atoms are not in the library (fallback): {detail}; library {os.path.abspath(library)}. "
                            f"The library lacks these environments: add a reference structure with them to the library")
    if low.any():
        warn.append(f"{int(low.sum())} atoms below similarity {threshold} kept the nearest library entry (min {sim[low].min():.3f}); "
                    f"see charges.csv (source soap-low)")
        log("WARNING " + warn[-1])
    q = np.where(unknown, q_fb, q)
    # ClayFF type (LJ, bonded terms) of every atom = type of the matched library entry; geometric typing only for the
    # atoms the library cannot describe
    types = []
    for i in range(len(at)):
        m = TYPE_ALIAS.get(str(match[i]), str(match[i]))
        types.append(m if (not unknown[i] and m in LJ) else types_geo[i])
    tarr = np.array(types)
    # charge of every atom from its type: ClayFF Table 1 and the water model (atoms of a type the table does not list keep the
    # library / fallback charge)
    tq = type_charges(charges, water, types=types)
    q = np.array([tq.get(t, q[i]) for i, t in enumerate(types)], float)
    # cross-check: the type found by SOAP against the bond-graph ClayFF type of the same atom (the rule used for
    # fallback).  A library without the environment can give a similar but wrong type (e.g. obts -> ob) at a
    # similarity above the threshold; such atoms are reported (not changed).
    # Only types with the same protonation are compared (bridging O: ob obts obos obss; OH O: oh ohs; cations): in
    # built initial structures water may overlap the mineral, and then the bond graph, not SOAP, misreads which O
    # carry H (ob -> oh, ow -> oaq, hw -> ho); substitution errors (obts -> ob, ohs -> oh, at -> st ...) are kept.
    geo = np.array([TYPE_ALIAS.get(str(t), str(t)) for t in types_geo])
    cls = lambda t: ("bridge" if t in ("ob", "obts", "obos", "obss") else "oh" if t in ("oh", "ohs") else
                     "water" if t in ("ow", "oaq", "hw") else "h" if t == "ho" else "cation")
    same_class = np.array([cls(a) == cls(b) and cls(a) in ("bridge", "oh", "cation") for a, b in zip(tarr, geo)])
    differ = (~unknown) & same_class & (tarr != geo)
    type_check = dict(compared=int(((~unknown) & same_class).sum()), differ=int(differ.sum()),
                      kinds=dict(collections.Counter(f"{g}->{t}" for g, t in zip(geo[differ], tarr[differ]))))
    if differ.any():
        warn.append(f"{int(differ.sum())} atoms: the SOAP type differs from the bond-graph ClayFF type "
                    f"({', '.join(f'{k} x{v}' for k, v in type_check['kinds'].items())}); the library may lack these "
                    "environments - see charges.csv (column bond_type)")
        log("WARNING " + warn[-1])
    role = np.where(np.isin(tarr, ["ow", "hw"]), "water", np.where(np.isin(tarr, list(FORMAL_ION) + ["Na", "K", "Ca", "Cs", "Cl"]), "ion", "mineral"))
    total = float(q.sum())
    # the charges are not changed to make the system neutral: the net charge is reported (ClayFF itself leaves +0.0002 e per
    # octahedral Mg and -0.0002 e per tetrahedral Al substitution; a larger one comes from fallback atoms, a charged
    # surface or a structure whose formal charge is not zero)
    if abs(total) > 1e-6:
        warn.append(f"net charge {total:+.4f} e (not neutralised; LAMMPS adds a uniform neutralising background for the "
                    f"long-range Coulomb)" + ("" if abs(total) < 0.01 else ": check the structure (counter-ions, protonation) "
                                              "and the fallback atoms"))
        log("WARNING " + warn[-1])
    problems = check_consistency(at, q, types)
    for p_ in problems:
        log("CHECK " + p_)
    if problems and check:
        json.dump(dict(structure=structure if isinstance(structure, str) else "Atoms", consistency=problems),
                  open(os.path.join(out_dir, "report.json"), "w"), indent=1)
        raise ValueError("the assigned types are inconsistent (likely mis-assigned atoms; see report.json): "
                         + "; ".join(problems) + ".  Use a library with more environments (e.g. MD frames), or "
                         "--no-check to write the input anyway")
    warn += [f"consistency: {p_}" for p_ in problems]
    gaps = {}
    if not all(periodic):
        at, gaps = center_nonperiodic(at, periodic)
        for ax, g in gaps.items():
            log(f"boundary {' '.join(bnd)}: widest empty gap along {ax} = {g:.2f} A, moved to the box faces")
            if g < cutoff:
                warn.append(f"the empty gap along the non-periodic {ax} axis is only {g:.2f} A (< cutoff {cutoff} A): atoms on "
                            f"the two faces interacted through the periodic boundary of the given structure, now they do not")
                log("WARNING " + warn[-1])
    info = write_system(at, q, types, out_dir, cutoff=cutoff, kspace=kspace, boundary=" ".join(bnd), reserve_water=reserve_water,
                        charges=charges, water=water)
    if info["kspace_note"]:
        log(f"boundary {' '.join(bnd)}: long-range Coulomb by {info['kspace_note']}")
    # atoms reaching a fixed ('f') face would be lost: reflect them there
    walls = " ".join(f"{ax}{side} EDGE" for ax, b in zip("xyz", bnd) if b != "p"
                     for side, c in zip(("lo", "hi"), b * 2 if len(b) == 1 else b) if c == "f")
    walls = f"fix wall all wall/reflect {walls}\n" if walls else ""
    # a shrink-wrapped ('s') face sits on the outermost atom; between rebuilds atoms move past it and MSM stops
    # ("Out of range atoms").  Rebuilding every step re-shrinks the box before each force call (on a 320-atom slab no
    # slower than 'm' or 'f': MSM dominates the cost).
    neigh_check = "yes"
    if any("s" in b for b in bnd):
        neigh_check = "no"
        log(f"boundary {' '.join(bnd)}: 's' faces with MSM -> neighbour list rebuilt every step (neigh_modify check no)")
    if info["unbonded"]:
        warn.append(f"{len(info['unbonded'])} water O/H could not be paired into molecules (atoms {[i + 1 for i in info['unbonded'][:10]]})")
        log("WARNING " + warn[-1])
    inp = f"""# soapff-cui: ClayFF + SOAP charges.  lmp -in in.lammps
include system.settings
group water type {' '.join(str(info['type_ids'][t]) for t in ('ow', 'hw') if t in info['type_ids'])}
neighbor 2.0 bin
neigh_modify every 1 delay 0 check {neigh_check}
thermo_style custom step temp press pe ke etotal vol density lx ly lz
thermo 1000
# minimisation with a small step limit first (built structures may contain overlapping water)
min_style cg
min_modify dmax 0.01
minimize 1.0e-4 1.0e-6 1000 10000
min_modify dmax 0.1
minimize 1.0e-6 1.0e-6 5000 50000
"""
    if md:
        inp += f"""reset_timestep 0
velocity all create {temp} 4928459 mom yes rot no dist gaussian
timestep {timestep_fs}
{walls}{md_protocol(temp, nvt_ps, npt_ps, press, baro, [t for t, _ in sorted(info['type_ids'].items(), key=lambda kv: kv[1])], dump_ps, analysis_lines(info['type_ids'], max(1, int(round(100.0 / timestep_fs)))), relax_ps, timestep_fs)}
write_data final.data
"""
    elif walls:
        inp += "# walls of the fixed faces (boundary)\n" + walls
    if not info["water_bond"]:
        inp = inp.replace("group water type \n", "")
    open(os.path.join(out_dir, "in.lammps"), "w").write(inp)
    el = at.get_chemical_symbols()
    with open(os.path.join(out_dir, "charges.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "element", "role", "clayff_type", "charge", "source", "matched_library_type", "similarity", "bond_type"])
        for i in range(len(at)):
            w.writerow([i + 1, el[i], role[i], types[i], f"{q[i]:+.6f}", "fallback" if unknown[i] else ("rule" if by_rule[i] else "soap-low" if low[i] else "soap"),
                        "" if unknown[i] else match[i], f"{sim[i]:.5f}", geo[i]])
    per = {}
    for t in sorted(set(types)):
        m = np.array(types) == t
        per[t] = dict(atoms=int(m.sum()), charges=sorted({round(float(x), 4) for x in q[m]}),
                      from_soap=int((m & ~unknown).sum()), fallback=int((m & unknown).sum()), min_similarity=float(sim[m].min()))
    rep = dict(structure=structure if isinstance(structure, str) else "Atoms", atoms=len(at), library=os.path.abspath(library),
               soap=lib["soap"], threshold=threshold, boundary=" ".join(bnd), kspace=info["kspace"], vacuum_gap=gaps,
               charges=charges, water_model=water, charge_sum_before=total, charge_sum=float(q.sum()),
               from_soap=int((~unknown).sum()), fallback=int(unknown.sum()), type_check=type_check, warnings=warn,
               roles={r: int((role == r).sum()) for r in ("mineral", "water", "ion")}, per_type=per,
               consistency=problems, md=dict(temp=temp, relax_ps=relax_ps, nvt_ps=nvt_ps, npt_ps=npt_ps, press=press, baro=baro, dump_ps=dump_ps, timestep_fs=timestep_fs))
    json.dump(rep, open(os.path.join(out_dir, "report.json"), "w"), indent=1)
    log(f"{len(at)} atoms ({rep['roles']}); charges from SOAP {rep['from_soap']}, fallback {rep['fallback']}; "
        f"net charge {total:+.4f} e; wrote {out_dir}/system.data, system.settings, in.lammps, charges.csv, report.json")
    return rep


# ------------------------------------------------------------------------------------------------- OVITO
def open_in_ovito(path):
    """open a soapff-cui output (a directory: traj.lammpstrj, else final.data / system.data; or a file) in OVITO.
    The trajectory carries the ClayFF type names and the SOAP charges, so OVITO colours atoms by type
    (add 'Create bonds' for the bonds).  macOS: the Ovito application; elsewhere: the `ovito` command (PATH, or ~/.soapff-gui/ovito/bin/ovito of install.sh)."""
    import shutil
    import subprocess
    import sys
    if os.path.isdir(path):
        for name in ("traj.lammpstrj", "final.data", "system.data"):
            if os.path.exists(os.path.join(path, name)):
                path = os.path.join(path, name)
                break
        else:
            raise FileNotFoundError(f"no traj.lammpstrj / final.data / system.data in {path}")
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if sys.platform == "darwin":
        cmd = ["open", "-a", "Ovito", path]
    elif shutil.which("ovito") or os.path.isfile(os.path.join(os.environ.get("SOAPFF_HOME") or os.path.expanduser("~/.soapff-gui"),
                                                              "ovito", "bin", "ovito")):
        cmd = [shutil.which("ovito") or os.path.join(os.environ.get("SOAPFF_HOME") or os.path.expanduser("~/.soapff-gui"),
                                                      "ovito", "bin", "ovito"), path]
    else:
        raise RuntimeError("OVITO not found (install it from https://www.ovito.org; on Linux the command `ovito` must be on the PATH)")
    subprocess.run(cmd, check=True)
    return path
