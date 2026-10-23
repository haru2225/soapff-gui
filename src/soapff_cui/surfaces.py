"""Surface (slab in vacuum) references for the SOAP library: the basal plane of a layered mineral.

A face (h k l) is given in the unit cell of the bulk reference (the reference itself may be a supercell of it).  The slab is
cut by pymatgen (SlabGenerator: the oriented cell, the slab of whole planes at least `min_thickness` thick and `vacuum` of
vacuum, every termination of the face), and only the terminations that break no cation-O (or OH O-H) bond are kept
(get_slabs with these bonds and max_broken_bonds=0).  A face that breaks bonds at every termination (an edge) is
refused: nothing is added to the atoms (no H).  Every atom keeps the ClayFF type of the bulk reference (its
charge: ClayFF Table 1, clayff.type_charges).  The interlayer ions are left as pymatgen cuts them: the slab has whole unit
cells, so the cleaved interlayer stays on one face with all its ions and the other face is bare (the slab is neutral
as the bulk, but the two faces differ).  Water of the reference is not used (the library leaves water out of the SOAP).
The net formal charge of the slab is reported."""
import itertools
import math

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list


# cation-O bond lengths (A, upper limits) by ClayFF type, and the formal charges of the ions
BOND_MAX = {"st": 2.0, "at": 2.05, "ao": 2.3, "mgo": 2.45, "mgh": 2.45, "cao": 2.75, "cah": 2.75, "feo": 2.35, "lio": 2.5}
OH_MAX = 1.25                     # O-H of a hydroxyl (ho)
ION_Z = {"Na": 1, "K": 1, "Cs": 1, "Li": 1, "Rb": 1, "Ca": 2, "Ba": 2, "Sr": 2, "Cl": -1}
CATION_Z = {"st": 4, "at": 3, "ao": 3, "mgo": 2, "mgh": 2, "cao": 2, "cah": 2, "feo": 3, "lio": 1}   # formal valences
MAX_INDEX = 1                     # only the faces |h|, |k|, |l| <= 1 of the unit cell (those of the library) can be cut


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


def bonds_of(at, types):
    """the bonds that a cut may not break, for pymatgen: {(element, "O"): length} of the framework cations of `types`
    (the longest length of their ClayFF types) and ("H", "O") of the hydroxyls.  An ion of the same element as a framework
    cation (e.g. Ca2+ next to octahedral Ca) would be counted too: such a structure is refused."""
    el = at.get_chemical_symbols()
    out = {}
    for e, t in zip(el, types):
        if t in BOND_MAX:
            out[(e, "O")] = max(out.get((e, "O"), 0.0), BOND_MAX[t])
    if "ho" in types:
        out[("H", "O")] = OH_MAX
    clash = {e for e, t in zip(el, types) if t == e and e in ION_Z and (e, "O") in out}
    if clash:
        raise ValueError(f"the ions {sorted(clash)} have the element of a framework cation: the bonds cannot be told apart")
    return out


def make_slab(bulk, types, miller_unit, min_thickness=16.0, vacuum=16.0):
    """slab of the bulk reference (Atoms, ClayFF type per atom; no water) for the face miller_unit of its unit cell, cut by
    pymatgen with the bonds of bonds_of() unbroken; returns (Atoms, types, info).  ValueError when every termination of the
    face breaks a bond."""
    from pymatgen.core import Lattice, Structure
    from pymatgen.core.surface import SlabGenerator
    if max(map(abs, miller_unit)) > MAX_INDEX:
        raise ValueError(f"face {tuple(miller_unit)}: only |h|, |k|, |l| <= {MAX_INDEX} can be cut")
    reps = unit_repeats(bulk, types)
    m = tuple(h * n for h, n in zip(miller_unit, reps))               # the same plane in the cell of the reference
    g = math.gcd(math.gcd(abs(m[0]), abs(m[1])), abs(m[2]))
    m = tuple(x // g for x in m)
    st = Structure(Lattice(np.array(bulk.cell)), bulk.get_chemical_symbols(), bulk.positions, coords_are_cartesian=True,
                   site_properties={"clayff_type": list(types)})
    gen = SlabGenerator(st, m, min_thickness, vacuum, center_slab=True, primitive=False)
    found = gen.get_slabs(bonds=bonds_of(bulk, types), max_broken_bonds=0)
    if not found:
        raise ValueError(f"face {tuple(miller_unit)} cuts cation-O bonds wherever it is placed: only faces that cut no bond (the "
                         "basal plane between the layers of a layered mineral, e.g. (0 0 1)) are supported; edges are not")

    def formal_of(sl):
        typ = sl.site_properties["clayff_type"]
        sym = [x.specie.symbol for x in sl]
        return (sum(CATION_Z.get(x, 0) for x in typ) + typ.count("ho")
                + sum(ION_Z[s] for s, x in zip(sym, typ) if x == s and s in ION_Z) - 2 * sym.count("O"))
    best = min(found, key=lambda sl: abs(formal_of(sl)))
    best = best.get_orthogonal_c_slab()                 # c along the normal (across the vacuum): packmol's box for the water
    typ = list(best.site_properties["clayff_type"])
    slab = Atoms([x.specie.symbol for x in best], positions=best.cart_coords, cell=best.lattice.matrix, pbc=True)
    slab.wrap()
    sym = slab.get_chemical_symbols()
    I, J, D = neighbor_list("ijd", slab, 1.5)
    bonded_oh = {(min(i, j), max(i, j)) for i, j, dd in zip(I, J, D) if {sym[i], sym[j]} == {"O", "H"} and dd < 1.1}
    close = [dd for i, j, dd in zip(I, J, D) if i < j and (i, j) not in bonded_oh]
    nrm = np.cross(slab.cell[0], slab.cell[1]); nrm /= np.linalg.norm(nrm)
    z = slab.positions @ nrm
    info = dict(miller=tuple(miller_unit), miller_cell=m, repeats=reps, terminations=len(found), cut_bonds=0, basal=True,
                atoms=len(slab), thickness=round(float(z.max() - z.min()), 2),
                ions=int(sum(1 for s, x in zip(sym, typ) if x == s and s in ION_Z)), formal_charge=int(round(formal_of(best))),
                close_contacts=len(close), min_contact=round(min(close, default=9.9), 2))
    return slab, typ, info


def hkl_name(m):
    return "".join(f"m{-x}" if x < 0 else str(x) for x in m)


def write_surfaces(data, out_dir, max_index=MAX_INDEX, min_thickness=16.0, vacuum=16.0, log=print):
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
            slab, typ, info = make_slab(at, types, m, min_thickness=min_thickness, vacuum=vacuum)
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
        log(f"{tag}: {info['atoms']} atoms (no bond cut; {info['terminations']} bond-free terminations), net formal charge {info['formal_charge']:+d}")
    return summary
