"""ClayFF pieces used for relaxation and for the LAMMPS output: atom typing, reference charges, O-H topology,
M-O-H / H-O-H angles, parameters and the LAMMPS data/settings writer.

Sources
  ClayFF LJ and charges : ClayCode ClayFF_Fe ffnonbonded.itp (Cygan et al. 2004 values), converted kJ/mol -> kcal/mol, nm -> A
  O-H bond, H-O-H angle : ClayFF flexible SPC (K = 554.1349 kcal/mol/A^2, r0 = 1 A; K = 45.7696 kcal/mol/rad^2, 109.47 deg)
  M-O-H angles          : Pouvreau et al. J. Phys. Chem. C 2017 (Al, Mg) and 2019 Table 5 (bulk/surface, Si-O-H),
                          E = k (theta - theta0)^2  (same convention as LAMMPS angle_style harmonic)
  Na+, K+, Ca2+, Cs+, Cl-: IOD set of Smith, Li, Landry, Merz, Li, JCTC 2023, 19, 2064 (as packaged in ClayCode Ions.ff)
  water                 : SPC (ClayFF, O -0.82 / H +0.41, flexible) or SPC/E (O -0.8476 / H +0.4238, rigid; Berendsen 1987)
"""
import numpy as np

from .bonds import Graph
from .chem import EXCHANGEABLE, FORMAL, MASS

KJ = 4.184
LJ = {  # type: (epsilon kcal/mol, sigma A)
    "st": (7.701e-6 / KJ, 3.302), "at": (7.701e-6 / KJ, 3.302), "ao": (5.564e-6 / KJ, 4.2712),
    "mgo": (3.778e-6 / KJ, 5.2643), "mgh": (3.778e-6 / KJ, 5.2643), "cao": (2.104e-5 / KJ, 5.5667),
    "Na": (0.12171955 / KJ, 2.6103332), "K": (0.71203622 / KJ, 3.1092365), "Ca": (0.34886029 / KJ, 2.8651303),
    "Cs": (1.62938558 / KJ, 3.5635949), "Cl": (2.22399118 / KJ, 3.8522461),
    "ho": (0.0, 0.0), "hw": (0.0, 0.0),
}
for _o in ("ob", "obos", "obts", "obss", "oh", "ohs", "ow", "oaq"):
    LJ[_o] = (0.6502 / KJ, 3.1655)
# alternative ion parameter sets (epsilon kcal/mol, sigma A), selected with write_lammps(..., ions=...)
ION_SETS = {
    "iod": {},                                                   # defaults above (ClayCode Ions.ff)
    # Smith & Dang, J. Chem. Phys. 100 (1994) 3757 (Na+, Cl-); Dang, J. Am. Chem. Soc. 117 (1995) 6954 (K+, Cs+);
    # the set used with SPC/E by Lammers et al. (2017)
    "smith-dang": {"Na": (0.1300, 2.350), "Cl": (0.1000, 4.400), "K": (0.1000, 3.334), "Cs": (0.1000, 3.831)},
}
ELEMENT = {"st": "Si", "at": "Al", "ao": "Al", "mgo": "Mg", "mgh": "Mg", "cao": "Ca", "Na": "Na", "K": "K", "Ca": "Ca",
           "Cs": "Cs", "Cl": "Cl",
           "ho": "H", "hw": "H", "ob": "O", "obos": "O", "obts": "O", "obss": "O", "oh": "O", "ohs": "O", "ow": "O", "oaq": "O",
           "oed": "O", "oab": "O", "omb": "O", "oaa": "O", "oma": "O", "alb": "Al", "ala": "Al", "mgb": "Mg", "mga": "Mg"}
CLAYFF_Q = {"st": 2.1, "at": 1.575, "ao": 1.575, "mgo": 1.36, "mgh": 1.05, "cao": 1.36, "Na": 1.0, "K": 1.0, "Ca": 2.0,
            "Cs": 1.0, "Cl": -1.0,
            "ob": -1.05, "obos": -1.1808, "obts": -1.1688, "obss": -1.2996, "oh": -0.95, "ohs": -1.0808,
            "ow": -0.82, "hw": 0.41, "oaq": -0.82, "ho": 0.425}
BOND_OH = (554.1349, 1.0)
WATER_Q = {"spc": (-0.82, 0.41)}
ANGLE_HOH = (45.7696, 109.47)
# (metal, "bulk" = OH bridging >= 2 framework cations, "surface" = terminal OH / OH2) -> (k kcal/mol/rad^2, theta0 deg)
ANGLE_MOH = {("Al", "bulk"): (15.0, 116.0), ("Al", "surface"): (15.0, 110.0),
             ("Mg", "bulk"): (6.0, 120.0), ("Mg", "surface"): (6.0, 110.0),
             ("Si", "bulk"): (15.0, 100.0), ("Si", "surface"): (15.0, 100.0)}


def ff_types(atoms, g=None):
    """ClayFF-style type of every atom (LJ / bonded typing; the charges used in MD come from the charge model)."""
    g = g or Graph(atoms)
    el = g.el
    role = atoms.arrays.get("role")
    out = [None] * len(el)
    cn_o = [sum(1 for j, _, _ in g.nb[i] if el[j] == "O") for i in range(len(el))]
    for i, e in enumerate(el):                              # pass 1: cations
        if e == "Si":
            out[i] = "st"
        elif e == "Al":
            out[i] = "at" if cn_o[i] <= 4 else "ao"
        elif e in EXCHANGEABLE or e == "Cl":
            out[i] = e
    for i, e in enumerate(el):                              # Mg replacing Al in a dioctahedral sheet = mgo, Mg-only sheet = mgh
        if e == "Mg":
            sub = any(el[k] == "Al" and out[k] == "ao" for j, _, _ in g.nb[i] for k, _, _ in g.nb[j])
            out[i] = "mgo" if sub else "mgh"
    for i, e in enumerate(el):                              # pass 2: O (ClayFF substitution-aware types), then H
        if e == "O":
            nbt = [out[j] for j, _, _ in g.nb[i]]
            nh = [el[j] for j, _, _ in g.nb[i]].count("H")
            fw = [t for t in nbt if t in ("st", "at", "ao", "mgo", "mgh")]
            mg, alt = "mgo" in fw, "at" in fw
            if nh >= 2:
                out[i] = "oaq" if fw else "ow"
            elif nh == 1:
                out[i] = "ohs" if mg else "oh"
            else:
                out[i] = "obss" if (mg and alt) else "obos" if mg else "obts" if alt else "ob"
    for i, e in enumerate(el):
        if e == "H":
            o = g.nb[i][0][0] if g.nb[i] else None
            out[i] = "hw" if o is not None and sum(1 for k, _, _ in g.nb[o] if el[k] == "H") >= 2 else "ho"
        elif out[i] is None:
            raise ValueError(f"no ClayFF type for {e}{i}")
    if role is not None:                                   # water added explicitly keeps SPC types
        for i, r in enumerate(role):
            if r == "water" and el[i] == "O":
                out[i] = "ow"
            elif r == "water" and el[i] == "H":
                out[i] = "hw"
    return out


def clayff_charges(types, target=0.0):
    """ClayFF reference charges; any residual (surfaces are not neutral with bulk types) is spread evenly over the
    non-ion, non-water atoms. Returns (charges, residual before correction)."""
    q = np.array([CLAYFF_Q[t] for t in types], float)
    free = np.array([t not in ("Na", "K", "Ca", "ow", "hw") for t in types])
    res = target - q.sum()
    if free.any():
        q[free] += res / free.sum()
    return q, float(res)


def topology(atoms, types, g=None):
    """O-H bonds and M-O-H / H-O-H angles. Returns bonds [(i, j, key)], angles [(i, j, k, key)] with j the O.
    Water molecules (role "water") get their own bond/angle types ("OHw", "HOHw") so that they can be kept rigid."""
    g = g or Graph(atoms)
    el = g.el
    role = atoms.arrays.get("role")
    is_w = (lambda i: role is not None and role[i] == "water")
    bonds, angles = [], []
    for h, e in enumerate(el):
        if e != "H" or not g.nb[h]:
            continue
        o = min(g.nb[h], key=lambda t: t[2])[0]
        bonds.append((o, h, "OHw" if is_w(o) else "OH"))
        fw = [j for j, _, _ in g.nb[o] if el[j] in ("Si", "Al", "Mg")]
        site = "bulk" if len(fw) >= 2 else "surface"
        for m in fw:
            angles.append((m, o, h, ("MOH", el[m], site)))
    for o, e in enumerate(el):
        if e != "O":
            continue
        hs = [j for j, _, _ in g.nb[o] if el[j] == "H"]
        for a in range(len(hs)):
            for b in range(a + 1, len(hs)):
                angles.append((hs[a], o, hs[b], ("HOHw",) if is_w(o) else ("HOH",)))
    return bonds, angles


def angle_params(key):
    if key[0] in ("HOH", "HOHw"):
        return ANGLE_HOH
    return ANGLE_MOH[(key[1], key[2])]


def lammps_frame(atoms, wrap=True):
    """LAMMPS-aligned cell (a along x, b in xy) and positions; returns (prism, lammps positions).
    wrap=False keeps molecules whole (read_data maps periodic atoms back into the box with consistent image flags)."""
    from ase.calculators.lammps import Prism
    p = Prism(np.array(atoms.cell))
    return p, p.vector_to_lammps(atoms.positions, wrap=wrap)


def write_lammps(atoms, charges, data_path, settings_path, types=None, boundary="p p p", kspace="pppm 1.0e-5",
                 cutoff=12.0, slab=None, title="mincharge", wrap=True, ions="iod"):
    """atom_style full data file + a settings file with pair/bond/angle coefficients (ClayFF + IOD + Pouvreau M-O-H)."""
    g = Graph(atoms)
    types = types or ff_types(atoms, g)
    bonds, angles = topology(atoms, types, g)
    tnames = sorted(set(types), key=lambda t: (ELEMENT[t], t))
    tid = {t: k + 1 for k, t in enumerate(tnames)}
    bkeys = sorted(set(b[2] for b in bonds))
    bid = {k: n + 1 for n, k in enumerate(bkeys)}
    akeys = sorted(set(a[3] for a in angles))
    aid = {k: n + 1 for n, k in enumerate(akeys)}
    # molecule ids: 1 = mineral, each water molecule and each ion its own
    role = atoms.arrays.get("role")
    mol = np.ones(len(atoms), int)
    nxt = 2
    if role is not None:
        wO = [o for o, h, k in bonds if k == "OHw"]
        seen = {}
        for o, h, k in bonds:
            if k == "OHw":
                if o not in seen:
                    seen[o] = nxt; nxt += 1
                mol[o] = mol[h] = seen[o]
        for i, r in enumerate(role):
            if r == "ion":
                mol[i] = nxt; nxt += 1
    p, xl = lammps_frame(atoms, wrap=wrap)
    lx, ly, lz, xy, xz, yz = p.get_lammps_prism()
    L = [title, "", f"{len(atoms)} atoms", f"{len(tnames)} atom types", f"{len(bonds)} bonds", f"{len(bkeys)} bond types",
         f"{len(angles)} angles", f"{len(akeys)} angle types", "",
         f"0.0 {lx:.8f} xlo xhi", f"0.0 {ly:.8f} ylo yhi", f"0.0 {lz:.8f} zlo zhi"]
    if any(abs(v) > 1e-8 for v in (xy, xz, yz)):
        L.append(f"{xy:.8f} {xz:.8f} {yz:.8f} xy xz yz")
    L += ["", "Masses", ""] + [f"{tid[t]} {MASS[ELEMENT[t]]:.5f}  # {t}" for t in tnames]
    L += ["", "Atoms  # full", ""]
    for i in range(len(atoms)):
        L.append(f"{i + 1} {mol[i]} {tid[types[i]]} {charges[i]:.8f} {xl[i, 0]:.8f} {xl[i, 1]:.8f} {xl[i, 2]:.8f}")
    if bonds:
        L += ["", "Bonds", ""] + [f"{k + 1} {bid[t]} {a + 1} {b + 1}" for k, (a, b, t) in enumerate(bonds)]
    if angles:
        L += ["", "Angles", ""] + [f"{k + 1} {aid[a[3]]} {a[0] + 1} {a[1] + 1} {a[2] + 1}" for k, a in enumerate(angles)]
    open(data_path, "w").write("\n".join(L) + "\n")
    S = ["# ClayFF (+ IOD ions, Pouvreau M-O-H bending) settings written by mincharge", "units real", "atom_style full",
         f"boundary {boundary}", "bond_style harmonic", "angle_style harmonic",
         f"pair_style lj/cut/coul/long {cutoff}", "pair_modify mix arithmetic tail no", "special_bonds lj/coul 0.0 0.0 0.0",
         f"read_data {data_path}"]
    for t in tnames:
        e, s = ION_SETS[ions].get(t, LJ[t])
        S.append(f"pair_coeff {tid[t]} {tid[t]} {e:.10g} {s:.6g}  # {t}" + (f" ({ions})" if t in ION_SETS[ions] else ""))
    for k in bkeys:
        S.append(f"bond_coeff {bid[k]} {BOND_OH[0]} {BOND_OH[1]}  # {'O-H (water)' if k == 'OHw' else 'O-H'}")
    for k in akeys:
        K, th = angle_params(k)
        S.append(f"angle_coeff {aid[k]} {K} {th}  # {'-'.join(map(str, k))}")
    S.append(f"kspace_style {kspace}")
    if slab:
        S.append(f"kspace_modify slab {slab}")
    open(settings_path, "w").write("\n".join(S) + "\n")
    return dict(types=tnames, type_ids=tid, angle_types=[str(k) for k in akeys], nbonds=len(bonds), nangles=len(angles),
                prism=p, water_bond=bid.get("OHw"), water_angle=aid.get(("HOHw",)))
