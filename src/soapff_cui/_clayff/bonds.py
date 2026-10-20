"""Bond graph (with periodic shifts), coordination, bond-valence sums and first-shell environment classes."""
import collections

import numpy as np
from ase.neighborlist import neighbor_list

from .chem import BV_B, BV_R0, CUT, EXCHANGEABLE, FORMAL, FRAMEWORK

ORDER = ("H", "Si", "Al", "Mg", "Ca", "Na", "K", "Cs", "Ti", "Fe")


def bond_list(atoms):
    """Directed bonds (i, j, S, d) between O and cations/H closer than CUT; r_j + S @ cell - r_i is the bond vector.
    Every bond appears twice (i->j and j->i)."""
    el = np.array(atoms.get_chemical_symbols())
    I, J, S, D = neighbor_list("ijSd", atoms, max(CUT.values()))
    ei, ej = el[I], el[J]
    cut = np.zeros(len(I))
    m1 = (ei == "O") & np.isin(ej, list(CUT))
    m2 = (ej == "O") & np.isin(ei, list(CUT))
    cut[m1] = [CUT[e] for e in ej[m1]]
    cut[m2] = [CUT[e] for e in ei[m2]]
    ok = (m1 | m2) & (D < cut)
    return I[ok], J[ok], S[ok], D[ok]


class Graph:
    """Neighbour lists nb[i] = [(j, shift(3 ints), distance)] for the M-O / H-O bonds of a periodic structure."""

    def __init__(self, atoms):
        self.atoms = atoms
        self.el = atoms.get_chemical_symbols()
        self.n = len(atoms)
        self.nb = [[] for _ in range(self.n)]
        for i, j, s, d in zip(*bond_list(atoms)):
            self.nb[i].append((int(j), (int(s[0]), int(s[1]), int(s[2])), float(d)))

    def cn(self, i):
        """number of bonded O (cation) or bonded cations + H (oxygen)"""
        return len(self.nb[i])

    def vec(self, i, j, s):
        cell = np.array(self.atoms.cell)
        return self.atoms.positions[j] + np.array(s) @ cell - self.atoms.positions[i]

    def n_h(self, i):
        return sum(1 for j, _, _ in self.nb[i] if self.el[j] == "H")

    def bvs(self):
        """geometric bond-valence sums (Brese & O'Keeffe); O-H bonds count 1 v.u. on the O side."""
        out = np.zeros(self.n)
        for i, e in enumerate(self.el):
            if e in BV_R0:
                out[i] = sum(np.exp((BV_R0[e] - d) / BV_B) for j, _, d in self.nb[i] if self.el[j] == "O")
            elif e == "O":
                out[i] = sum((np.exp((BV_R0[self.el[j]] - d) / BV_B) if self.el[j] in BV_R0 else 1.0) for j, _, d in self.nb[i])
            elif e == "H":
                out[i] = 1.0 if self.nb[i] else 0.0
        return out

    def pauling(self, i, cn_ref=None):
        """Pauling bond strength z/CN of cation i (cn_ref overrides the coordination number, e.g. the bulk one)."""
        e = self.el[i]
        if e == "H":
            return 1.0
        cn = cn_ref if cn_ref is not None else self.cn(i)
        return FORMAL[e] / max(cn, 1)


def env_classes(at):
    """first-shell classes; neighbours from an explicit periodic neighbour list (correct also for cells smaller than the cutoffs)

    Ported from test4/minerals_dict.py (the DFT environment dictionary used these labels); the only changes are "Cs",
    then "Ti" and "Fe" appended to the neighbour order (absent from the training data, so existing labels are
    unchanged)."""
    el = at.get_chemical_symbols(); N = len(el)
    I, Jn, Dn = neighbor_list("ijd", at, max(3.4, max(CUT.values())))       # 3.4 in test4; widened for Cs-O (3.6)
    nb = [[] for _ in range(N)]
    for i, j, d in zip(I, Jn, Dn):
        if el[i] == el[j]: continue
        c = CUT.get(el[j]) if el[i] == "O" else (CUT.get(el[i]) if el[j] == "O" else None)
        if c is not None and d < c: nb[i].append(int(j))
    cls = []
    for i, e in enumerate(el):
        if e == "O":
            c = collections.Counter(el[j] for j in nb[i]); cls.append("O[" + " ".join(f"{k}{c[k]}" for k in ("H", "Si", "Al", "Mg", "Ca", "Na", "K", "Cs", "Ti", "Fe") if c[k]) + "]")
        elif e == "H":
            o = nb[i][0] if nb[i] else None
            if o is None: cls.append("H[?]")
            else:
                c = collections.Counter(el[j] for j in nb[o] if j != i); cls.append("H-O[" + " ".join(f"{k}{c[k]}" for k in ("H", "Si", "Al", "Mg", "Ca", "Na", "K", "Cs", "Ti", "Fe") if c[k]) + "]")
        elif e == "Si":
            br = sum(1 for o in nb[i] if any(el[j] in ("Si",) and j != i for j in nb[o]))
            sial = sum(1 for o in nb[i] if any(el[j] == "Al" for j in nb[o]))
            cls.append(f"Si[CN{len(nb[i])} OSi{br} OAl{sial}]")
        else:
            cls.append(f"{e}[CN{len(nb[i])}]")
    return cls


def _parse(label):
    head, _, body = label.partition("[")
    items = body.rstrip("]").split()
    return head, items


def class_levels(label):
    """Progressively coarser labels for an environment class, most specific first.

    O[H1 Si1 Al2 Na1] -> O[H1 Si1 Al2] -> O[H1 M3] -> O
    H-O[Al2 Na1]      -> H-O[Al2]      -> H-O[M2]   -> H
    Si[CN4 OSi3 OAl1] -> Si[CN4]                    -> Si
    Al[CN6]                                         -> Al
    """
    head, items = _parse(label)
    levels = [label]
    if head in ("O", "H-O"):
        keep = [it for it in items if not any(it.startswith(x) for x in EXCHANGEABLE)]
        l1 = f"{head}[{' '.join(keep)}]"
        nh = sum(int(it[1:]) for it in keep if it.startswith("H"))
        nm = sum(int(it[len(x):]) for it in keep for x in FRAMEWORK if it.startswith(x))
        l2 = f"{head}[" + " ".join(([f"H{nh}"] if nh else []) + ([f"M{nm}"] if nm else [])) + "]"
        levels += [l1, l2, "O" if head == "O" else "H"]
    elif head == "Si":
        levels += [f"Si[{items[0]}]" if items else "Si", "Si"]
    else:
        levels += [head]
    out = []
    for l in levels:
        if l not in out:
            out.append(l)
    return out


def element_of(label):
    h = label.partition("[")[0]
    return "H" if h == "H-O" else h
