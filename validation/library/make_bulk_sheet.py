#!/usr/bin/env python3
"""Bulk sheet silicates from /Users/harutokono/clay/noencrypt (POSCAR; elements, positions, cell only), as supercells at
least ~10 A wide along each axis (twice the 5 A SOAP cut-off), like make_bulk.py.
Three files are VASP 4 POSCARs without the element line: their element order was taken from the counts and checked by
the coordination numbers (Si: 4 O; octahedral Al and Mg: 6 O, tetrahedral Al: 4 O; H: 1 O within 1.2 A):
    muscovite_1M  12 12 4 48 8  -> Al Si K O H   (Al: 8 octahedral + 4 tetrahedral, Si: 12 tetrahedral)
    nanpingite     8 24 24 96 16 -> Cs Al Si O H  (Al: 16 octahedral + 8 tetrahedral)
    lizerdite      2 3 4 9      -> Si Mg H O
annite (K Fe2+3 AlSi3O10(OH)2) is not built: ClayFF (Cygan et al. 2004, Table 1) has octahedral Fe3+ only.
    python make_bulk_sheet.py -> <name>_bulk.xyz"""
import os
import numpy as np
import ase.io
SRC = "/Users/harutokono/clay/noencrypt"
HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = {"muscovite_1M": ["Al", "Si", "K", "O", "H"], "nanpingite": ["Cs", "Al", "Si", "O", "H"],
         "lizerdite": ["Si", "Mg", "H", "O"]}
NAMES = ("muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "kaolinite", "lizerdite", "talc", "pyrophyllite")


def read(name):
    f = os.path.join(SRC, name + ".vasp")
    at = ase.io.read(f)
    if name in ORDER:
        counts = [int(x) for x in open(f).read().splitlines()[5].split()]
        at.set_chemical_symbols(sum(([e] * n for e, n in zip(ORDER[name], counts)), []))
    return at


if __name__ == "__main__":
    for name in NAMES:
        at = read(name)
        V = at.get_volume()
        widths = [V / np.linalg.norm(np.cross(at.cell[(k + 1) % 3], at.cell[(k + 2) % 3])) for k in range(3)]
        rep = [max(1, int(np.ceil(10.0 / w - 1e-6))) for w in widths]
        big = at.repeat(rep)
        big.wrap()
        out = name if name != "pyrophyllite" else "pyrophyllite_2"      # pyrophyllite_bulk.xyz (build_ase) is kept as is
        ase.io.write(os.path.join(HERE, f"{out}_bulk.xyz"), big, format="extxyz")
        print(f"{out:16s} {len(big):4d} atoms  {big.get_chemical_formula():28s} repeat {rep}  widths",
              np.round(np.array(widths) * rep, 1))
