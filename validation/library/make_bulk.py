#!/usr/bin/env python3
"""Bare bulk structures without water (elements, positions, cell only), built with ASE exactly like the bulk of the
reference slabs (../minerals/build_ase.py: CIF via ase.io.read, pyrophyllite hydroxyl H added, brucite from
ase.spacegroup.crystal), to test whether the library built from surfaces in water also covers the dry crystals.
Supercells are at least ~10 A along each axis (twice the 5 A SOAP cut-off).
    python make_bulk.py -> {quartz,albite,brucite,pyrophyllite}_bulk.xyz"""
import os, sys
import ase.io
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "minerals"))
import build_ase as B
for name, rep in (("quartz", (3, 3, 2)), ("albite", (2, 1, 2)), ("brucite", (4, 4, 3)), ("pyrophyllite", (2, 2, 2))):
    at = B.bulk(name).repeat(rep)
    ase.io.write(os.path.join(HERE, f"{name}_bulk.xyz"), at, format="extxyz")
    print(name, len(at), at.get_chemical_formula(), at.cell.lengths().round(2))
