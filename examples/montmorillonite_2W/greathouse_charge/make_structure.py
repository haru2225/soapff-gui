"""Na-montmorillonite 2W with the layer charge and water content of the high-charge model of Greathouse et al. (2016),
made from ../initial.xyz (JAEA, 0.83 e / O20(OH)4, 9.6 H2O / O20(OH)4):
  layer charge 0.75 e / O20(OH)4: 27 Mg of 36 unit cells (14 in one layer, 13 in the other) -> 3 octahedral Mg become Al
      (the Mg closest to another Mg of its layer, so that the rest stay spread), 3 Na+ removed (1 and 2 from the interlayers)
  water 7.9 H2O / O20(OH)4: 284 molecules -> 30 of the 172 of each interlayer removed (random, seed 1)
The substitution stays octahedral as in the JAEA model (where Greathouse et al. put theirs was not checked).
    python make_structure.py -> initial.xyz (elements, positions and box only)"""
import os
import numpy as np
from ase.io import read, write

HERE = os.path.dirname(os.path.abspath(__file__))
at = read(os.path.join(HERE, "..", "initial.xyz"))
el = np.array(at.get_chemical_symbols())
z = at.positions[:, 2] % at.cell[2, 2]
Lz = at.cell[2, 2]
rng = np.random.default_rng(1)
D = at.get_all_distances(mic=True)

layer = (z > Lz / 2).astype(int)                 # octahedral sheets at z ~ 3.3 and ~ 16.1 A; interlayers at ~ 9.7, ~ 22.4 A
mg = np.nonzero(el == "Mg")[0]
change = []
for lay, n in ((0, 1), (1, 2)):                  # 15 -> 14 and 15 -> 13
    idx = [i for i in mg if layer[i] == lay]
    for _ in range(n):
        near = {i: min(D[i, j] for j in idx if j != i) for i in idx}
        i = min(near, key=near.get)
        change.append(i); idx.remove(i)
el[change] = "Al"

# water molecules: O of water = O with two H within 1.2 A; interlayer by height
O = np.nonzero(el == "O")[0]; H = np.nonzero(el == "H")[0]
water = []
for o in O:
    hs = [h for h in H if D[o, h] < 1.2]
    if len(hs) == 2 and z[o] % Lz > 0 and not any(el[j] in ("Si", "Al", "Mg") and D[o, j] < 2.3 for j in range(len(el))):
        water.append((o, hs))
inter = lambda i: 0 if z[i] < Lz / 2 + 0.0 and z[i] > 6 else 1
na = np.nonzero(el == "Na")[0]
remove = []
for k, n in ((0, 1), (1, 2)):
    cand = [i for i in na if inter(i) == k]
    remove += list(rng.choice(cand, n, replace=False))
for k in (0, 1):
    mols = [w for w in water if inter(w[0]) == k]
    for j in rng.choice(len(mols), 30, replace=False):
        remove += [mols[j][0]] + mols[j][1]
keep = np.setdiff1d(np.arange(len(el)), remove)
at.set_chemical_symbols(list(el))
new = at[keep]
write(os.path.join(HERE, "initial.xyz"), new, format="extxyz")
e = np.array(new.get_chemical_symbols())
print("atoms", len(new), {x: int((e == x).sum()) for x in ("Si", "Al", "Mg", "Na", "O", "H")},
      "water", len(water) - 60, "of", len(water), "changed Mg", change)
