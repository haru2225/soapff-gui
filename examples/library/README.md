# Example: the SOAP libraries and their references (paper: references of the library, Fig. "minerals", Table "one-to-one")

The reference structures of soapff-gui's libraries, as LAMMPS data files with the ClayFF type name of each atom type in the Masses comments, and `run.py`, which rebuilds the libraries from them and checks them.

| Folder | Contents |
|---|---|
| `refs_bulk/` | the 11 bulk layered minerals: Na montmorillonite (dry; Cygan et al. composition Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆, 163 atoms, see `../montmorillonite_swelling/`), brucite, pyrophyllite (COD 9000809), and 8 sheet silicates (muscovite 2M1 and 1M, phlogopite 1M, nanpingite, kaolinite, lizardite, talc, pyrophyllite from a second source); supercells with every edge at least 10 Å; ClayFF types from the bond graph, Table 1 charges |
| `refs_basal/` | the basal (0 0 1) slab of each, cut by pymatgen (16 Å of mineral at least, 16 Å of vacuum) |
| `refs_hydrate/` | the Na montmorillonite 2-layer hydrate (JAEA structure, 2,502 atoms, water and Na⁺; ClayFF charges and SPC water): a reference of the initial 12-structure library of the paper, not of the current one |
| `refs_edge/` | the montmorillonite particle with edges on every side (`clay_montmorillonite.data`, `.settings`, and `clay_montmorillonite_2.data`, the particle twice along c): the optional edge library |

The bulk references were built from CIF files and the POSCAR files of the coauthors, typed by the bond graph (soapff_cui._clayff.ff.ff_types); the scripts that did it used an external package and are not included: the data files here are the result.

## What run.py does

```
~/.soapff-gui/venv/bin/python run.py        # about 10 minutes -> results.json, library_clayff_bulk.npz, library_clayff_edge.npz
```

1. The default library (`src/soapff_gui/data/library_clayff_bulk.npz`): `refs_bulk/` + `refs_basal/`, water left out of the SOAP, the slab rows pruned at similarity 0.999 (the bulk rows kept whole).
2. The basal slabs cut again from `refs_bulk/` with `soapff_cui.surfaces.make_slab`, compared with `refs_basal/`.
3. The optional edge library (`library_clayff_edge.npz`): `refs_edge/` only, with its force field read from the data and settings files (`core.read_forcefield`).
4. The initial 12-structure library of the paper (the hydrate, the dry Na montmorillonite and the 10 other bulk crystals, water and ions in the SOAP), and the hydrate alone, at l_max = 3, 4, 5.

Each library is checked by `core.check_library`: for every row, the nearest other row of the same element must carry the same charge.

## Results (`results.json`)

| Library | Structures | Rows | Nearest other row with the same charge | Same as the bundled file |
|---|---|---|---|---|
| default (ClayFF) | 22 | 3,265 | 3,265 / 3,265 | yes |
| Montmorillonite edge | 2 | 2,292 | 2,292 / 2,292 | yes |
| initial, hydrate alone (l_max 3, 4, 5) | 1 | 2,502 | 2,502 / 2,502 at each l_max | — |
| initial, 12 structures (l_max 3, 4, 5) | 12 | 5,671 | 5,671 / 5,671 at each l_max | — |

The 11 basal slabs cut again have the same atoms and types as `refs_basal/` (formal charge 0; the interlayer ions on one face: the top face, except nanpingite, bottom).
