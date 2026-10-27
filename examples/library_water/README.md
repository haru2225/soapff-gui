# Example: the libraries with water and aqueous ions (solution libraries) and the two-step typing

The default solution libraries are built from **initial structures only** (what packmol builds: not relaxed, no MD frame):

| File (also in `src/soapff_gui/data/`) | For | Rows | Contents |
|---|---|---|---|
| `library_clayff_water.npz` | ClayFF library | 3,996 | the 3,265 rows of the water-free mineral library (11 bulk layered minerals and their basal (0 0 1) slabs) + 731 rows of the packmol structures below |
| `library_clayff_edge_water.npz` | edge library (Underwood & Bourg 2020) | 3,146 | the two particle files of the edge library + the packmol structures below |

| System | Composition | Atoms |
|---|---|---|
| `water` | 512 H2O, 24.83 Å cube | 1,536 |
| `NaCl`, `KCl`, `CsCl` | 496 H2O + 8 M+ + 8 Cl- | 1,504 |
| `CaCl2`, `BaCl2` | 494 H2O + 6 M2+ + 12 Cl- | 1,500 |
| `if_<mineral>` (11; ClayFF library only) | the library's basal slab with water packed above it (packmol, 1 g/cm3, seed 2) | 692-1,659 |
| particle with a water film (edge library only) | the edge particle, a film 2-12 Å above it (1 g/cm3), 4 Na+ + 4 Cl- (packmol, seed 1; `../montmorillonite_edge/out_particle_water/`) | 9,980 |

Water: flexible SPC; ions: ClayFF Table 1 (also in the edge library: the model of the particle has neither Cl- nor water). Rows of the systems with water are pruned at similarity 0.99, slab rows at 0.999. The 2-layer hydrate (its structure comes from an NPT run) is not in these libraries.

Typing (`soapff_cui.core.make(..., solution_library=...)`, `soapff_gui.runner.prepare` chooses the solution library of the library): step 1 types the mineral atoms with the water-free library on the structure without the water molecules (and without the atoms marked "added"); step 2 types the water, the added ions and the step-1 atoms below the threshold with the solution library on the whole structure (water and ions may only match water or ion types; a mineral atom may match any type except water).

The two new libraries are not copied into this directory (they are in `src/soapff_gui/data/`; `build_library_initial.py` and `build_library_clayff_edge_water.py` write them here again).

The libraries built with the frames of a 300 K MD (the previous default, 7,435 rows; `build_library.py`, `build_refs.py`, `test.py`, `soap_only.py`, `mineral_in_water.py`) are kept as `library_clayff_water_md.npz`; its MD frames are not used by the new libraries.

## Files

- `build_library_initial.py` -> `library_clayff_water.npz`, `build_library_initial.json`; `build_library_clayff_edge_water.py` -> `library_clayff_edge_water.npz`, `build_library_clayff_edge_water.json` (rows per type, one-to-one check).
- `build_refs.py`: the reference systems (inputs; the MD of the previous library); `run_remote.py`: the same MD on another computer over SSH (`soapff_gui.remote`; the password is read from the environment variable `SOAPFF_PW` only). `refs/<system>/`: inputs (the packmol structure of each system is its `system.data`; the MD output is not used by the new libraries).
- `water_common.py`: the nearest row of a library for every atom of a structure (whole-structure SOAP), the true types of a structure with water.
- `initial_only.py` -> `initial_only.json`: the 10 packmol interfaces (`../interfaces_001_water/`) and the 6 shapes (`../interface_shapes/`), tested with the MD library (inwater: the three ways to type the mineral atoms and ions of the interfaces).
- `initial_only_lib.py` -> `initial_only_lib.json`: the same 16 structures against the new `library_clayff_water.npz` (whole-structure SOAP).
- `test_edge_water.py` -> `test_edge_water.json`: the edge particle and the particle twice with a water film (seeds 2 and 3, not in the library) against `library_clayff_edge_water.npz`.

## Results

- One-to-one check (nearest other row of the same element): `library_clayff_water.npz` 3,996 / 3,996; `library_clayff_edge_water.npz` 3,146 / 3,146 rows with the same charge and type. (The previous library with MD frames: 7,433 / 7,435.)
- `initial_only_lib.json`, 16 initial structures (28,311 atoms; 4,037 water O, 422 ions) with the whole-structure SOAP: 0 wrong types; 10 interfaces 0 below 0.95 (minimum water and ions 0.965, mineral 0.986); 6 shapes 5 atoms below 0.95 (cylinder 2, ions shape 3; minimum 0.891, a Cl- that packmol put 2.21 Å from another), all with the right nearest type.
- Mineral atoms and ions of the 10 interfaces (4,328 atoms) with the whole-structure SOAP: water-free library 23 below 0.95 (minimum 0.885) and 1 wrong type; with the water removed 0 and 0 (minimum 0.999); solution library 0 and 0 (minimum 0.974) (`initial_only.json`). The atoms below 0.95 are the hydroxyl H of kaolinite (11) and brucite (7) that point to the water (the nearest water atom 2.2-3.0 Å) and 5 interlayer K+ and Cs+.
- `test_edge_water.json`: soapff-gui's two steps give 0 wrong types and 0 fallback on both structures (9,980 and 12,280 atoms); the whole-structure SOAP against the solution library has 0 wrong and 0 below 0.95 (minimum 0.960 water, 0.973 ions, 0.957 particle); against the water-free edge library 1 and 2 atoms of the particle fall below 0.95.
