# Example: the NaCl droplet of Meng et al. 2022 on montmorillonite (0 0 1) (paper: test (6), Fig. "Meng")

The initial structure of Meng et al., Appl. Clay Sci. 224, 106513 (2022), their Droplet-1: a droplet of 9,171 water molecules with 112 Na⁺ and 112 Cl⁻ (3.83 wt%), about 80 Å across and about 5 Å above the (0 0 1) surface of montmorillonite, in a 209.2 x 217.4 x 201 Å box. Made with soapff-gui twice:

- `paper/`: as in the paper. The (0 0 1) slab of dry Na montmorillonite (`../library/refs_bulk/montmorillonite_dry_Na.data`, pymatgen, the vacuum above it) with c set to 201 Å, repeated 20 x 12 (206.9 x 215.5 Å, 78,240 atoms); the water, Na⁺ and Cl⁻ inside one sphere (packmol `inside sphere`) that holds 9,171 molecules at 1 g/cm³ (radius 40.3 Å), its bottom 5 Å above the surface; the force field (bundled library).
- `gui/`: the same with the GUI's steps only: Surface (0 0 1) with the vacuum above the slab, "Vacuum thickness" 183 Å (pymatgen rounds the cell up to whole layers: c = 204.4 Å), Supercell 20 x 12 x 1, Add water with the same sphere.

```
~/.soapff-gui/venv/bin/python run.py [paper|gui]     # packmol must be installed; ~5 minutes each
```

Each folder keeps the structure `structure.extxyz` (105,977 atoms), the files soapff-gui writes for it (`system.data`, `system.settings`, `in.lammps`, `charges.csv`, `report.json`; `in.lammps` is the GUI's, up to the minimisation) and packmol's input and log (`packmol/`); packmol's xyz files are written by `run.py` but not kept in the repository. The structure is large, so no LAMMPS calculation was run for this test.

## Results (`results.json`)

| | Cell (Å) | Atoms | Placed | Nearest atom to the droplet (Å) | Droplet bottom above the surface (Å) | Fallback | Type check differ | Net charge (e) |
|---|---|---|---|---|---|---|---|---|
| `paper/` | 206.93 x 215.53 x 201.00 | 105,977 | 9,171 water, 112 Na⁺, 112 Cl⁻ | 6.38 | 5.28 | 0 | 0 | +0.096 |
| `gui/` | 206.93 x 215.53 x 204.41 | 105,977 | 9,171 water, 112 Na⁺, 112 Cl⁻ | 5.58 | 5.28 | 0 | 0 | +0.096 |

packmol converged in both (about 30 s). The net charge is the rounding of ClayFF Table 1 over the 480 units of the slab (+0.0002 e each). The times in `results.json` for `gui/` include a period when the computer was asleep.
