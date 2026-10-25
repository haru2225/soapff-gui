# Example: water on the basal (0 0 1) slabs (paper: test (4), Fig. "interfaces")

The library slab of each of the 10 layered minerals of the figure (`../library/refs_basal/`; the second pyrophyllite source is left out, as in the paper) with its widest empty slice (the vacuum above the slab) filled with water by packmol at 1 g/cm³: `soapff_cui.water`, the same as the GUI's "Add water" (tolerance 2 Å, seed 1). Then soapff-gui's force field (bundled library, flexible SPC water) and the LAMMPS minimization of `in.lammps`.

```
~/.soapff-gui/venv/bin/python run.py [lmp]     # lmp: the LAMMPS executable (default lmp_serial); packmol must be installed
```

Each `<mineral>/` holds `structure.extxyz` (slab + water), the force field (`system.data`, `system.settings`, `in.lammps`, `charges.csv`, `report.json`) and `log.min`. The minimization uses Ewald (1e-4) instead of PPPM, as in the paper: with the LAMMPS on the machine used there, PPPM stopped at the start of the minimization.

## Results (`results.json`)

| Mineral | Slab atoms | Water molecules | Atoms | Nearest atom to the water (Å) | Fallback | Type check differ | Net charge (e) | Minimization stopped by |
|---|---|---|---|---|---|---|---|---|
| Na montmorillonite (dry) | 326 | 122 | 692 | 2.33 | 0 | 0 | +0.0004 | energy tolerance |
| muscovite 2M1 | 336 | 134 | 738 | 2.87 | 0 | 0 | -0.0032 | energy tolerance |
| muscovite 1M | 336 | 141 | 759 | 2.2 | 0 | 0 | -0.0032 | energy tolerance |
| phlogopite 1M | 352 | 144 | 784 | 2.2 | 0 | 0 | -0.0032 | energy tolerance |
| nanpingite | 336 | 146 | 774 | 2.36 | 0 | 0 | -0.0032 | energy tolerance |
| pyrophyllite | 320 | 130 | 710 | 2.2 | 0 | 0 | 0 | energy tolerance |
| kaolinite | 816 | 281 | 1659 | 2.16 | 0 | 0 | 0 | energy tolerance |
| lizardite | 648 | 229 | 1335 | 2.68 | 0 | 0 | 0 | energy tolerance |
| talc | 378 | 157 | 849 | 2.36 | 0 | 0 | 0 | energy tolerance |
| brucite | 480 | 135 | 885 | 2.19 | 0 | 0 | 0 | energy tolerance |

packmol converged in every case. The net charges are those of ClayFF Table 1 (the rounding of the substituted sites; -0.0002 e per tetrahedral Al and +0.0002 e per octahedral Mg: 16 tetrahedral Al in the mica slabs; 4 octahedral Mg and 2 tetrahedral Al in the montmorillonite slab).
