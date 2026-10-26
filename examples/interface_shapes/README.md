# Example: initial structures made with packmol's constraints (paper: test (5), Fig. "shapes")

Six initial structures of water on a basal (0 0 1) surface, each made with the constraints of the GUI's "Add water" (`soapff_cui.water.pack`, tolerance 2 Å, seed 1). The substrate is the (0 0 1) slab of dry Na montmorillonite or muscovite 2M1 as the GUI's Surface cuts it (pymatgen, at least 16 Å of mineral, 40 Å of vacuum above it), repeated 3 x 2 in the a-b plane. The constraints are given in packmol's frame (x' = x - y cot γ, z normal to the slab); the number of water molecules is 1 g/cm³ in the free volume of the region (`water.estimate`, the GUI's "Number from density"). Then soapff-gui's force field (bundled library, flexible SPC water) and a LAMMPS minimization (Ewald 1e-4, as in `../interfaces_001_water/`). The paper's figure shows the last one (`ions/`).

```
~/.soapff-gui/venv/bin/python run.py [lmp]     # lmp default lmp_serial; packmol must be installed
```

Each folder holds `structure.extxyz`, the files soapff-gui writes for it, as the GUI writes them (`system.data`, `system.settings`, `in.lammps` up to the minimization, `charges.csv`, `report.json`), the LAMMPS input and log of the minimization (`in.min`: `in.lammps` with `kspace_style ewald 1.0e-4` after the `include`; `log.min`), and packmol's input and log (`packmol/pack.inp`, `packmol/packmol.log`; the large xyz files packmol writes are not kept). The constraints of each case are in `results.json`.

## Results (`results.json`)

| Folder | Region | Substrate | Placed | Atoms | Nearest atom (Å) | Fallback | Type check differ | Net charge (e) | Minimization stopped by |
|---|---|---|---|---|---|---|---|---|---|
| `sphere/` | droplet: inside sphere (r 10 Å) | Na montmorillonite | 139 water | 2373 | 4.92 | 0 | 0 | +0.0024 | energy tolerance |
| `hemisphere/` | droplet on the surface: inside sphere (r 13 Å) + over plane | muscovite 2M1 | 151 water | 2469 | 2.23 | 0 | 0 | -0.0192 | energy tolerance |
| `cylinder/` | water column to the slab above (periodic image): inside cylinder (r 8 Å) | Na montmorillonite | 291 water | 2829 | 3.12 | 0 | 0 | +0.0024 | energy tolerance |
| `cavity/` | water film 2-22 Å above the surface with a cavity: inside box + outside sphere (r 7 Å) | Na montmorillonite | 697 water | 4047 | 3.3 | 0 | 0 | +0.0024 | energy tolerance |
| `xygauss/` | water film with a mound: inside box + below xygauss | Na montmorillonite | 401 water | 3159 | 3.13 | 0 | 0 | +0.0024 | energy tolerance |
| `ions/` | water film 2-22 Å above the surface with 8 Na⁺ + 8 Cl⁻: inside box | muscovite 2M1 (the paper's figure) | 739 water, 8 Na, 8 Cl | 4249 | 2.21 | 0 | 0 | -0.0192 | energy tolerance |

packmol converged in every case. The net charges are those of ClayFF Table 1 (+0.0002 e per octahedral Mg, -0.0002 e per tetrahedral Al of the substrate).
