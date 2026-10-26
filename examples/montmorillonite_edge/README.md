# Example: the optional library "Montmorillonite edge (Underwood & Bourg 2020)" (paper: test (7))

The montmorillonite particle with edges on every side (`../library/refs_edge/`: `clay_montmorillonite.data` and `.settings`, a single layer of 2,292 atoms; `clay_montmorillonite_2.data`, the particle twice along c, 4,584 atoms) typed with the optional edge library (`src/soapff_gui/data/library_clayff_edge.npz`), whose force field is taken whole from those files (`core.read_forcefield`).

```
~/.soapff-gui/venv/bin/python run.py [lmp]     # lmp default lmp_serial; packmol must be installed
```

1. `clay_montmorillonite.extxyz`, `clay_montmorillonite_2.extxyz`: the particle files reduced to elements, coordinates, and box; their force field (`out_*/`) against the files: types, charges, every pair_coeff, the bond and angle coefficients, the numbers of bonds and angles.
2. The same particle with the default library (ClayFF): it stops, since its edge atoms are not in that library.
3. `single_point/`: a LAMMPS single point (Ewald 1e-5, cutoff 10 Å, special_bonds 0 0 0) of the original files (`original.data`: a copy with a blank line added after "Angle Coeffs", which LAMMPS needs) and of the written `system.data` / `system.settings`.
4. Minimizations (Ewald 1e-4) of the particle, the particle twice, and the particle with a water film 2-12 Å above it with 4 Na⁺ and 4 Cl⁻ (`particle_water.extxyz`, packmol): the water and Cl⁻ take ClayFF (flexible SPC), the particle and Na⁺ the edge model.

## Results (`results.json`)

| | Atoms | Wrong type | Wrong charge | Fallback | Smallest similarity | pair_coeff differing | Bonds, angles | Minimization stopped by |
|---|---|---|---|---|---|---|---|---|
| particle | 2,292 | 0 | 0 | 0 | 1.000 | none | 288, 624 (as in the file) | energy tolerance |
| particle twice | 4,584 | 0 | 0 | 0 | 1.000 | none | 576, 1,248 | energy tolerance |

Written bonded coefficients: O-H 553.935 kcal/mol/Å², 1.0 Å; M-O-H 30.0191 kcal/mol/rad², 109.47°; edge metal-O-metal 800 kcal/mol/rad², 100° (all as in the files).

Single point (kcal/mol), the original files and the written ones: pe -509528.707277, bonds 7.090368, angles 5917.865731, LJ 38869.102638, real-space Coulomb -335112.681592, k-space -219210.084421, the same in every term to 1e-6.

With the default library the particle stops: 237 atoms are not in the library (the edge O and metals and their neighbours; the first listed, edge O, have similarity 0.78-0.80).

The particle with water: 2,560 water molecules, 4 Na⁺, 4 Cl⁻ (9,980 atoms; 2,268 mineral atoms, 32 ions with the 24 Na⁺ of the particle), fallback 0, type check differ 0, net charge 0; the minimization stopped by the energy tolerance.

`results.json` of this run was written from `run.log` (the last step of `run.py` failed on a numpy integer, since fixed; the values are those of the run).
