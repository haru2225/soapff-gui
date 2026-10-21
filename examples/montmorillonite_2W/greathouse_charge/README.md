# Diffusion coefficients: Na montmorillonite with a 2-layer hydrate, layer charge and water content matched to Greathouse et al. (2016) (calculation for the paper)

Starting from `../initial.xyz` (JAEA, layer charge 0.83 e / O₂₀(OH)₄, 9.6 water molecules / O₂₀(OH)₄), the structure was adjusted to the composition of the high-charge model of Greathouse et al. (J. Phys. Chem. C 120, 1640, 2016) (0.75 e, 7.9 water molecules), and the calculation was done using only soapff-gui.

1. `python make_structure.py` -> `initial.xyz` (2,319 atoms): 3 octahedral Mg atoms (the ones closest to other Mg atoms) are turned back into Al (leaving 14 and 13 Mg per layer), 3 interlayer Na⁺ ions are removed (leaving 14 and 13 per layer), and 30 water molecules are removed from each interlayer (random, seed 1). The substituted sites stay octahedral (the substitution sites of Greathouse et al. were not checked).
2. Build the force field from `initial.xyz` with soapff-gui (all atoms by SOAP, fallback 0, 0 mismatches in the type cross-check).
   `python check_types.py`: compared with the correct answer (the ClayFF types and standard charges of the JAEA structure; Mg changed to Al becomes ao, and O bonded only to those Mg atoms changes from obos to ob and from ohs to oh), the types and charges of all 2,319 atoms agree (0 errors, minimum similarity 0.952, fallback 0). `kept_atoms.txt` and `mg_to_al.txt` are the indices of the atoms kept by `make_structure.py` (indices in the JAEA structure) and of the Mg atoms changed to Al.
3. Build the input from the lists (temperature T, pressure 0): + other > soapff procedures "soapff relaxation" (initial velocities, 10 ps NVT), + fix "NVT" + run 40000 (20 ps) + unfix, + fix "NPT (only the direction perpendicular to the plane)" + run 200000 (100 ps), + other "trajectory" (every 10 ps), then four times "diffusion coefficient MSD" (change the ID and file name to `msd_s1_...` through `msd_s4_...`) + run 400000 (200 ps) + unfix / uncompute, and finally write_data.
4. Set "Run on" to an SSH computer (8 physical cores, 16 threads) and run 300, 330, and 366 K at the same time with 5, 5, and 6 MPI processes (3 hours each).
5. `python analyze_diffusion.py` -> `diffusion_summary.json`: for each segment, fits the in-plane MSD with a straight line over 50-200 ps to get D_xy, then gives the mean and standard error over the 4 segments and the activation energy.

`T300/`, `T330/`, and `T366/` contain the inputs, logs, MSD files, and last structures (trajectories are not included).

| D_xy (10⁻¹⁰ m²/s) | 300 K | 330 K | 366 K | Activation energy (kJ/mol) |
|---|---|---|---|---|
| Water O (this work) | 13.9 ± 1.1 | 24.7 ± 2.8 | 45.6 ± 2.6 | 16.4 ± 0.4 |
| Water O (Greathouse et al.) | 13.1 | 24.3 | 39.4 | 15.2 |
| Na⁺ (this work) | 5.9 ± 1.4 | 9.0 ± 1.9 | 15.8 ± 1.8 | 13.7 ± 1.1 |
| Na⁺ (Greathouse et al.) | 5.1 | 10.4 | 17.0 | 16.5 |
| Interlayer distance d₀₀₁ (Å) | 14.83 | 14.92 | 15.04 | |

The values are 0.86-1.16 times those of Greathouse et al.; water O agrees within 2.5 standard errors and Na⁺ within 0.8 times. With the JAEA composition unchanged (`../diffusion/`) the ratio was 0.75-0.99, so the main cause of the difference was the system composition (layer charge and water content).
