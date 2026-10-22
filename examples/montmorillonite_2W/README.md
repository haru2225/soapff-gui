# Example: Na montmorillonite with a 2-layer hydrate (2,502 atoms)

`initial.xyz` is the JAEA initial structure of Na montmorillonite with a 2-layer hydrate (two 2:1 layers with Mg substitution, 30 interlayer Na⁺ ions and 344 water molecules, 30.96 x 26.90 x 25.55 Å), reduced to elements, coordinates, and box only. Some of the water overlaps the clay (water H and clay O are as close as 0.5 Å). It contains neither charges nor types.

Note: `initial.xyz` has two Si-O pairs only 0.83 Å apart (and 24 water molecules inside the clay layer). The calculations here were made with an earlier library that contained this structure itself; the current library leaves it out, and soapff-gui stops on it ("overlapping atoms"). Correct the structure, or start from a dry montmorillonite and use "Add water".

## Procedure (soapff-gui)

1. Start soapff-gui and choose `initial.xyz` with "New calculation". The force field is created, and `in.lammps` goes up to the minimization.
2. In the "MD conditions" section of Settings, set temperature 300 K, relaxation NVT 10 ps, NVT 10 ps, NPT 100 ps, pressure 1 atm, and output interval 1 ps (0.1-0.5 ps for diffusion coefficients).
3. Select the last line of `in.lammps` and add the following from the lists, in this order (each new item is inserted after the selected line).
   - + other > soapff procedures "soapff relaxation (initial velocities, then relaxation NVT)"
   - + fix > "NVT (Nosé-Hoover constant temperature)", + other > "trajectory", + other > soapff procedures "diffusion coefficient MSD", + other > "run MD" (run)
   - Then + other > "unfix", + fix > "NPT (only the direction perpendicular to the plane: interlayer, interface)", + other > "run MD", + other > "write the last structure"
   - Correct the number of steps of each run in the field below (the time step is 0.5 fs, so 10 ps = 20000 and 100 ps = 200000).
4. Press "Start calculation" (a little under 1 hour with about 3 MPI processes; on a laptop, using 4 or more can actually be slower).
5. When it finishes, use "Analyze" to see the diffusion coefficients of water O and Na⁺, the SOVA PDF, RDF, coordination numbers, and S(Q). "Open in OVITO" shows the structure colored by type.

The force field is ClayFF (Tables 1 and 2 of Cygan et al. 2004) with flexible SPC water, and the time step is 0.5 fs.

## Expected results (`expected/`)

Force-field generation (`expected/report.json`):

| Item | Value |
|---|---|
| Atoms | 2,502 (clay 1,440, water 1,032, ions 30) |
| Atoms whose charge and type were assigned by SOAP | 2,502 (fallback 0), minimum similarity 1.000 |
| Types and charges | st +2.1 (288), ao +1.575 (114), mgo +1.36 (30), ob -1.05 (600), obos -1.1808 (120), oh -0.95 (84), ohs -1.0808 (60), ho +0.425 (144), Na +1.0 (30), ow -0.82 (344), hw +0.41 (688) |
| Bonds / angles | 832 / 632 (water H-O-H 344, M-O-H 288) |
| Total charge | +0.006 e, set to 0 by a uniform correction on the clay atoms |

MD (`expected/md.png`, 300 K):
- In the relaxation NVT the temperature is 272-312 K (second-half average 299 K), and in the NVT 293-312 K.
- In NPT the interlayer distance d₀₀₁ (half the box height) grows from 12.8 Å to about 15.3 Å in about 5 ps, and then stays constant at 15.2-15.4 Å. The potential energy decreases over roughly the first 40 ps and is almost constant afterwards.
- The x and y box lengths are fixed, so the overall pressure does not reach 1 atm; it is about 6x10³ atm in the second half.

## Diffusion coefficients (calculation for the paper, `diffusion/`)

This calculation is compared with Greathouse et al. (J. Phys. Chem. C 120, 1640, 2016; high-charge model with a layer charge of 0.75 e), who obtained the diffusion coefficients of the 2-layer hydrate with the same ClayFF and flexible SPC water.

1. `sh diffusion/make_start.sh`: builds the force field from `initial.xyz`, then minimization, relaxation NVT 10 ps, NVT 10 ps, and NPT (perpendicular to the plane, P = 0) 100 ps at 300 K (`diffusion/start/`).
2. `sh diffusion/run_all_T.sh` (`make_T.sh 300`, `330`, `366`): at each temperature, assigns new velocities, equilibrates with NVT 20 ps and NPT (perpendicular to the plane, P = 0) 100 ps, and then continues NPT for four 200 ps segments. `compute msd` is redefined for each segment (the time origin is the start of the segment). Time step 0.5 fs.
3. `python diffusion/analyze_diffusion.py [figure.png]`: for each segment, fits the in-plane MSD with a straight line over 50-200 ps to get D_xy, then computes the mean and standard error over the 4 segments and the Arrhenius activation energy -> `diffusion_summary.json`.

The trajectory (`traj.lammpstrj`) is large and is not included. The MSD files (`msd_Ow.dat`, `msd_Na.dat`), the inputs, the logs, and the last structure (`final.data`) are included.

| | 300 K | 330 K | 366 K | Activation energy |
|---|---|---|---|---|
| Interlayer distance d₀₀₁ (Å) | 15.31 | 15.43 | 15.63 | |
| D_xy of water O (10⁻¹⁰ m²/s) | 11.7 ± 0.9 | 18.3 ± 0.5 | 32.5 ± 1.5 | 14.1 ± 1.1 kJ/mol |
| Greathouse et al. | 13.1 | 24.3 | 39.4 | 15.2 kJ/mol |
| D_xy of Na⁺ | 4.3 ± 0.9 | 10.3 ± 0.7 | 15.3 ± 4.7 | 17.8 ± 3.9 kJ/mol |
| Greathouse et al. | 5.1 | 10.4 | 17.0 | 16.5 kJ/mol |

The values here are 0.75-0.99 times those of Greathouse et al., and the activation energies agree within the errors. Reasons for the differences:

- **Different systems**: the layer charge here is 0.83 e / O₂₀(OH)₄ (octahedral Mg substitution), which has more Na⁺ than the high-charge model of Greathouse et al. (0.75 e). Na⁺ attracts water strongly and slows it down, which pushes both water and Na⁺ toward slower diffusion. The amount of water also differs, 9.6 and 7.9 molecules per O₂₀(OH)₄ (more water pushes toward faster diffusion).
- **Averaging**: here there are four 200 ps segments (800 ps in total) with one time origin per segment, fitted with a straight line over 50-200 ps. Greathouse et al. saved coordinates every 2 fs over 2,000 ps and fitted, over 150-250 ps, an MSD averaged with every frame as a time origin. The statistics here are poorer, and in particular the segment-to-segment scatter of Na⁺ is large (7.2-28.1 at 366 K).
- **Finite size**: neither applies a box-size correction, but the box sizes differ, so the size of the correction would also differ.

The force field (ClayFF, flexible SPC, time step 0.5 fs, NPT conditions) is the same, so the difference is mainly due to the system and the statistics.

When the calculation was redone using only soapff-gui on a structure whose layer charge and water content were matched to Greathouse et al., the results agreed to within 0.86-1.16 times (`greathouse_charge/`). The main cause of the difference was the difference in system composition.
