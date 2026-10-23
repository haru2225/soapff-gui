# Example: Na montmorillonite with a 2-layer hydrate (1,554 atoms) and its diffusion coefficients

`initial.xyz` is a 3x2x1 supercell (31.12 x 35.95 x 14.76 Å) of the 32-water state of `../montmorillonite_swelling/` (the last structure of the NPT run with 32 water molecules, `gui_run/swelling_states.lammpstrj`), reduced to elements, coordinates, and box only: one 2:1 layer with the composition of Cygan et al. (2004), Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆ per unit (layer charge 0.75 e per O₂₀(OH)₄, 3 Na⁺ per unit), and 8.0 water molecules per O₂₀(OH)₄ (2-layer hydrate, d₀₀₁ ≈ 14.8 Å). It has 960 clay atoms, 18 Na⁺ ions, and 192 water molecules (576 atoms); it contains neither charges nor types.

The earlier version of this example used the JAEA initial structure (2 layers, 2,502 atoms, layer charge 0.83 e). That structure has water overlapping the clay (24 water molecules inside the layer, two Si-O pairs 0.83 Å apart) and is no longer in the library: soapff-gui stops on it (52 atoms below the similarity threshold), so it was replaced by this structure, which the current soapff-gui handles with fallback 0.

## Procedure (soapff-cui)

1. `sh diffusion/make_start.sh`: builds the force field from `initial.xyz` and runs the minimization, a relaxation NVT of 10 ps, NVT 10 ps, and NPT (perpendicular to the plane only, 1 atm) 100 ps at 300 K (`diffusion/start/`). The time step is 0.5 fs, and the force field is ClayFF (Tables 1 and 2 of Cygan et al. 2004: the charges of Table 1 as they are, unrounded LJ parameters) with the flexible SPC water.
2. `sh diffusion/run_all_T.sh` (`make_T.sh 300`, `330`, `366`): at each temperature, assigns new velocities to `start/final.data`, equilibrates with NVT 20 ps and NPT (perpendicular to the plane) 100 ps, and then continues NPT for four 200 ps segments. `compute msd` of the water O and Na⁺ is redefined for each segment (the time origin is the start of the segment).
3. `python diffusion/analyze_diffusion.py [figure.png]`: for each segment, fits the in-plane MSD with a straight line over 50-200 ps to get D_xy, then computes the mean and standard error over the 4 segments and the Arrhenius activation energy -> `diffusion_summary.json`.

Environment variables of the scripts: `PY`, `LMP`, `MPIRUN`, `NP` (number of MPI processes), `MPIOPTS` (extra options of `mpirun`, for example `--bind-to none` when several jobs share a computer), and `RUN_ONLY=1` (`make_start.sh` only: skip the force-field generation). The runs here were made on an SSH computer (Linux, Open MPI) with 3 MPI processes per temperature, the three temperatures at the same time.

The trajectory (`traj.lammpstrj`) is large and is not included. The MSD files (`msd_Ow.dat`, `msd_Na.dat`), the inputs, the logs, and the last structure (`final.data`) are included.

## Force-field generation and the start run (`expected/`)

`expected/report.json`:

| Item | Value |
|---|---|
| Atoms | 1,554 (clay 960, water 576, ions 18) |
| Atoms whose charge and type were assigned by SOAP | 1,554 (fallback 0); the type cross-check differs for 0 of the 882 atoms compared |
| Types and charges | st +2.1 (186), at +1.575 (6), ao +1.575 (84), mgo +1.36 (12), ob -1.05 (408), obts -1.1688 (24), obos -1.1808 (48), oh -0.95 (72), ohs -1.0808 (24), ho +0.425 (96), Na +1.0 (18), ow -0.82 (192), hw +0.41 (384) |
| Bonds / angles | 480 / 384 (water O-H 384 and H-O-H 192; clay O-H 96 and M-O-H 192) |
| Total charge | +0.0012 e (the rounding of Table 1: +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al, 6 units; not shifted; LAMMPS adds a uniform neutralizing background) |

The 18 Na⁺ ions are reported as "not described by the library", and they get the ClayFF ion type and charge by rule. Water is typed by rule, as always.

MD of the start run (`expected/md.png`, 300 K):
- In the relaxation NVT the temperature is 264-316 K (second-half average 298 K), and in the NVT 291-310 K (second-half average 299 K).
- In NPT the interlayer distance d₀₀₁ (the box height) stays at 14.7-15.0 Å (average 14.88 Å over the last 50 ps), and the potential energy is constant within the fluctuation. The x and y box lengths are fixed, so the pressure of the whole system is not controlled in-plane (average -150 atm in the second half).

## Diffusion coefficients

Compared with Greathouse et al. (J. Phys. Chem. C 120, 1640, 2016, Table 3: Na montmorillonite, layer charge 0.75 e per O₂₀(OH)₄, 2-layer hydrate), who used the same ClayFF and flexible SPC water. The layer charge (0.75 e) and the water content (8.0 H₂O per O₂₀(OH)₄) of this system are the same as theirs. `python diffusion/analyze_diffusion.py expected/diffusion.png` gives `diffusion/diffusion_summary.json` and the figure (`expected/diffusion.png`: in-plane MSD of the water O and Na⁺, and the Arrhenius plot).

| | 300 K | 330 K | 366 K | Activation energy |
|---|---|---|---|---|
| Interlayer distance d₀₀₁ (Å) | 14.87 | 14.96 | 15.08 | |
| D_xy of water O (10⁻¹⁰ m²/s) | 12.3 ± 0.6 | 21.4 ± 1.9 | 36.8 ± 1.8 | 15.1 ± 0.04 kJ/mol |
| Greathouse et al. | 13.1 | 24.3 | 39.4 | 15.2 kJ/mol |
| D_xy of Na⁺ | 2.5 ± 0.7 | 4.1 ± 0.8 | 7.3 ± 1.7 | 14.8 ± 0.6 kJ/mol |
| Greathouse et al. | 5.1 | 10.4 | 17.0 | 16.5 kJ/mol |

(mean ± standard error of the four 200 ps segments; the activation energy is the Arrhenius fit of the three means, its error the fit error.)

- **Water**: 0.88-0.94 times the values of Greathouse et al., and the same activation energy (15.1 and 15.2 kJ/mol).
- **Na⁺**: 0.40-0.49 times their values; the activation energy is close (14.8 and 16.5 kJ/mol). The scatter between the segments is large (0.8-3.7 at 300 K, 2.1-9.6 at 366 K): there are only 18 Na⁺ ions and one time origin per 200 ps segment, while Greathouse et al. averaged over every frame of 2,000 ps. The composition also differs in where the charge is: this layer (Cygan et al. 2004) has 6 tetrahedral Al per 18 Na⁺, and a tetrahedral substitution, close to the interlayer, holds Na⁺ more strongly than an octahedral one; the substitution sites of their model were not checked here. The Na⁺ values are therefore not a test of the force field at the level of the water values.
- The force field (ClayFF Table 1 charges as they are, flexible SPC, 0.5 fs, NPT along z) is the one of Greathouse et al.; neither applies a finite-size correction.

The earlier version of this example (JAEA structure, layer charge 0.83 e, 9.6 H₂O per O₂₀(OH)₄, the earlier ClayFF settings) gave 0.75-0.99 times their values for both water and Na⁺; it was removed together with the JAEA structure.
