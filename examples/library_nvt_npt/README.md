# Example: NVT and NPT of the 13 library structures

For each of the 13 structures used to build the bundled SOAP library (`src/soapff_gui/data/library_clayff_bulk.npz`), soapff-gui built a force field from the structure file and ran NVT and NPT.

| Folder | Mineral |
|---|---|
| `brucite/` | brucite (Mg(OH)₂) |
| `kaolinite/` | kaolinite |
| `lizardite/` | lizardite |
| `montmorillonite_2W/` | Na montmorillonite with a 2-layer hydrate (JAEA structure, including water and Na⁺) |
| `montmorillonite_dry_Na/` | dry Na montmorillonite (composition of Cygan et al. 2004, Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆, no water; built in `examples/montmorillonite_swelling/`) |
| `muscovite_1M/`, `muscovite_2M1/` | muscovite 1M, 2M1 |
| `nanpingite/` | nanpingite (Cs mica) |
| `phlogopite_1M/` | phlogopite 1M |
| `pyrophyllite/`, `pyrophyllite_2/` | pyrophyllite (two sources) |
| `quartz/` | quartz |
| `talc/` | talc |

In each folder, `<name>.extxyz` is the structure file (elements, coordinates, and box only; written from the library's reference LAMMPS data with the charges and types removed), and `run/` is the GUI output.

## Procedure (GUI)

1. Choose `<name>.extxyz` with "New calculation" (the bundled library; this creates the force field and an `in.lammps` that goes up to the minimization).
2. Add the following from the lists (temperature 300 K, pressure 1 atm, time step 0.5 fs).
   - + other > soapff relaxation (initial velocities, then relaxation NVT 10 ps)
   - + other > trajectory (every 1 ps)
   - + fix > NVT, `run 100000` (50 ps), `unfix`
   - + fix > NPT (everything including box tilt: crystals) or NPT (only the direction perpendicular to the plane: montmorillonite_2W), `run 100000` (50 ps), `unfix`
   - + other > write the last structure
3. Press "Start calculation" (on this computer, 4 MPI processes (2 for 240 atoms or fewer), 4-35 minutes). When it finishes, the analysis runs automatically (`run/analysis/`).

`run/` contains the inputs (`in.lammps`, `system.data`, `system.settings`, `charges.csv`, `report.json`), the outputs (`log.lammps`, `final.data`), and the analysis (`analysis/`). The trajectory (`traj.lammpstrj`) is large and is not included.

## Results

`python summarize.py` produces `summary.csv` and the table below. Temperature, pressure, volume, and density are averages over the second half of each 50 ps run. The volume change is the change of the NPT second-half average relative to the input box. The PE change is the difference between the averages of the first and last 12.5 ps of the NPT second half.

| Mineral | Atoms | fallback | Type mismatches | T NVT (K) | T NPT (K) | P NPT (atm) | Volume change (%) | Density (g/cm³) | PE change (kcal/mol/atom) |
|---|---|---|---|---|---|---|---|---|---|
| brucite | 240 | 0 | 0 | 299.8 | 304.1 | 84 | +9.17 | 2.177 | +0.0051 |
| kaolinite | 408 | 0 | 0 | 300.9 | 299.1 | 263 | +0.67 | 2.591 | -0.0197 |
| lizardite | 324 | 0 | 0 | 301.3 | 297.9 | 25 | +3.58 | 2.492 | +0.0128 |
| montmorillonite_2W | 2502 | 0 | 0 | 299.6 | 301.6 | 5825 | +19.88 | 2.132 | +0.0006 |
| montmorillonite_dry_Na | 163 | 0 | 0 | 300.4 | 300.3 | 8591 | +0.61 | 2.815 | -0.0116 |
| muscovite_1M | 168 | 0 | 0 | 297.1 | 298.8 | -513 | -7.09 | 2.916 | -0.0048 |
| muscovite_2M1 | 336 | 0 | 0 | 300.4 | 301.0 | 310 | -1.96 | 2.904 | -0.0252 |
| nanpingite | 336 | 0 | 0 | 299.6 | 298.7 | -35 | -4.49 | 3.417 | -0.0031 |
| phlogopite_1M | 176 | 0 | 0 | 302.0 | 301.6 | 609 | +0.71 | 2.769 | +0.0112 |
| pyrophyllite | 320 | 0 | 0 | 301.3 | 297.6 | 353 | +1.96 | 2.761 | -0.0001 |
| pyrophyllite_2 | 320 | 0 | 0 | 303.3 | 301.2 | -522 | +2.09 | 2.757 | -0.0057 |
| quartz | 162 | 0 | 0 | 299.6 | 298.2 | 1381 | +4.74 | 2.530 | -0.0133 |
| talc | 378 | 0 | 0 | 299.5 | 303.8 | 91 | +2.85 | 2.699 | -0.0039 |

- For every mineral, all atoms were assigned by SOAP (fallback 0), and there were 0 mismatches between the SOAP types and the types from the bonding connectivity.
- For every mineral, the temperature stayed within 300 K +/- 5 K, and the potential energy was almost constant in the second half of the NPT run (within 0.03 kcal/mol/atom).
- The NPT volume change is the difference between the input structure (experimental lattice or DFT structure) and the ClayFF equilibrium lattice. Brucite (+9 %) and muscovite 1M (-7 %) change a lot (the same trend as the 300 K MD in the paper).
- In montmorillonite_2W, the interlayer of the input structure is narrow and water overlaps the clay, so the interlayer expands during NPT and the volume increases by 20 %. montmorillonite_2W and montmorillonite_dry_Na keep the x and y box lengths fixed (NPT only in the direction perpendicular to the plane), so the pressure does not reach 1 atm (the dry one stays at the box relaxed at 0 K).
- Albite is outside the scope of ClayFF, so it was removed from the library and the examples.
- The average pressure (P NPT) fluctuates by about +/- 1500 atm over 50 ps (more for smaller crystals).
