# Example: NVT and NPT of the 11 bulk library structures

For each of the 11 bulk structures used to build the bundled SOAP library (`src/soapff_gui/data/library_clayff_bulk.npz`; the library also holds the (001) slab of each layered mineral), soapff-gui built a force field from the structure file (ClayFF Table 1 charges as they are, LJ parameters unrounded, flexible SPC water) and ran NVT and NPT. The Na montmorillonite with a 2-layer hydrate (JAEA structure) that was in the earlier 13-structure library is no longer a library structure and was removed from this example; see `examples/montmorillonite_2W/` for a 2-layer hydrate.

| Folder | Mineral |
|---|---|
| `brucite/` | brucite (Mg(OH)₂) |
| `kaolinite/` | kaolinite |
| `lizardite/` | lizardite |
| `montmorillonite_dry_Na/` | dry Na montmorillonite (composition of Cygan et al. 2004, Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆, no water; built in `examples/montmorillonite_swelling/`) |
| `muscovite_1M/`, `muscovite_2M1/` | muscovite 1M, 2M1 |
| `nanpingite/` | nanpingite (Cs mica) |
| `phlogopite_1M/` | phlogopite 1M |
| `pyrophyllite/`, `pyrophyllite_2/` | pyrophyllite (two sources) |
| `talc/` | talc |

In each folder, `<name>.extxyz` is the structure file (elements, coordinates, and box only; written from the library's reference LAMMPS data with the charges and types removed), and `run/` is the soapff-cui output.

## Procedure (soapff-cui; the same as in the GUI)

```
python -m soapff_cui make <name>/<name>.extxyz -o <name>/run --library src/soapff_gui/data/library_clayff_bulk.npz \
    --relax-ps 10 --nvt-ps 50 --npt-ps 50 --baro tri --dump-ps 1      # force field + in.lammps
cd <name>/run && mpirun -np 4 lmp -in in.lammps -log log.lammps        # (2 processes for 240 atoms or fewer)
python -m soapff_cui analyze <name>/run --what msd,sova                # run/analysis/
```

The input is: minimization, initial velocities (300 K), a relaxation NVT of 10 ps (not recorded), NVT 50 ps (the trajectory is written every 1 ps), NPT 50 ps (1 atm, all box components and tilts), and the last structure. The time step is 0.5 fs. The runs were made on an SSH computer (4 MPI processes, 2 for 240 atoms or fewer, Open MPI) with PPPM (1e-5), a few minutes to about ten minutes each.

`run/` contains the inputs (`in.lammps`, `system.data`, `system.settings`, `charges.csv`, `report.json`), the outputs (`log.lammps`, `final.data`), and the analysis (`analysis/`). The trajectory (`traj.lammpstrj`) is large and is not included.

## Results

`python summarize.py` produces `summary.csv` and the table below. Temperature, pressure, volume, and density are averages over the second half of each 50 ps run. The volume change is the change of the NPT second-half average relative to the input box. The PE change is the difference between the averages of the first and last 12.5 ps of the NPT second half.

| Mineral | Atoms | fallback | Type mismatches | T NVT (K) | T NPT (K) | P NPT (atm) | Volume change (%) | Density (g/cm³) | PE change (kcal/mol/atom) |
|---|---|---|---|---|---|---|---|---|---|
| brucite | 240 | 0 | 0 | 299.5 | 292.5 | 115 | +9.42 | 2.172 | -0.0062 |
| kaolinite | 408 | 0 | 0 | 298.9 | 300.9 | -1669 | +0.94 | 2.584 | -0.0305 |
| lizardite | 324 | 0 | 0 | 295.5 | 301.7 | 2388 | +3.21 | 2.501 | -0.0116 |
| montmorillonite_dry_Na | 163 | 0 | 0 | 300.6 | 302.8 | 887 | +0.87 | 2.808 | -0.0159 |
| muscovite_1M | 168 | 0 | 0 | 299.3 | 301.5 | -476 | -7.12 | 2.917 | -0.0289 |
| muscovite_2M1 | 336 | 0 | 0 | 301.2 | 301.1 | -153 | -1.73 | 2.898 | -0.0126 |
| nanpingite | 336 | 0 | 0 | 299.1 | 299.3 | -298 | -4.59 | 3.420 | -0.0132 |
| phlogopite_1M | 176 | 0 | 0 | 300.4 | 301.8 | 64 | +0.85 | 2.766 | -0.0006 |
| pyrophyllite | 320 | 0 | 0 | 301.1 | 297.0 | 687 | +2.03 | 2.759 | +0.0158 |
| pyrophyllite_2 | 320 | 0 | 0 | 301.4 | 302.8 | -555 | +2.17 | 2.755 | -0.0028 |
| talc | 378 | 0 | 0 | 302.5 | 301.4 | -949 | +2.86 | 2.699 | +0.0037 |

- For every mineral, all atoms were assigned by SOAP (fallback 0), and there were 0 mismatches between the SOAP types and the types from the bonding connectivity.
- For every mineral, the temperature stayed within 300 K +/- 5 K, and the potential energy was almost constant in the second half of the NPT run (within about 0.03 kcal/mol/atom).
- The NPT volume change is the difference between the input structure (experimental lattice or DFT structure) and the ClayFF equilibrium lattice. Brucite (+9 %) and muscovite 1M (-7 %) change a lot (the same trend as the 300 K MD in the paper).
- Albite is outside the scope of ClayFF, so it was removed from the library and the examples. Quartz (a 3D framework: every face cuts Si-O bonds, so it has no interface that soapff-gui can make) was removed from the library and from this example too.
- The charges are the values of Table 1 of Cygan et al. as they are, so the net charge is not exactly 0 for the minerals with substitutions (+0.0002 e for the dry Na montmorillonite, -0.0016 e for muscovite 1M and phlogopite, -0.0032 e for muscovite 2M1 and nanpingite; it is the rounding of the table, and LAMMPS neutralizes it with a uniform background).
- The average pressure (P NPT) fluctuates by about +/- 1500 atm over 50 ps (more for smaller crystals).
