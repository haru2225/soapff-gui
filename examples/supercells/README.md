# Example: crystals in supercells larger than the library's (paper: Fig. 1, Table "larger crystals", "MD with LAMMPS")

Two crystals in supercells that differ from the library references, as elements, coordinates, and box only:

- `brucite_5x5x4.extxyz`: brucite Mg(OH)₂, 5 x 5 x 4 unit cells, 500 atoms (the structure of Fig. 1, the workflow);
- `pyrophyllite_3x2x3.extxyz`: pyrophyllite Al₂Si₄O₁₀(OH)₂ (COD 9000809), 3 x 2 x 3 unit cells, 720 atoms.

`run.py` builds soapff-gui's force field for each (bundled library; `out_<name>/`: `structure.extxyz` and the files the GUI writes for it, `system.data`, `system.settings`, `in.lammps` up to the minimization, `charges.csv`, `report.json`) and compares it with the truth: the ClayFF type from the bond graph (`soapff_cui._clayff.ff.ff_types`) and its Table 1 charge ("wrong": a different type or a charge more than 0.01 e off). `md.py` runs soapff-gui's MD input with LAMMPS (minimization, 10 ps NVT relaxation, 10 ps NVT, 40 ps NPT with all cell lengths and tilts at 1 atm; 300 K, 0.5 fs; `md_<name>/`: the input `in.lammps` with `system.data`, `system.settings`, and `log.lammps`, `final.data`; `traj.lammpstrj` is not kept in the repository).

```
~/.soapff-gui/venv/bin/python run.py                      # -> out_<name>/, results.json
~/.soapff-gui/venv/bin/python md.py [mpirun lmp np]       # default: mpirun lmp 2  -> md_<name>/, md_results.json
```

## Results

`results.json`:

| Crystal | Atoms | Wrong | Fallback | Smallest similarity | Type check differ | Net charge (e) |
|---|---|---|---|---|---|---|
| brucite 5x5x4 | 500 | 0 | 0 | 1.000 | 0 | 0 |
| pyrophyllite 3x2x3 | 720 | 0 | 0 | 1.000 | 0 | 0 |

`md_results.json` (Open MPI, 2 processes, PPPM 1e-5):

| Crystal | T during NVT and NPT (K) | Mean T, NVT (K) | Mean pressure, NVT in the experimental box (GPa) | Volume, last 20 ps of NPT, vs the input (%) |
|---|---|---|---|---|
| brucite 5x5x4 | 272-325 | 297.2 | +5.83 | +9.2 |
| pyrophyllite 3x2x3 | 278-327 | 299.9 | +1.90 | +1.9 |

The MD was run again for this example (the paper's run was not kept); the paper gives +5.7 GPa, +9.1 % (brucite) and +1.9 GPa, +2.0 % (pyrophyllite), within the scatter of two runs. These pressures and volume changes are those of ClayFF for these crystals, not of the SOAP assignment.
