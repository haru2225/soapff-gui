# Example: 8 sheet silicates, supercells and structures after MD (paper: Table "sheet silicates", "MD with LAMMPS")

For each of the 8 sheet silicates of the library (muscovite 2M1 and 1M, phlogopite 1M, nanpingite, kaolinite, lizardite, talc, pyrophyllite from a second source):

- `<mineral>.extxyz`: a supercell one unit cell longer along a than the library reference (elements, coordinates, and box only);
- `md_<mineral>/`: soapff-gui's force field for it and a 300 K MD with LAMMPS (`in.lammps`: minimization, 10 ps NVT relaxation, 10 ps NVT, 40 ps NPT with all cell lengths and tilts at 1 atm, 0.5 fs, frames every 1 ps; `log.lammps`, `traj.lammpstrj`, `final.data`; run with Open MPI, 2 processes, PPPM 1e-5); `supercell.json`: the repeats of the unit cell and the input cell.

`run.py` types the supercell and the last MD frame (the same atoms, disordered by 300 K and with the box changed by the NPT) with the bundled library and compares them with the truth: the ClayFF type from the bond graph of the supercell (`soapff_cui._clayff.ff.ff_types`, same atom order) and its Table 1 charge ("wrong": a different type or a charge more than 0.01 e off). It also averages the lattice parameters over the last 30 ps of the NPT (per unit cell, against the input), and reads the temperatures from `log.lammps`.

```
~/.soapff-gui/venv/bin/python run.py        # a few minutes -> out_<mineral>/, out_<mineral>_md/, results.json
```

## Results (`results.json`)

| Mineral | Types | Atoms | Wrong, supercell | Wrong, last MD frame | Smallest similarity, last frame | Fallback, last frame | a, b, c vs input (%) | α, β, γ vs input (°) | Volume (%) | T, production (K) | Mean T, NVT (K) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| muscovite 2M1 | K ao at ho ob obts oh st | 672 | 0 | 0 | 0.953 | 0 | +0.01 -0.26 -1.49 | +0.02 +0.27 -0.02 | -1.79 | 280-332 | 302.8 |
| muscovite 1M | K ao at ho ob obts oh st | 336 | 0 | 0 | 0.990 | 0 | -1.59 -1.05 -4.82 | -0.04 -0.59 -0.01 | -7.12 | 264-337 | 300.2 |
| phlogopite 1M | K at ho mgh ob obts oh st | 352 | 0 | 0 | 0.977 | 0 | -0.76 -0.58 +2.31 | -0.04 +0.62 -0.01 | +0.75 | 274-336 | 297.7 |
| nanpingite | Cs ao at ho ob obts oh st | 504 | 0 | 0 | 0.969 | 0 | -1.14 -0.82 -2.41 | -0.04 -0.04 +0.00 | -4.30 | 273-323 | 301.1 |
| kaolinite | ao ho ob oh st | 544 | 0 | 0 | 0.962 | 0 | +0.81 +0.01 -0.19 | +4.59 -2.17 +0.50 | +0.93 | 279-327 | 302.2 |
| lizardite | ho mgh ob oh st | 432 | 0 | 0 | 0.976 | 0 | +1.42 +1.40 +0.81 | -0.16 +0.13 -0.00 | +3.68 | 271-330 | 300.3 |
| talc | ho mgh ob oh st | 504 | 0 | 0 | 0.980 | 0 | +0.22 +0.27 +2.70 | +1.84 +0.70 -0.07 | +2.86 | 277-328 | 297.6 |
| pyrophyllite (second source) | ao ho ob oh st | 480 | 0 | 0 | 0.994 | 0 | +0.60 +0.58 +0.20 | -0.15 -2.25 +0.20 | +2.05 | 273-320 | 299.1 |

Every atom of the supercells has similarity 1.000. The layer spacing d001 of muscovite goes from 9.99 Å (input) to 9.83 Å in 2M1 and from 10.27 Å to 9.8 Å in 1M: the -4.8 % of c in 1M comes from its input d001, larger than that of 2M1. These lattice changes are those of ClayFF, not of the SOAP assignment (the types and charges are the truth).

The paper's table gives the smallest similarity of the last frame as typed with the library of that time (muscovite 2M1 0.962, phlogopite 0.976, pyrophyllite 0.995); with the current bundled library it is 0.953, 0.977 and 0.994 (the other five the same), all above the threshold 0.95 and with no wrong atom. The temperatures here are the thermo output every 0.5 ps of the NVT and NPT, and the NVT mean of the 10 ps NVT.
