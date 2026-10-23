# soapff-GUI

A GUI for clay minerals. You pick a structure file (elements, coordinates, and box only); SOAP assigns the ClayFF type and charge of every atom, the GUI writes the LAMMPS input, runs the MD, and analyzes it. The force field is ClayFF (as tabulated in Cygan et al. 2004) with flexible SPC or rigid SPC/E water.

## Features

- **Automatic force-field generation**: each atom gets the type and charge of the most similar atom in the bundled library (12 bulk crystals, dry Na montmorillonite and the JAEA montmorillonite particle with edges, and slabs in vacuum of every face (h k l) with |h|, |k|, |l| ≤ 1 of the 12 crystals and the dry montmorillonite: 169 structures, 9,843 rows after pruning duplicates). The GUI writes `system.data`, `system.settings` (bonds, angles, and LJ parameters), and an `in.lammps` that goes up to the minimization. If any atom has an environment that is not in the library (similarity below 0.95), it stops without building a force field.
- **Charges and water model** (Calculation settings › Library and details, or `--charges` / `--water`): SOAP gives the type of every atom, the chosen scheme its charge. *ClayFF 2004* (Cygan et al., Table 1 cations; every O by the bond-valence rule from them, unrounded, so that it differs from the 4-decimal Table 1 by < 4e-5 e and each substitution carries exactly its formal charge) or *Underwood & Bourg* (octahedral Mg 1.361 and the O next to it -1.081 / -1.181, so that each Mg substitution carries exactly -1 e; with their bond-valence edge charges). Water: flexible SPC (ClayFF) or rigid SPC/E (its charges, and `fix shake1 all shake ...` on the water bond and angle written after the minimization; the swelling block redefines it for the water it adds). The charges are not shifted to make the system neutral: the net charge is reported; it is zero to machine precision whenever the formal charge is (the rounded 4-decimal Table 1 would leave +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al substitution).
- **Adding MD steps**: after the minimization, choose from the lists of LAMMPS fix, compute, and other commands to add them (temperature, pressure, time step, and IDs are filled in automatically). Commands can also be entered by **copy and paste** (select a line and press Cmd+V, or use "Open as text (large, paste)" in the editor to edit the whole input in a large window and press "Apply"). The editor can be enlarged by dragging the border at the bottom or with "Expand".
- **Supercell**: replicates the structure along a, b, and c with OVITO's Replicate, then builds the force field.
- **Surface**: cuts a slab in vacuum of a face (h k l) from the current calculation (it uses the ClayFF types of the bulk), with the vacuum you choose; the basal plane of a layered mineral is cleaved in the interlayer (its ions shared by the two faces), other faces (edges) are healed by bond valence (see below). Then builds the force field.
- **Add water**: lists the empty slices of the structure along z (the vacuum above a surface, an interlayer) with their thickness and area; the chosen ones are filled by packmol at a density (default 1.0 g/cm³) or a number of molecules, at least 2 Å from every atom. Then builds the force field (the water typed by rule, with the chosen water model). The recommended way to make hydrated or interfacial systems: give the dry structure and add the water here.
- **Swelling curve**: runs NPT while adding water molecules to the interlayer, and records the box of each state in `swelling.dat` (same procedure as Cygan et al. 2004; the plot is made by `examples/montmorillonite_swelling/analyze_swelling.py`).
- **Running**: LAMMPS runs on this computer or on another computer registered over SSH (just enter `user@host` to register it).
- **Analysis**: diffusion coefficients (when a `compute msd` is present), and SOVA PDF, RDF, coordination numbers, and S(Q). The SOVA dr, maximum r, Q range and step, X-ray or neutron, number of frames, and the pairs to compute, including pairs split by type (example: K-ob K-obts), can be set on the analysis screen (Settings > Analysis). Structures can be displayed in OVITO ("OVITO" (or the menu item "Open in OVITO")).

## Installation

| OS | Steps |
|---|---|
| macOS | `git clone https://github.com/haru2225/soapff-gui.git && sh soapff-gui/install.sh`, then launch from `~/Applications/soapff-GUI.app` |
| Linux | `sudo apt install python3-venv python3-tk` (Ubuntu example), then the same as macOS. Launch with `~/.soapff-gui/venv/bin/soapff-gui` |
| Windows | Install Python 3.10-3.12 from python.org, run `git clone ...`, then `powershell -ExecutionPolicy Bypass -File soapff-gui\install.ps1`, and launch from the desktop shortcut (not tested. Running over SSH is not available, so to use SSH, install in WSL2 following the Linux steps) |

- **`install.sh` installs everything that is needed automatically** (it does not use `sudo`): the Python libraries (numpy, scipy, ASE, dscribe, matplotlib, ovito), **LAMMPS and Open MPI** (the MPI build from conda-forge, in `~/.soapff-gui/lammps`, about 1 GB; skipped if `lmp` is already installed), **SOVA** (`pip install sovapy`, Python 3.11 or later), and the **OVITO application** (on macOS `brew install --cask ovito`; on Linux from conda.ovito.org into `~/.soapff-gui/ovito`, about 1.5 GB). `SOAPFF_NO_LAMMPS=1` and `SOAPFF_NO_OVITO=1` skip the respective installs.
- The Windows `install.ps1` installs only the Python libraries and SOVA (get LAMMPS from https://packages.lammps.org/windows.html and OVITO from https://www.ovito.org). If they are not found, specify their locations in "Settings > Run on".
- If no Python with Tk is found, run `brew install python-tk@3.12` or `brew install uv` and then run `install.sh` again. To update, run `git pull` (if dependencies were added, run `~/.soapff-gui/venv/bin/python -m pip install -e soapff-gui`).

## Usage

1. Choose a structure file with "New calculation" (or "Supercell" to make a supercell). The force field and an input that goes up to the minimization are created.
2. Add MD steps with "+ fix", "+ compute", and "+ Other" (example: soapff relaxation, then NVT, then NPT, then write the last structure). For the swelling curve, use "+ other > soapff procedures > swelling curve".
3. Choose where to run ("Run on") and press "Start calculation". When the run finishes, the GUI analyzes it and shows the plots and numbers (`analysis/summary.txt`).

## Running on another computer (SSH)

Force-field generation and analysis run on this computer; only LAMMPS runs on the remote side (Linux or macOS).

1. **Remote side**: you can log in over SSH, and `lmp` (with KSPACE and MOLECULE) and `mpirun` are available.
2. **Register**: next to "Run on" in the "Run" section, press **+ Add**, enter the `ssh` destination (example: `taro@192.168.0.10`, or `taro@192.168.0.10:10022` if the port is not 22), and press "Add". When asked for the password, enter it; the GUI connects, checks the remote CPU and LAMMPS, and sets the number of MPI processes (physical cores) and `mpirun` automatically.
3. **Run**: choose that computer under "Run on" and press "Start calculation". The GUI sends the input and starts it with `nohup` (it keeps running if the connection drops); when it finishes, the whole result folder is copied back.
4. **Remove**: choose it under "Run on" and press **- Remove**. Detailed settings (the remote folder, the locations of `lmp` and `mpirun`, commands to run first such as `module load lammps`, and the number of processes) can be changed in "Settings > Run on".

Passwords are not saved (you are asked each time you register, run, or press "Test connection"). If you can log in with a key (`ssh-copy-id`), leave the password empty.

## Examples (`examples/`)

| Folder | Contents |
|---|---|
| `library_nvt_npt/` | NVT and NPT of each of the 13 library structures (inputs, outputs, and analysis made with the GUI) |
| `montmorillonite_2W/` | Diffusion coefficients of Na montmorillonite with a 2-layer hydrate (300, 330, and 366 K; compared with Greathouse et al. 2016). The JAEA initial structure used there has two Si-O pairs 0.83 Å apart; the current library leaves it out, and the GUI stops on it (overlapping atoms) |
| `montmorillonite_swelling/` | How to build the same Na montmorillonite structure as Cygan et al. 2004, and its swelling curve |
| `muscovite_2M1/` | Elastic constants, tension, and shear of muscovite 2M1 (calculations for the paper) |

## Supported structures and limits

- Library: Na montmorillonite (dry; the JAEA particle with edges), quartz, brucite, pyrophyllite, muscovite (2M1, 1M), phlogopite, nanpingite, kaolinite, lizardite, and talc. For these crystals, the types and charges of all atoms matched the correct answer, both for supercells of different sizes and for structures after MD at 300 K.
- Water is typed by rule (an O with two H in the shape of a water molecule; SPC or SPC/E charges) and left out of the SOAP, in the library and in the structures: the ClayFF type of a mineral atom does not depend on the water around it, so the same library serves dry crystals, hydrated interlayers and mineral-water interfaces. The ions stay in the SOAP (their positions belong to the clay).
- Surfaces and edges: the slabs (`soapff-cui make-surfaces`) keep every Si and tetrahedral Al whole; an O that lost cations is healed by its bond valence (White & Zelazny at near-neutral pH, as in the edge procedure of Lammers et al. and Underwood & Bourg): one H on Si-O, Al(tet)-O, Al2-O, Mg2-O; an OH2 left on a single octahedral Al or Mg is removed and its charge given to the metal (one OH stays on a metal with several). The new O types `oe_<cations>` and edge metals `<metal>_e<n>` get their charge by the bond-valence rule of Lammers et al. (2017; Table 2 there for the illite edge), q_O = -2 + Σ (Z_i - q_i)/CN_i, with the ClayFF charges (ClayFF 2004) or with Mg 1.361 (Underwood & Bourg, who also add the M-O-M restraint, 800 kcal/mol/rad², 100°); this rule gives every O type of ClayFF Table 1. It moves charge between atoms but never creates it: the slab is neutral when its formal charge (valences + H - 2 per O) is zero, so after the White & Zelazny protonation the surface proton charge is brought to zero (as Lammers et al.): an H is added to the most undersaturated O (or taken from the most oversaturated one), alternating between the faces. On faces not in the library ((210), (120), (201), (021), (112); 73,643 atoms) 0.05 % of the atoms fall back and 0.06 % get another type (mostly an edge Al with one OH2 more or less, 0.24 e).
- Other minerals: make their bulk reference (LAMMPS data file with ClayFF type names), `soapff-cui make-surfaces` for the slabs, and `soapff-cui build-library` with all of them.
- The force field is ClayFF (Tables 1 and 2 of Cygan et al. 2004; ions as in Table 1) with flexible SPC or rigid SPC/E water.

## Layout

```
src/soapff_gui/  GUI (app.py), command lists (catalog.py), SSH (remote.py), execution (runner.py), library (data/)
src/soapff_cui/  force-field generation (core.py), analysis (analysis.py). The command is soapff-cui
examples/        the examples above
```

## Interfaces in the library

Each of the 12 minerals (pyrophyllite, pyrophyllite_2, talc, kaolinite, lizardite, brucite, muscovite 1M, muscovite 2M1, phlogopite 1M, nanpingite, montmorillonite (dry, Na), quartz) has the 13 faces with |h|,|k|,|l| <= 1 (indices of the cell of its own reference file, signs removed: 100 010 001 1-10 110 10-1 101 01-1 011 1-1-1 1-11 11-1 111) as slabs in vacuum, plus the 13 bulk references and the JAEA montmorillonite particle (hexagonal edges, no face index): 156 slabs, 10,017 library rows.

| kind | faces | bonds cut | charges |
|---|---|---|---|
| basal (layered minerals) | 001 | none; the interlayer ions shared by the two faces | the ClayFF types of the bulk |
| edges (layered minerals) | the other 12 | cation-O bonds cut; O healed by bond valence (H added, OH2 on a single octahedral metal removed) | oe_ / _e types by the bond-valence rule (Lammers et al. 2017, Eq. 1) |
| quartz | all 13 | Si-O bonds cut; Si-OH silanol surfaces | oe_sh |

Faces outside this list are cut by `soapff-cui make-surfaces` / the GUI "Surface" button in the same way.

Slab cell: `soapff-cui make-surfaces --keep-cell` (GUI Surface: "Keep the cell", on by default) adds the vacuum along the lattice vector that leaves the plane, so the in-plane vectors and the shape of the input cell are kept (a (1 0 0) cut lengthens a, b and c stay); without it the third vector is the surface normal. `--repeat NA NB NC` (GUI: "Repeat the cell") repeats the bulk cell before the cut. `--periods N|auto` (GUI: "Cut inside the cell", on by default) cuts inside the input cell without changing it: the slab keeps N unit cells (auto: all but one) along the direction out of the plane and the rest of the cell is empty; `--vacuum` is then only an extra length added to the cell.
