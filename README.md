# soapff-GUI

A GUI for clay minerals. You pick a structure file (elements, coordinates, and box only); SOAP assigns the ClayFF type and charge of every atom, the GUI writes the LAMMPS input, runs the MD, and analyzes it. The force field is ClayFF (as tabulated in Cygan et al. 2004) with flexible SPC or rigid SPC/E water.

## Features

- **Automatic force-field generation**: each atom gets the type and charge of the most similar atom in the bundled library (12 bulk crystals including dry Na montmorillonite, and the basal (0 0 1) slab in vacuum of the 11 layered ones: 23 structures, 3,453 rows after pruning duplicates; no edges). The GUI writes `system.data`, `system.settings` (bonds, angles, and LJ parameters), and an `in.lammps` that goes up to the minimization. If any atom has an environment that is not in the library (similarity below 0.95), it stops without building a force field.
- **Charges and water model** (Calculation settings › Library and details, or `--water`): SOAP gives the type of every atom; its charge is the published value of ClayFF Table 1 for that type (Cygan et al. 2004; 4 decimals, the ions as in the table), with no other charge scheme. Water: flexible SPC (ClayFF) or rigid SPC/E (its charges, and `fix shake1 all shake ...` on the water bond and angle written after the minimization; the swelling block redefines it for the water it adds). The charges are not shifted to make the system neutral: the net charge is reported; Table 1 leaves +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al substitution.
- **Adding MD steps**: after the minimization, choose from the lists of LAMMPS fix, compute, and other commands to add them (temperature, pressure, time step, and IDs are filled in automatically). Commands can also be entered by **copy and paste** (select a line and press Cmd+V, or use "Open as text (large, paste)" in the editor to edit the whole input in a large window and press "Apply"). The editor can be enlarged by dragging the border at the bottom or with "Expand".
- **Supercell**: replicates the structure along a, b, and c with OVITO's Replicate, then builds the force field.
- **Surface**: cuts the basal plane of a layered mineral from the current calculation (it uses the ClayFF types of the bulk): the cut goes through the interlayer, no cation-O bond is cut, and the interlayer ions are shared by the two faces. A face that would cut bonds (an edge, any face of quartz) is refused: edges, whose O would need H, are not supported. Then builds the force field.
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

- Library: Na montmorillonite (dry), quartz, brucite, pyrophyllite, muscovite (2M1, 1M), phlogopite, nanpingite, kaolinite, lizardite, and talc. For these crystals, the types and charges of all atoms matched the correct answer, both for supercells of different sizes and for structures after MD at 300 K.
- Water is typed by rule (an O with two H in the shape of a water molecule; SPC or SPC/E charges) and left out of the SOAP, in the library and in the structures: the ClayFF type of a mineral atom does not depend on the water around it, so the same library serves dry crystals, hydrated interlayers and mineral-water interfaces. The ions stay in the SOAP (their positions belong to the clay).
- Surfaces: only interfaces that cut no bond, i.e. the basal plane of a layered mineral (cleaved in the interlayer, the ions shared by the two faces; the atoms keep the ClayFF types of the bulk). Edges (cation-O bonds cut, O protonated by bond valence) are not supported: the GUI and `soapff-cui make-surfaces` refuse those faces, the library has no edge references (no JAEA particle, no edge slabs), and a structure with edges is stopped by the similarity check. The edge procedure (bond-valence protonation, `oe_`/`_e` types, Lammers et al. 2017) is still in `soapff_cui/surfaces.py` (`make_slab(..., allow_bond_cuts=True)`) and `clayff.py`, unused.
- Other minerals: make their bulk reference (LAMMPS data file with ClayFF type names), `soapff-cui make-surfaces` for the basal slab, and `soapff-cui build-library` with all of them.
- The force field is ClayFF (Tables 1 and 2 of Cygan et al. 2004; ions as in Table 1) with flexible SPC or rigid SPC/E water.

## Layout

```
src/soapff_gui/  GUI (app.py), command lists (catalog.py), SSH (remote.py), execution (runner.py), library (data/)
src/soapff_cui/  force-field generation (core.py), analysis (analysis.py). The command is soapff-cui
examples/        the examples above
```

## Interfaces in the library

Only interfaces that cut no bond: the basal plane (0 0 1) of the 11 layered minerals (pyrophyllite, pyrophyllite_2, talc, kaolinite, lizardite, brucite, muscovite 1M, muscovite 2M1, phlogopite 1M, nanpingite, montmorillonite (dry, Na)), as slabs in vacuum, with the 12 bulk references: 23 structures, 3,453 library rows. Every other face with |h|, |k|, |l| <= 1 (and every face of quartz) cuts cation-O bonds wherever it is placed and is refused.

| kind | faces | bonds cut | charges |
|---|---|---|---|
| basal (layered minerals) | 001 | none; the interlayer ions shared by the two faces | the ClayFF types of the bulk |
| edges, quartz surfaces | the others | cation-O bonds | not supported |

A mineral-water interface is the basal slab with water added in the empty part ("Add water").

Slab cell: `soapff-cui make-surfaces --keep-cell` (GUI Surface: "Keep the cell", on by default) adds the vacuum along the lattice vector that leaves the plane, so the in-plane vectors and the shape of the input cell are kept; without it the third vector is the surface normal. `--repeat NA NB NC` (GUI: "Repeat the cell") repeats the bulk cell before the cut. `--periods N|auto` (GUI: "Cut inside the cell", on by default) cuts inside the input cell without changing it: the slab keeps N unit cells (auto: all but one) along the direction out of the plane and the rest of the cell is empty; `--vacuum` is then only an extra length added to the cell. The automatic cut keeps all unit cells but one (and a slab at least 10 A thick); the cell is lengthened by whole periods along that direction only when the empty part would be narrower than 6 A (wider than the SOAP cut-off, so that the two faces do not see each other), and the result says so.
