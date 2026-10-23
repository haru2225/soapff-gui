# soapff-GUI

A GUI for clay minerals. You pick a structure file (elements, coordinates, and box only); SOAP assigns the ClayFF type and charge of every atom, the GUI writes the LAMMPS input, runs the MD, and analyzes it. The force field is ClayFF (as tabulated in Cygan et al. 2004) with its flexible SPC water.

## Features

- **Automatic force-field generation**: each atom gets the type and charge of the most similar atom in the bundled library (11 bulk layered minerals including dry Na montmorillonite, and the basal (0 0 1) slab in vacuum of each: 22 structures, 3,265 rows after pruning duplicates; no edges). The GUI writes `system.data`, `system.settings` (bonds, angles, and LJ parameters), and an `in.lammps` that goes up to the minimization. If any atom has an environment that is not in the library (similarity below 0.95), it stops without building a force field.
- **Charges and water model** (Calculation settings › Library and details, or `--water`): SOAP gives the type of every atom; its charge is the published value of ClayFF Table 1 for that type (Cygan et al. 2004; 4 decimals, the ions as in the table), with no other charge scheme. Water: the flexible SPC water of ClayFF (Tables 1 and 2), the only model. The charges are not shifted to make the system neutral: the net charge is reported; Table 1 leaves +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al substitution.
- **Adding MD steps**: after the minimization, choose from the lists of LAMMPS fix, compute, and other commands to add them (temperature, pressure, time step, and IDs are filled in automatically). Commands can also be entered by **copy and paste** (select a line and press Cmd+V, or use "Open as text (large, paste)" in the editor to edit the whole input in a large window and press "Apply"). The editor can be enlarged by dragging the border at the bottom or with "Expand".
- **Supercell**: replicates the structure along a, b, and c with OVITO's Replicate, then builds the force field.
- **Surface** (valid only for the (0 0 1) face): cuts the basal plane of a layered mineral from the current calculation (it uses the ClayFF types of the bulk): you give only the face h k l; pymatgen (SlabGenerator) cuts a slab at least 16 Å thick with 30 Å of vacuum and keeps a termination that cuts no cation-O bond. The interlayer ions stay where the cut leaves them, on one face. A face that would cut bonds (an edge) is refused: edges, whose O would need H, are not supported. Then builds the force field.
- **Add water** (packmol, as PACKMOL-GUI): a list of items, each a molecule (the water, or an ion of ClayFF Table 1: Na⁺, K⁺, Cs⁺, Ca²⁺, Ba²⁺, Cl⁻; only what ClayFF types), a number and the region where packmol puts it, given by packmol's constraints: inside/outside box, cube, sphere, ellipsoid, cylinder, over/below plane, xygauss, any number of them at once. The regions are drawn over the structure (side x'–z and top x'–y views, Monte Carlo points outside the van der Waals spheres of the atoms). The number of water molecules can be set from a density and the free volume of the region (Monte Carlo). The empty slices of the structure along z (the vacuum above a surface, an interlayer) are listed and can be used as the box of an item (then the number of water molecules is the density × area × thickness of the slice, as before); at the start every slice at least 6 Å thick is an item of water. A button sets the number of an ion to cancel the net charge (an item of that ion, or a new one in the region of the water). The packmol input can be shown before the run; the run keeps its files in `<structure>_water_packmol/` (pack.inp, the xyz files, packmol.log) and shows the log. Tolerance (default 2 Å, the smallest distance to every other atom) and seed can be set. The coordinates are those of packmol's frame: x' = x − y·cot γ (x' = x when γ = 90°), y, z normal to the a–b plane (standard orientation); a region may run past the top of the cell (all the regions are moved together along c so that none crosses the cell boundary). Then builds the force field (the water and the ions typed by rule, with the chosen water model). The recommended way to make hydrated or interfacial systems: give the dry structure and add the water here. The same from the command line: `soapff-cui add-water STRUCTURE --list`, then e.g. `soapff-cui add-water slab.extxyz --region 30 50 --x 0 0.5 --ion Ca:3 --neutralize Cl -o slab_water.extxyz` (boxes, the ions shared among them in proportion to the water), or the general form `soapff-cui add-water slab.extxyz --spec items.json [--show-input] [--workdir DIR]` with items such as `[{"molecule": "water", "density": 1.0, "constraints": [["inside", "box", [0, 0, 48, 31, 27, 60]]]}, {"molecule": "Na", "number": 4, "constraints": [["inside", "sphere", [15, 13, 54, 5]]]}]`.
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
| `library_nvt_npt/` | NVT and NPT of each of the 11 bulk library structures (soapff-cui `make` inputs, LAMMPS outputs, and analysis; summary table in its README) |
| `montmorillonite_2W/` | Diffusion coefficients of Na montmorillonite with a 2-layer hydrate (1,554 atoms made from the 32-water state of the swelling example; 300, 330, and 366 K; compared with Greathouse et al. 2016) |
| `montmorillonite_swelling/` | How to build the same Na montmorillonite structure as Cygan et al. 2004, and its swelling curve (10 water contents) |
| `muscovite_2M1/` | 0 K elastic constants, tension, compression, and shear of muscovite 2M1 (the property commands of an earlier soapff-gui version, rerun with the current force field) |

All the examples use the current force field: the charges of Table 1 of Cygan et al. 2004 as they are, the LJ parameters unrounded, and the flexible SPC water. The runs were made on an SSH computer (Linux, Open MPI, PPPM 1e-5).

## Supported structures and limits

- Library: Na montmorillonite (dry), brucite, pyrophyllite, muscovite (2M1, 1M), phlogopite, nanpingite, kaolinite, lizardite, and talc. For these crystals, the types and charges of all atoms matched the correct answer, both for supercells of different sizes and for structures after MD at 300 K.
- Water is typed by rule (an O with two H in the shape of a water molecule; ClayFF SPC charges) and left out of the SOAP, in the library and in the structures: the ClayFF type of a mineral atom does not depend on the water around it, so the same library serves dry crystals, hydrated interlayers and mineral-water interfaces. The ions stay in the SOAP (their positions belong to the clay).
- Surfaces: only interfaces that cut no bond, i.e. the basal plane of a layered mineral (cleaved in the interlayer by pymatgen, the interlayer ions left on one face; the atoms keep the ClayFF types of the bulk). Edges (cation-O bonds cut, O protonated by bond valence) are not supported: the GUI and `soapff-cui make-surfaces` refuse those faces, the library has no edge references (no JAEA particle, no edge slabs), and a structure with edges is stopped by the similarity check. There are no edge types or edge charges: every atom has a ClayFF Table 1 type and its charge.
- Other minerals: make their bulk reference (LAMMPS data file with ClayFF type names), `soapff-cui make-surfaces` for the basal slab, and `soapff-cui build-library` with all of them.
- The force field is ClayFF (Tables 1 and 2 of Cygan et al. 2004; ions as in Table 1) with its flexible SPC water.

## Layout

```
src/soapff_gui/  GUI (app.py), command lists (catalog.py), SSH (remote.py), execution (runner.py), library (data/)
src/soapff_cui/  force-field generation (core.py), analysis (analysis.py). The command is soapff-cui
examples/        the examples above
```

## Interfaces in the library

Only interfaces that cut no bond: the basal plane (0 0 1) of the 11 layered minerals (pyrophyllite, pyrophyllite_2, talc, kaolinite, lizardite, brucite, muscovite 1M, muscovite 2M1, phlogopite 1M, nanpingite, montmorillonite (dry, Na)), as slabs in vacuum, with the 11 bulk references: 22 structures, 3,265 library rows. Every other face with |h|, |k|, |l| <= 1 cuts cation-O bonds wherever it is placed and is refused.

| kind | faces | bonds cut | charges |
|---|---|---|---|
| basal (layered minerals) | 001 | none; the interlayer ions on one face | the ClayFF types of the bulk |
| edges | the others | cation-O bonds | not supported |

Which mineral and which face (|h|, |k|, |l| <= 1; each was tested by cutting all 13 faces of the bulk reference):

| Mineral | Interlayer ion | (001) slab and its water interface | Other faces |
|---|---|---|---|
| montmorillonite (dry, Na; Cygan et al. composition) | Na⁺ | yes | refused |
| muscovite 2M1, muscovite 1M | K⁺ | yes | refused |
| phlogopite 1M | K⁺ | yes | refused |
| nanpingite (Cs mica) | Cs⁺ | yes | refused |
| pyrophyllite (two sources), talc | none | yes | refused |
| kaolinite, lizardite | none (1:1 layers) | yes | refused |
| brucite | none | yes | refused |

Of the 143 slabs (11 minerals x 13 faces, cut by pymatgen) only the 11 (001) slabs can be made (formal charge 0, no cation-O bond cut); the other 132 are refused. The 11 (001) slabs were also cut through the GUI Surface dialog from the examples in `library_nvt_npt/`, and cut again from the bulk doubled along a with 30 Å of vacuum: fallback 0 and no wrong type. Water was added to all 11 (122-281 molecules) and the minimization converged (fallback 0). Any other ClayFF mineral can be added with its own bulk reference (see above), and an interface with it follows the same rule.

A mineral-water interface is the basal slab with water added in the empty part ("Add water").

How the slab is cut (valid only for the (0 0 1) face: every other face cuts bonds and is refused, and the library has no references for it): the Surface dialog asks only for the face h k l (of the unit cell; the structure may be a supercell of it, with ions or substitutions). pymatgen's `SlabGenerator` makes the slab of whole planes, at least 16 Å thick, with 30 Å of vacuum along the normal, and lists every termination of that face; `get_slabs(bonds=..., max_broken_bonds=0)` keeps those that break no cation-O bond (Si-O, Al-O, Mg-O, ...) and no O-H of a hydroxyl. If none is left the face is refused. The c axis of the slab cell is then made normal to the slab (`get_orthogonal_c_slab`), so that "Add water" can fill the vacuum. Where the vacuum is can be chosen (GUI: "Vacuum"; command line: `--vacuum-side`): above and below the slab (`both`, the default: the slab in the middle of the cell, as pymatgen centres it), above it only (`above`: the slab at the bottom of the cell, 0.5 Å above c = 0), or below it only (`below`). The cell is periodic along c, so the vacuum above and below is one region either way; the choice matters for the picture and for a non-periodic z (`p p f`). The interlayer ions are on the top face of the slab as pymatgen cuts the bundled minerals; the result line says on which face they are. Nothing is added or moved: no H, no splitting of the interlayer ions (the cut goes through the interlayer, so its ions stay on one face and the other face is bare; the slab has the formula of the bulk and is neutral). The command line does the same: `soapff-cui make-surfaces bulk.data -o out/ [--thickness 16] [--vacuum 16] [--vacuum-side both|above|below]` (every face with |h|, |k|, |l| <= 1 is tried; those that cut bonds are skipped).
