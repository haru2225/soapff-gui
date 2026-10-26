# soapff-GUI

A GUI for clay minerals. You pick a structure file (elements, coordinates, and box only); SOAP assigns the ClayFF type and charge of every atom, the GUI writes the LAMMPS input, runs the MD, and analyzes it. The force field is ClayFF (as tabulated in Cygan et al. 2004) with its flexible SPC water. Every calculation of the paper is in `examples/` with its inputs, a script that redoes it, and its results (see "Examples and the paper").

## Features

- **Automatic force-field generation**: each atom gets the type and charge of the most similar atom in the bundled library (11 bulk layered minerals including dry Na montmorillonite, and the basal (0 0 1) slab in vacuum of each: 22 structures, 3,265 rows after pruning duplicates; no edges). An optional library, a montmorillonite particle with edges and its own force field, can be used instead ("Optional library: montmorillonite edges" below); it is chosen in Calculation settings › Library and details › Bundled library, or in the Structure library window. The GUI writes `system.data`, `system.settings` (bonds, angles, and LJ parameters), and an `in.lammps` that goes up to the minimization. If any atom has an environment that is not in the library (similarity below 0.95), it stops without building a force field.
- **Charges and water model** (Calculation settings › Library and details, or `--water`): SOAP gives the type of every atom; its charge is the published value of ClayFF Table 1 for that type (Cygan et al. 2004; 4 decimals, the ions as in the table), with no other charge scheme. Water: the flexible SPC water of ClayFF (Tables 1 and 2), the only model. The charges are not shifted to make the system neutral: the net charge is reported; Table 1 leaves +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al substitution.
- **Adding MD steps**: after the minimization, choose from the lists of LAMMPS fix, compute, and other commands to add them (temperature, pressure, time step, and IDs are filled in automatically). Commands can also be entered by **copy and paste** (select a line and press Cmd+V, or use "Open as text (large, paste)" in the editor to edit the whole input in a large window and press "Apply"). The editor can be enlarged by dragging the border at the bottom or with "Expand".
- **Supercell**: replicates the structure along a, b, and c with OVITO's Replicate, then builds the force field.
- **Surface** (valid only for the (0 0 1) face): cuts the basal plane of a layered mineral from the current calculation with pymatgen, with a chosen vacuum thickness ("Surface" below). A face that would cut bonds (an edge) is refused. Then builds the force field.
- **Add water** (packmol, as PACKMOL-GUI): water and ions in regions given by packmol's constraints (box, sphere, cylinder, plane, ...), drawn over the structure ("Add water" below). Then builds the force field.
- **Swelling curve**: runs NPT while adding water molecules to the interlayer, and records the box of each state in `swelling.dat` (same procedure as Cygan et al. 2004; the plot is made by `examples/montmorillonite_swelling/analyze_swelling.py`).
- **Running**: LAMMPS runs on this computer or on another computer registered over SSH (just enter `user@host` to register it).
- **Analysis**: diffusion coefficients (when a `compute msd` is present), and SOVA PDF, RDF, coordination numbers, and S(Q). The SOVA dr, maximum r, Q range and step, X-ray or neutron, number of frames, and the pairs to compute, including pairs split by type (example: K-ob K-obts), can be set on the analysis screen (Settings > Analysis). Structures can be displayed in OVITO ("OVITO" (or the menu item "Open in OVITO")).

## Installation

| OS | Steps |
|---|---|
| macOS | `git clone https://github.com/haru2225/soapff-gui.git && sh soapff-gui/install.sh`, then launch from `~/Applications/soapff-GUI.app` |
| Linux | `sudo apt install python3-venv python3-tk` (Ubuntu example), then the same as macOS. Launch with `~/.soapff-gui/venv/bin/soapff-gui` |
| Windows | Install Python 3.10-3.12 from python.org, run `git clone ...`, then `powershell -ExecutionPolicy Bypass -File soapff-gui\install.ps1`, and launch from the desktop shortcut (not tested. Running over SSH is not available, so to use SSH, install in WSL2 following the Linux steps) |

- **`install.sh` installs everything that is needed automatically** (it does not use `sudo`): the Python libraries (numpy, scipy, ASE, dscribe, matplotlib, ovito, pymatgen), **LAMMPS and Open MPI** (the MPI build from conda-forge, in `~/.soapff-gui/lammps`, about 1 GB; skipped if `lmp` is already installed), **packmol**, **SOVA** (`pip install sovapy`, Python 3.11 or later), and the **OVITO application** (on macOS `brew install --cask ovito`; on Linux from conda.ovito.org into `~/.soapff-gui/ovito`, about 1.5 GB). `SOAPFF_NO_LAMMPS=1` skips LAMMPS and packmol, `SOAPFF_NO_OVITO=1` skips OVITO.
- The Windows `install.ps1` installs only the Python libraries and SOVA (get LAMMPS from https://packages.lammps.org/windows.html and OVITO from https://www.ovito.org). If they are not found, specify their locations in "Settings > Run on".
- If no Python with Tk is found, run `brew install python-tk@3.12` or `brew install uv` and then run `install.sh` again. To update, run `git pull` (if dependencies were added, run `~/.soapff-gui/venv/bin/python -m pip install -e soapff-gui`).

## Usage

1. Choose a structure file with "New calculation" (or "Supercell" to make a supercell). The force field and an input that goes up to the minimization are created.
2. Add MD steps with "+ fix", "+ compute", and "+ Other" (example: soapff relaxation, then NVT, then NPT, then write the last structure). For the swelling curve, use "+ other > soapff procedures > swelling curve".
3. Choose where to run ("Run on") and press "Start calculation". When the run finishes, the GUI analyzes it and shows the plots and numbers (`analysis/summary.txt`).

The same force field from the command line: `soapff-cui make STRUCTURE -o OUT` (and `soapff-cui --help` for the other commands: `add-water`, `make-surfaces`, `build-library`, `check-library`, `analyze`, `open`).

### Running on another computer (SSH)

Force-field generation and analysis run on this computer; only LAMMPS runs on the remote side (Linux or macOS).

1. **Remote side**: you can log in over SSH, and `lmp` (with KSPACE and MOLECULE) and `mpirun` are available.
2. **Register**: next to "Run on" in the "Run" section, press **+ Add**, enter the `ssh` destination (example: `taro@192.168.0.10`, or `taro@192.168.0.10:10022` if the port is not 22), and press "Add". When asked for the password, enter it; the GUI connects, checks the remote CPU and LAMMPS, and sets the number of MPI processes (physical cores) and `mpirun` automatically.
3. **Run**: choose that computer under "Run on" and press "Start calculation". The GUI sends the input and starts it with `nohup` (it keeps running if the connection drops); when it finishes, the whole result folder is copied back.
4. **Remove**: choose it under "Run on" and press **- Remove**. Detailed settings (the remote folder, the locations of `lmp` and `mpirun`, commands to run first such as `module load lammps`, and the number of processes) can be changed in "Settings > Run on".

Passwords are not saved (you are asked each time you register, run, or press "Test connection"). If you can log in with a key (`ssh-copy-id`), leave the password empty.

## Force field and libraries

### Default library: ClayFF (Cygan et al. 2004)

`src/soapff_gui/data/library_clayff_bulk.npz`, rebuilt and checked by `examples/library/run.py`.

- References: 11 bulk layered minerals (Na montmorillonite (dry, Cygan et al. composition), brucite, pyrophyllite (two sources), muscovite 2M1 and 1M, phlogopite 1M, nanpingite, kaolinite, lizardite, talc) and the basal (0 0 1) slab in vacuum of each: 22 structures, 3,265 rows (the slab rows pruned at similarity 0.999). For every row the nearest other row of the same element carries the same charge (3,265 / 3,265).
- Water is typed by rule (an O with two H in the shape of a water molecule; ClayFF SPC charges) and left out of the SOAP, in the library and in the structures: the ClayFF type of a mineral atom does not depend on the water around it, so the same library serves dry crystals, hydrated interlayers and mineral-water interfaces. The ions stay in the SOAP (their positions belong to the clay). Ions in solution that the library does not describe (e.g. Cl⁻) get their formal charge and ion type by rule.
- The force field is ClayFF (Tables 1 and 2 of Cygan et al. 2004; ions as in Table 1) with its flexible SPC water; LJ parameters are written unrounded (ε = D₀, σ = R₀ / 2^(1/6)).

### Optional library: montmorillonite edges

`src/soapff_gui/data/library_clayff_edge.npz` (Calculation settings › Library and details › Bundled library: "Montmorillonite edge (Underwood & Bourg 2020)"; or the Structure library window). The default library ("ClayFF (Cygan et al. 2004)") is unchanged and stays the default. Rebuilt by `examples/library/run.py`, checked by `examples/montmorillonite_edge/run.py`.

- References: only the montmorillonite particle with edges on every side (`clay_montmorillonite.data`, a single layer, 2,292 atoms; `clay_montmorillonite_2.data`, the same particle twice along c, adds no new environment): 2,292 rows. No ClayFF structure: other minerals and bulk structures stop (not in this library).
- Force field: taken whole from `clay_montmorillonite.data` and `clay_montmorillonite.settings` and stored in the library (`core.read_forcefield`, the `forcefield` key of the npz): the charge of every type, the LJ of every type (pair_coeff), the hydroxyl O-H bond (553.935 kcal/mol/Å², 1.0 Å), the M-O-H angle (30.0191 kcal/mol/rad², 109.47°) and the edge metal-O-octahedral metal angle (800 kcal/mol/rad², 100°). This holds also for the types shared with ClayFF, whose values differ slightly from Table 1 (charges ohs -1.081, obos -1.181, mgo +1.361; σ rounded to 4 decimals; Na σ 2.6378 Å; O-H bond 553.935 instead of 554.1349; M-O-H 30.0191 instead of 30.0). Edge types: oed, oab, omb (edge OH), oaa, oma (edge O bound to two metals), alb, ala (edge Al), mgb, mga (edge Mg).
- Not from those files: the masses (the element masses: the files give 13.0 and 12.0, the atomic numbers, for the edge Al and Mg), the pair style, cutoff and long-range Coulomb (the run settings of the GUI), and the water and ions added with "Add water" that the model does not have (flexible SPC water, Cl⁻, ...: ClayFF).
- Source: T. R. Underwood, I. C. Bourg, J. Phys. Chem. C 124, 3702 (2020) (after Lammers et al.).
- Checks: leave-one-out 2,292 of 2,292 rows. The particle and the particle twice (4,584 atoms), from elements and coordinates only: no wrong type, no wrong charge, no fallback, smallest similarity 1.000; every pair_coeff, bond and angle coefficient equal to the files; 288 bonds and 624 angles as in the file. A LAMMPS single point (Ewald 1e-5, 10 Å) of the written system.data / system.settings equals that of the original files in every term (pe -509528.707277 kcal/mol, the same to 1e-6). When this option was added, the output of the default library was compared before and after the change: unchanged byte for byte, with and without water and ions.

### Other minerals

Make their bulk reference (LAMMPS data file with ClayFF type names in the Masses comments), `soapff-cui make-surfaces` for the basal slab, and `soapff-cui build-library` with all of them (as `examples/library/run.py` does for the bundled library).

## Surface

How the slab is cut (valid only for the (0 0 1) face: every other face cuts bonds and is refused, and the library has no references for it): the Surface dialog asks only for the face h k l (of the unit cell; the structure may be a supercell of it, with ions or substitutions) and the vacuum. It uses the ClayFF types of the bulk of the current calculation. pymatgen's `SlabGenerator` makes the slab of whole planes, at least 16 Å thick, with the chosen vacuum along the normal (GUI: "Vacuum thickness (Å)", 5-1000 Å, default 30 Å; the file name gets `_vac<thickness>A` when it is not 30; command line: `--vacuum`). pymatgen rounds the cell up to whole layers, so the vacuum made can be larger and is reported: e.g. for dry Na montmorillonite 175 Å gives c = 195.1 Å and 183 Å gives c = 204.4 Å, room for a droplet of 80 Å diameter. `get_slabs(bonds=..., max_broken_bonds=0)` keeps the terminations that break no cation-O bond (Si-O, Al-O, Mg-O, ...) and no O-H of a hydroxyl; if none is left the face is refused. The c axis of the slab cell is then made normal to the slab (`get_orthogonal_c_slab`), so that "Add water" can fill the vacuum.

Where the vacuum is can be chosen (GUI: "Vacuum"; command line: `--vacuum-side`): above and below the slab (`both`, the default: the slab in the middle of the cell, as pymatgen centres it), above it only (`above`: the slab at the bottom of the cell, 0.5 Å above c = 0), or below it only (`below`). The cell is periodic along c, so the vacuum above and below is one region either way; the choice matters for the picture and for a non-periodic z (`p p f`). The interlayer ions stay on one face of the slab (the top face as pymatgen cuts most of the bundled minerals; the bottom one for nanpingite); the result line says on which face they are. Nothing is added or moved: no H, no splitting of the interlayer ions (the cut goes through the interlayer, so its ions stay on one face and the other face is bare; the slab has the formula of the bulk and is neutral). The command line does the same: `soapff-cui make-surfaces bulk.data -o out/ [--thickness 16] [--vacuum 16] [--vacuum-side both|above|below]` (every face with |h|, |k|, |l| <= 1 is tried; those that cut bonds are skipped).

Only interfaces that cut no bond are made: the basal plane (0 0 1) of the 11 layered minerals. Every other face with |h|, |k|, |l| <= 1 cuts cation-O bonds wherever it is placed and is refused (of the 143 slabs, 11 minerals x 13 faces, only the 11 (001) slabs can be made; `examples/surfaces_001/`).

| kind | faces | bonds cut | charges |
|---|---|---|---|
| basal (layered minerals) | 001 | none; the interlayer ions on one face | the ClayFF types of the bulk |
| edges | the others | cation-O bonds | not supported (cutting refused; a montmorillonite particle with edges: the optional edge library) |

| Mineral | Interlayer ion | (001) slab and its water interface | Other faces |
|---|---|---|---|
| montmorillonite (dry, Na; Cygan et al. composition) | Na⁺ | yes | refused |
| muscovite 2M1, muscovite 1M | K⁺ | yes | refused |
| phlogopite 1M | K⁺ | yes | refused |
| nanpingite (Cs mica) | Cs⁺ | yes | refused |
| pyrophyllite (two sources), talc | none | yes | refused |
| kaolinite, lizardite | none (1:1 layers) | yes | refused |
| brucite | none | yes | refused |

## Add water

A list of items, each a molecule (the water, or an ion of ClayFF Table 1: Na⁺, K⁺, Cs⁺, Ca²⁺, Ba²⁺, Cl⁻; only what ClayFF types), a number and the region where packmol puts it, given by packmol's constraints: inside/outside box, cube, sphere, ellipsoid, cylinder, over/below plane, xygauss, any number of them at once (e.g. a sphere and a plane for a hemispherical droplet). The regions are drawn over the structure (side x'–z and top x'–y views, Monte Carlo points outside the van der Waals spheres of the atoms). The number of water molecules can be set from a density and the free volume of the region (Monte Carlo); a number that does not fit (above 1.2 g/cm³ of the free volume) is refused before packmol runs. The empty slices of the structure along z (the vacuum above a surface, an interlayer) are listed and can be used as the box of an item (then the number of water molecules is the density × area × thickness of the slice); at the start every slice at least 6 Å thick is an item of water. A button sets the number of an ion to cancel the net charge (an item of that ion, or a new one in the region of the water). The packmol input can be shown before the run; the run keeps its files in `<structure>_water_packmol/` (pack.inp, the xyz files, packmol.log), shows the log, and can be stopped. Tolerance (default 2 Å, the smallest distance to every other atom) and seed can be set. The coordinates are those of packmol's frame: x' = x − y·cot γ (x' = x when γ = 90°), y, z normal to the a–b plane (standard orientation); a region may run past the top of the cell (all the regions are moved together along c so that none crosses the cell boundary). Then builds the force field (the water and the ions typed by rule, with the chosen water model). The recommended way to make hydrated or interfacial systems: give the dry structure and add the water here.

A droplet on the basal plane is made in the GUI alone: Surface with the vacuum above the slab only and a vacuum thicker than the droplet, Supercell, then Add water with a sphere (or ellipsoid, cylinder, ...) constraint (`examples/meng2022_droplet/`).

The same from the command line: `soapff-cui add-water STRUCTURE --list`, then e.g. `soapff-cui add-water slab.extxyz --region 30 50 --x 0 0.5 --ion Ca:3 --neutralize Cl -o slab_water.extxyz` (boxes, the ions shared among them in proportion to the water), or the general form `soapff-cui add-water slab.extxyz --spec items.json [--show-input] [--workdir DIR]` with items such as `[{"molecule": "water", "density": 1.0, "constraints": [["inside", "box", [0, 0, 48, 31, 27, 60]]]}, {"molecule": "Na", "number": 4, "constraints": [["inside", "sphere", [15, 13, 54, 5]]]}]`.

## Examples and the paper

Every example has a `README.md`, its inputs, a `run.py` that redoes the calculation with the Python of soapff-gui (`~/.soapff-gui/venv/bin/python run.py`), and its results (`results.json`). The table follows the order of the paper.

| Paper | Folder | What it holds |
|---|---|---|
| Fig. 1 (workflow), Table "larger crystals", "MD with LAMMPS" | `supercells/` | brucite 5×5×4 and pyrophyllite 3×2×3 against the bond-graph ClayFF types; their MD (`md.py`) |
| References of the library, Fig. "minerals", Table "one-to-one" | `library/` | the reference structures (11 bulk crystals, 11 basal slabs, the 2-layer hydrate, the particle with edges); the default, edge and initial 12-structure libraries rebuilt and checked (l_max 3, 4, 5) |
| Table "sheet silicates", "MD with LAMMPS" | `sheet_silicates/` | 8 sheet silicates: supercells and the last frame of a 300 K NVT/NPT MD against the truth; lattice parameters, temperatures, d001 |
| "MD procedure" | `library_nvt_npt/` | NVT and NPT of the 11 bulk library structures with analysis |
| Basal planes, tests (2) and (3) | `surfaces_001/` | 13 faces of the 11 minerals cut with pymatgen (only (001) made); the (001) slabs cut again from the bulk doubled along a and typed |
| Test (4), Fig. "interfaces" | `interfaces_001_water/` | water on the (001) slabs of the 10 minerals of the figure, force field and minimization |
| Test (5), Fig. "shapes" | `interface_shapes/` | 6 initial structures made with packmol's constraints (sphere, hemisphere, cylinder, cavity, xygauss, NaCl film), force field and minimization |
| Test (6), Fig. "Meng" | `meng2022_droplet/` | the NaCl droplet of Meng et al. 2022 on montmorillonite (001), 105,977 atoms; as in the paper and by the GUI's steps only |
| Test (7) | `montmorillonite_edge/` | the optional edge library on the particle: types, charges, coefficients, single point against the original files, minimizations |
| Swelling curve, Fig. "swelling" | `montmorillonite_swelling/` | how to build the Na montmorillonite of Cygan et al. 2004 and its swelling curve (10 water contents) |

Other examples:

| Folder | Contents |
|---|---|
| `montmorillonite_2W/` | Diffusion coefficients of Na montmorillonite with a 2-layer hydrate (1,554 atoms made from the 32-water state of the swelling example; 300, 330, and 366 K; compared with Greathouse et al. 2016) |
| `muscovite_2M1/` | 0 K elastic constants, tension, compression, and shear of muscovite 2M1 (the property commands of an earlier soapff-gui version, rerun with the current force field) |

All the examples use the current force field: the charges of Table 1 of Cygan et al. 2004 as they are, the LJ parameters unrounded, and the flexible SPC water (the edge example: the particle's own force field). The MD of `library_nvt_npt/`, `montmorillonite_2W/`, `montmorillonite_swelling/` and `muscovite_2M1/` was run on an SSH computer (Linux, Open MPI, PPPM 1e-5); the minimizations of the interface examples use Ewald 1e-4 (PPPM stopped at the start of a minimization with the LAMMPS on the machine used for the paper): their input is `in.min` beside the GUI's `in.lammps`, `system.settings` is kept as the GUI writes it. Each example folder holds the input structure (`structure.extxyz`, or `<name>.extxyz` / `initial.xyz` in the older examples), the files soapff-gui writes for it, and the LAMMPS input of every calculation run for the test (`in.lammps`, `in.min`, `in.original` / `in.written`) with its log; trajectories (`*.lammpstrj`) and packmol's xyz files are not kept.

## Supported structures and limits

- Library: Na montmorillonite (dry), brucite, pyrophyllite, muscovite (2M1, 1M), phlogopite, nanpingite, kaolinite, lizardite, and talc. For these crystals, the types and charges of all atoms matched the correct answer, both for supercells of different sizes and for structures after MD at 300 K (`supercells/`, `sheet_silicates/`).
- Surfaces: only interfaces that cut no bond, i.e. the basal plane of a layered mineral (cleaved in the interlayer by pymatgen, the interlayer ions left on one face; the atoms keep the ClayFF types of the bulk). Edges (cation-O bonds cut, O protonated by bond valence) cannot be cut: the GUI and `soapff-cui make-surfaces` refuse those faces, the default library has no edge references, and a structure with edges is stopped by the similarity check. With the default library there are no edge types or edge charges: every atom has a ClayFF Table 1 type and its charge. The optional edge library types a montmorillonite particle with edges, with the parameters of its own model.
- Not tested, and not in the library: other faces and edges (with the default library), corners where two faces meet, other interlayer cations (Ca²⁺, ...) and numbers of layers, framework silicates (every face cuts Si-O bonds), minerals with elements ClayFF has no type for (Ti, S, ...) or that are not neutral with ClayFF (e.g. Fe²⁺ in annite).

## Layout

```
src/soapff_gui/  GUI (app.py), command lists (catalog.py), SSH (remote.py), execution (runner.py), libraries (data/)
src/soapff_cui/  force-field generation (core.py), ClayFF parameters (clayff.py), surfaces (surfaces.py), water and packmol
                 (water.py), analysis (analysis.py). The command is soapff-cui
examples/        the examples above
```
