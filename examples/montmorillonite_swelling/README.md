# Example: swelling curve of Na montmorillonite (same calculation as Cygan et al. 2004)

For the Na montmorillonite Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆·nH₂O of Cygan et al. (2004, J. Phys. Chem. B 108, 1255) (4 unit cells, one layer), NPT was run while water was added to the interlayer, and the relation between the basal spacing d001 and the amount of water was obtained.

## 1. Structure (`build_structure.py`)

A pyrophyllite layer (the library reference) is replicated 2x2x1, 1 tetrahedral Si is changed to Al, 2 octahedral Al atoms are changed to Mg, and 3 Na⁺ ions are placed in the interlayer. The Mg positions are the ones with the lowest static ClayFF energy (with the box also relaxed) among all 120 ways of choosing 2 of the 16 octahedral sites (4 placements of Na⁺ for each; the 29 arrangements in which an O bonds to both a tetrahedral Al and a Mg (the ClayFF obss type) are excluded because the library does not have them). All results are written to `screening.csv`, and the chosen structure (163 atoms, dry) to `dry.extxyz`.

Na⁺ in a dry interlayer is an environment that no other reference has, so this structure (with ClayFF types and standard charges, `dry_reference.data`) has been added to the library.

## 2. Calculation (soapff-gui)

1. Choose `dry.extxyz` with "New calculation" (fallback 0, 0 mismatches in the type cross-check).
2. With "+ other > soapff procedures > swelling curve", specify the numbers of water molecules as "4 8 12 16 20 25 32 40 48 59" (each state is equilibrated with NPT for 30 ps, and the following 20 ps are recorded). The structure has no water type, so a water type is added to `system.data` and `water.mol` is written.
3. Press "Start calculation" (SSH computer, 4 MPI processes, about 15 minutes).

For each state, water molecules are inserted into the interlayer (the region between the clay and the box edge, on the side across the periodic boundary) with `create_atoms ... mol` (at least 2 Å from other atoms; if they do not fit, the box is stretched by 0.5 Å in z), the clay is moved back to the center of the box, and then a minimization, a 2 ps relaxation with damping, and NPT (300 K, 1 atm) are run. The NPT moves x, y, z, and the xy tilt, and keeps the tilts across the layer, xz and yz, fixed (when water enters, the layer slides and the tilt becomes large, which makes LAMMPS stop; d001 = lz does not depend on this).

`gui_run/` contains the inputs (`in.lammps`, `system.data`, `system.settings`, `water.mol`) and the outputs (`log.lammps`, `swelling.dat`, and the last structure of each state `swelling_states.lammpstrj` (large, so it is not in git)).

## 3. Analysis (`analyze_swelling.py`)

From `swelling.dat`, the mean and standard deviation of lz (= d001) over the last 20 ps of each state, and the water content (g H₂O / g clay; the clay is all atoms other than water, 2,945 g/mol), are obtained (`gui_run/swelling.csv`, `swelling.png`).

| Water (molecules) | g/g | d001 (Å) |
|---|---|---|
| 4 | 0.024 | 11.35 |
| 8 | 0.049 | 12.03 |
| 12 | 0.073 | 12.14 |
| 16 | 0.098 | 12.27 |
| 20 | 0.122 | 12.43 |
| 25 | 0.153 | 14.17 |
| 32 | 0.196 | 14.92 |
| 40 | 0.245 | 15.95 |
| 48 | 0.294 | 17.34 |
| 59 | 0.361 | 18.92 |

The 1-layer hydrate (about 12.0-12.4 Å at 0.05-0.12 g/g), the expansion to the 2-layer hydrate (about 14-15 Å at 0.15-0.2 g/g), and the continued swelling with more water were the same as in the MD of Fig. 8 of Cygan et al.
