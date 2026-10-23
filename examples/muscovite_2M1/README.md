# Example: elastic constants, tension, and shear of muscovite 2M1 (calculation for the paper)

A crystal made by repeating the unit cell of muscovite 2M1 (KAl₂(AlSi₃)O₁₀(OH)₂) 2x2x1 (672 atoms, 20.71 x 17.97 x 20.07 Å, β = 95.76°) was calculated with "+ Properties" of soapff-gui. The coordinate axes are: x is the a axis, y is the b axis (the unique axis of the monoclinic cell), and z is perpendicular to the layers.

## Procedure (GUI)

1. Choose `muscovite_2M1_221.extxyz` with "New calculation" (the bundled library; all atoms are assigned by SOAP, fallback 0, 0 mismatches in the type cross-check).
2. With "+ Properties", select the elastic constants at 0 K (strain 0.005), tension / compression (x y z, both, 10⁹ /s, maximum strain 0.16), and shear (xy xz yz, solid, + and -, 10⁹ /s, maximum strain 0.16). Common settings: 300 K, 1 atm, relaxation NVT 10 ps, NPT 20 ps, NPT direction tri.
3. Press "Start calculation" (SSH computer, 8 MPI processes, 35 minutes). When it finishes, the analysis runs automatically, and the moduli combining both directions (obtained in the same way as Teich-McGoldrick et al.) are also given.

"+ Properties" does not exist in the current version of soapff-gui. The calculation was first made with the version before it was removed (commit 2c43d80). It was redone with the current force field: the force field (`system.data`, `system.settings`, `charges.csv`) was rebuilt by the current soapff-gui from `muscovite_2M1_221.extxyz` (ClayFF Table 1 charges as they are: net -0.0064 e, not shifted to 0; LJ sigma unrounded), and the `in.lammps` of the earlier version was reused unchanged (it holds the "+ Properties" commands). The analysis was made with the analysis code of that version (`analysis/summary.txt`, translated to English). Compared with the earlier run, only the charges (the earlier version shifted every atom by 0.00001 e to make the net charge 0) and the number of digits of the LJ sigma differ.

| Folder | Contents |
|---|---|
| `gui_run_remote/` | The inputs of the procedure above (`in.lammps`, `system.data`, `system.settings`, `charges.csv`, `report.json`), the outputs (`elastic0.dat`, `mech_*.dat`, `log.lammps`), and the GUI analysis (`analysis/summary.txt` and plots). The trajectories (`traj_*.lammpstrj`) and dumps are large and are not included |
| `elastic0_strain_check/` | Checks used to decide the strain size for the 0 K elastic constants (`e005_np4`, `e005_np2`: strain 0.005 with 4 and 2 MPI processes; `e010_np4`: strain 0.01), and `tm_cell`: a calculation with the box fixed to the lattice of Teich-McGoldrick et al. (`change_box` to their box, then the atoms only minimized). The force field (`system.data`, `system.settings`) is the same as in `gui_run_remote/`; each `in.lammps` holds only the 0 K elastic-constant block |

## Results (GPa)

| | This work 0 K | This work 300 K (both directions) | This work 300 K (one direction) | Teich-McGoldrick et al. 0 K | Same, 298 K | Experiment (Vaughan and Guggenheim) |
|---|---|---|---|---|---|---|
| C11 / C22 / C33 | 262 / 233 / 66.4 | | | 236 / 198 / 68 | | 181 / 178 / 59 |
| C44 / C55 / C66 | 8.4 / 16.2 / 95 | | | 15.4 / 14.1 / 94 | | 16.5 / 19.5 / 72 |
| C12 / C13 / C23 | 67 / 23.5 / 21.6 | | | 42 / 25 / 40 | | 49 / 26 / 21 |
| E_x / E_y / E_z | 233 / 212 / 62 | 239 / 225 / 62.4 | 220 / 221 / 53.8 | 222 / 178 / 56 | 218 / 190 / 69 | 149 / 161 / 54 |
| G_yz / G_xz / G_xy | 8.4 / 15.5 / 95 | 14.0 / 19.3 / 88.5 | 8.6 / 13.2 / 92.4 | 15.2 / 12.6 / 94 | 15.2 / 19.6 / 84 | 16.1 / 18.1 / 70 |

The 300 K "both directions" value follows Teich-McGoldrick et al.: the stress-strain data of both directions (tension and compression, + and - shear) are combined, and the slope is obtained by a linear fit for |strain| up to 0.02 in plane, 0.01 perpendicular to the layers, and 0.1 for shear. The differences from Teich-McGoldrick et al. at 298 K are E_x +10, E_y +18, E_z -10, G_xy +6, G_xz -1, and G_yz -8 %. The "one direction" value is the slope up to strain 0.02 for tension and + shear; E_z, G_xz, and G_yz come out small (because the curves for tension perpendicular to the layers and for shear across the layers are curved at the origin).

300 K strength (maximum of the 0.5 ps moving average) and failure strain (the strain at which the stress drops below 70 % of its maximum): tension x 15.1 (0.104), y 13.3 (0.078), z 3.6 (0.106); shear xy 11.6, xz 3.0 (0.142), yz 2.3 (0.133). In compression the stress keeps increasing up to a strain near 0.16 (28-31 GPa). The MD of DelRio et al. (2026) gives about 12 GPa for tension in x and y and about 3 GPa in z, and failure strains of about 9-10 % in x, about 7 % in y, and about 10-11 % in z.

Compared with the earlier run (the same structure and `in.lammps`, the charges shifted to a net 0 and the LJ sigma rounded to 5 digits), the 0 K elastic constants agree within 1.6 GPa, the moduli for both directions within 1.4 %, and the strengths within 0.4 GPa (MD fluctuations): the unshifted Table 1 charges (net -0.0064 e) change nothing that can be seen here.

## Differences from previous work

**When the way the moduli are obtained is matched, the difference from Teich-McGoldrick et al. is within 18 %** (the table above). At first (one direction, up to strain 0.02) the shear across the layers came out 33-43 % smaller, but the main cause was the difference in how the moduli were obtained. For the remaining differences (the in-plane E_x and E_y are 10-18 % larger), the following were considered.

**The difference from experiment (30-50 % stiffer in plane, softer shear across the layers) is a property of ClayFF.** ClayFF is a force field fitted to reproduce crystal structures; the cation-O interaction has no bonded term and consists only of Coulomb and LJ forces. The elastic moduli were not fitted. Teich-McGoldrick et al. also report that it "consistently overestimates the in-plane elastic response". Between the layers there is only the Coulomb force between K⁺ and the layer O atoms, so the shear across the layers (C44, G_yz) comes out weak. The experimental values are from Brillouin scattering (adiabatic, natural samples containing a little Fe and Mg), so the conditions also differ somewhat from the MD (isothermal, ideal composition).

**Differences from Teich-McGoldrick et al., who also used ClayFF** (the 0 K C11, C22, and C12 are 25-35 GPa larger, C23 is 18 GPa smaller, and the 300 K G_xz and G_yz from one direction are 33-43 % smaller) come from the following differences.

1. **Reference lattice (checked and ruled out)**: the zero-stress box here has a, b, c of 5.151, 8.955, and 19.789 Å, which are -1.0, -0.8, and -1.4 % smaller than the optimized box of Teich-McGoldrick et al. (5.202, 9.024, and 20.079 Å, with β kept at the input value of 95.7566°). When the box was fixed to their lattice, only the atoms were minimized, and the constants were obtained by the same method (strain +/-0.005) (`elastic0_strain_check/tm_cell/`), that box was under pressures of -3.5, -2.9, and -1.3 GPa in x, y, and z (in the tensile direction) with our force field; the in-plane C11, C22, and C66 were 269, 248, and 88 GPa, almost unchanged, and C33, C44, and C55 dropped to 45.9, 2.7, and 10.7 GPa because the layers were separated. Therefore the in-plane stiffness being higher than theirs is not due to the reference lattice.
2. **Arrangement of tetrahedral Al and positions of K⁺**: they used a model built from the Kuwahara structure with Al placed according to Loewenstein's rule (equal numbers of Si₄Al₂ and Si₅Al₁ rings); here the arrangement is that of the structure provided by a co-author. The O next to a tetrahedral Al (obts, -1.1688 e) has a larger charge than ob (-1.05 e), and the interlayer K⁺ is attracted near it, so the arrangement of Al changes the interaction across the layers (C44, C55, C23).
3. **How the 0 K elastic constants are obtained**: they used a strain of 5x10⁻⁷, and C44 has a large uncertainty, 15.4 ± 12.2 GPa. Here, the minimization stops leaving a stress of up to about 0.1 GPa because of the PPPM forces (in this run -0.11 GPa along z at the zero-stress box), so a very small strain gives scattered values; strain 0.005 (2 and 4 MPI processes) and 0.01 agree within 3 GPa, and 0.005 was chosen (`elastic0_strain_check/`). In the earlier version, strain 0.001 gave C44 of 16-19 GPa, close to theirs, but scattered by 5-10 GPa (that check was not redone).
4. **System size and calculation conditions**: 672 and 5,376 atoms, PPPM (10⁻⁵) and Ewald (10⁻⁴), time steps of 0.5 fs and 1 fs, and strain rates of 10⁹ and 10⁸ /s. The effect on the moduli is small, but the strength comes out higher at faster rates.

The difference in strength from DelRio et al. (ours is 10-20 % higher) is thought to be mainly due to the differences in strain rate and system size (their MD conditions are in the supplementary material and have not been checked here). The failure strains and the relation that the in-plane strength is about 4 times the strength perpendicular to the layers are the same.
