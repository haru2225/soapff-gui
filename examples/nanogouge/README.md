# Nano-gouge: shearing an assemblage of edge-bearing Na-montmorillonite nanoparticles (all-atom, LAMMPS, PBS)

A disordered assemblage of N montmorillonite nanoparticles with edges (the particle of the optional library "Montmorillonite edge",
`../montmorillonite_edge/out_clay_montmorillonite/`, 2,292 atoms, ~6 nm disk, ClayFF + Lammers-type edges, 24 Na⁺) in flexible-SPC water,
built and compacted as in Underwood & Bourg, J. Phys. Chem. C 124, 3702 (2020), then sheared at constant normal stress to get the friction
coefficient of a gouge-like assemblage (grains of random orientation, edges, water, ions), to be compared with the coarse-grained results of
Lin & Hatano (2024; σn = 5.6 and 56 MPa, rate ~1e8 /s) and with the two-plate results of this repository.

| Stage | LAMMPS input | What it does |
|---|---|---|
| build | `build_gouge.py` | N particles, random position and orientation, no overlap; 300 random configurations, the lowest order parameter S is kept (paper: S ≈ 0.1); water at m_w/m_c ≈ 3 (paper) |
| 01 | `lammps/in.01_equil` | minimisation, 100 ps NVT, NPzT (only z moves) at 1 atm |
| 02 | `lammps/in.02_dehydrate` | every `step_ps` delete 2 % of the initial water (1 % below m_w/m_c = 1) and run NPzT: compaction, particles stack (paper: 1 % every 1 ns, 134 ns); stops at `final_mwmc` |
| 03 | `lammps/in.03_consolidate` | hold at the normal stress σn (z only; x, y fixed) |
| 04 | `lammps/in.04_shear` | SLLOD shear (box tilt xz, rate 1e-7 /fs = 1e8 /s) at constant σn (`fix press/langevin` in z); `shear_<σ>.dat` = strain, shear stress, normal stress, μ, Lz, ρb; restart files every 0.5 ns |
| 05 | `lammps/in.05_heal` | stop the shear, hold (static contact time): stress relaxation, re-stacking |

Everything ran on a small test system (3 particles, 27,000 atoms; every stage, the restart of stage 04, `analysis/`) with LAMMPS 7 Feb 2024; the real-size runs have not been done yet.
The force field is the one of the repository (ClayFF, Cygan et al. 2004; the edge types of Lammers et al. as taken whole from `clay_montmorillonite.settings`; flexible SPC water).
The paper used SPC/E with SHAKE and GROMACS; here the water is flexible, so the time step is 1 fs.

## From GitHub to the supercomputer

```
# on the supercomputer (login node)
git clone https://github.com/haru2225/soapff-gui.git
cd soapff-gui/examples/nanogouge

# 1. look at the centre's sample LAMMPS script and copy its settings into pbs/config.env
cat /home/app/lammps/run.sh          # "#PBS" lines -> QUEUE, SELECT, WALLTIME_*;  module line -> MODULES;  mpirun / lmp line -> MPIRUN, LMP
vi pbs/config.env

# 2. build the system (numpy and scipy are needed; or build it on your own computer and copy gouge.* to RUN)
python3 build_gouge.py -n 12 --mwmc 3.0 --nconf 300 --seed 1 -o RUN_N12/gouge

# 3. submit the whole chain: 01 -> 02 -> (03 -> 04 x3 -> 05) for sigma_n = 5.6 and 56 MPa
pbs/submit_chain.sh RUN_N12 "5.6 56"
qstat -u $USER
```

`pbs/submit_chain.sh RUN_N12 "5.6 56" "01 02"` runs only the compaction. Jobs are chained with `qsub -W depend=afterany`; a stage that is already finished is skipped, and
stages 02 and 04 continue from their last data/restart file, so a job stopped by the wall time is simply followed by the next one (`N_DEHYD_JOBS`, `N_SHEAR_JOBS` in `config.env`).
If the centre does not allow `depend`, submit the jobs by hand in the order of the table: `cd RUN_N12; qsub -v STAGE=01 job.sh`, etc. (`job.sh` reads `STAGE`, `SIGMA`).

Size (N = 12, m_w/m_c = 3): about 12 × 2,292 + 3 × 80,000 ≈ 268,000 atoms at the start (the builder takes seconds; S ≈ 0.13 of 300 configurations) and ≈ 40,000 after the dehydration to m_w/m_c = 0.2.
Defaults are small enough for a laboratory-scale first try; for a larger system use `-n 30` and the 24 × 24 × 40 nm box of the paper is approached.
Shear to γ = 3 at 1e8 /s is 30 ns (30 million steps of 1 fs); for a quick look use `-var rate 5.0e-7 -var gmax 1.0` through `EXTRA` (`qsub -v STAGE=04,SIGMA=5.6,EXTRA="-var rate 5.0e-7 -var gmax 1.0" job.sh`).

## Analysis (any computer; numpy and matplotlib)

```
python analysis/order_parameter.py RUN_N12/01_equil.data RUN_N12/02_dehydrate.data RUN_N12/03_consolidate_5.6.data RUN_N12/gouge.json   # S, m_w/m_c, dry bulk density
python analysis/shear_curve.py RUN_N12/shear_5.6.dat RUN_N12/shear_56.dat                                                              # mu(gamma), shear_curve.png
```

`RUN_N12/dehydrate.dat` lists the water number, m_w/m_c, dry bulk density and Lz of every dehydration step; `RUN_N12/shear_<σ>.dump` holds the clay atoms (unwrapped, every 0.2 ns) for the pore-size and stacking analysis.

## Notes and limits

- σn in the input is in MPa (1 MPa = 9.869 atm); the box in x and y is fixed during stages 03-05, z is the free direction (NPzT, as in the paper).
- 6 nm particles are smaller than natural grains and carry many edges per area: read the result as the limit in which edges matter most.
- A real gouge compacted at 5.6 or 56 MPa is far from the dry bulk density reached by a short dehydration; `final_mwmc` (default 0.2) and the water removal speed decide what is sheared, so report them with μ.
- Not done: edge protonation states other than the one of the library, particles of other minerals, Si-O bond breaking.
