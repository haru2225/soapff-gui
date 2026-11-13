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
| 04 | `lammps/in.04_shear` | SLLOD shear (box tilt xz, rate 1e-7 /fs = 1e8 /s) at constant σn (`fix press/langevin` in z); `shear_<σ>.dat` = strain, shear stress, normal stress, μ, Lz, ρb; restart files every 10 ps and at each job stop |
| 05 | `lammps/in.05_heal` | stop the shear, hold (static contact time): stress relaxation, re-stacking |

Everything ran on a small test system (3 particles, 27,000 atoms; every stage, the restart of stage 04, `analysis/`) with LAMMPS 7 Feb 2024; the real-size runs have not been done yet.
The force field is the one of the repository (ClayFF, Cygan et al. 2004; the edge types of Lammers et al. as taken whole from `clay_montmorillonite.settings`; flexible SPC water).
The paper used SPC/E with SHAKE and GROMACS; here the water is flexible, so the time step is 1 fs.

## On the supercomputer (no git, no soapff-gui, nothing to edit: one self-contained pack of 18 files)

```
curl -L https://github.com/haru2225/soapff-gui/raw/main/examples/nanogouge/nanogouge_pack.tar.gz | tar xz     # or scp nanogouge_pack.tar.gz and tar xzf it
cd nanogouge_pack
./go.sh                          # reads /home/app/lammps/run.sh, builds the system (N=12) and submits the whole chain
qstat -u $USER
```

`go.sh` needs no settings: `autoconfig.sh` reads the centre's sample script (`/home/app/lammps/run.sh`, or `CONFIG_SRC=/path/run.sh ./go.sh`) and takes from it every `#PBS` line (queue, resources, wall time), the `module` lines,
the Singularity image and the MPI launcher, and writes `RUN_N12/config.auto` (it prints what it found). `NODES=8 ./go.sh` rescales `select=<n>`. `pbs/config.env` holds optional overrides only (all commented out).
numpy and scipy are installed with `pip --user` if missing. `./go.sh 12 "5.6 56" "01 02"` = N particles, normal stresses in MPa, stages (here only the compaction).
`./go.sh 67 "5.6 56" all 4 RUN_N67_b` (5th argument = run folder) runs the same calculation again in another folder: the initial structure, `config.env` and `lmp.auto`
of `RUN_N67` are copied. `./go.sh 67 "5.6 56" all 4 RUN_N67_s2 2` (6th argument = seed, default 1) is an independent run: another random initial structure,
and LAMMPS uses the seed too (velocities, water removal, thermostat; `RUN_N67_s2/lammps_seed`). Each folder has its own chain and can run at the same time as the others.
Wall time of the jobs: `env WALLTIME=3 ./go.sh 67 "5.6 56" all 4` (hours, or `3:00:00`) replaces every `walltime` of the qsub options for this folder and keeps it in
`RUN_N67/walltime` (later jobs and `./go.sh` runs use it). The waiting jobs of this folder (old wall time) are deleted and new ones queued after the job working now;
the running job goes on. The queue must allow that wall time: `qstat -Qf <queue> | grep walltime`.
`./pbs/jobs.sh` lists your PBS jobs with the run folder each one works in; `./pbs/jobs.sh keep RUN_N67` deletes all your other jobs (it asks first).
The pack is made by `make_pack.sh` from this folder (the particle is `particle.data` / `particle.settings`, a copy of `../montmorillonite_edge/out_clay_montmorillonite/system.*`);
from a clone of the repository the same steps are `python3 build_gouge.py -n 12 -o RUN_N12/gouge` and `pbs/submit_chain.sh RUN_N12 "5.6 56"`.

Singularity (detected from the sample script, or set `SIF="/path/to/lammps.sif"` in `config.env`): every job then runs `mpirun -np <all cores> --hostfile ... singularity exec --bind $PWD <image> lmp ...`.
All cores of every allocated node are used (the number of MPI ranks = nodes × `nproc`, whatever `mpiprocs` the `#PBS` line says; `CORES_PER_NODE` overrides it; `HOSTFILE_FMT="mpich"` for Intel MPI / MPICH).
To run one stage by hand like the centre's sample: edit the marked lines of `run.sh` (queue, `select`, `STAGE`, `SIGMA`) and `qsub run.sh` in a run folder made by `./go.sh`.

`./go.sh 12 "5.6 56" "01 02"` = N particles, normal stresses in MPa, stages (here only the compaction). The pack is made by `make_pack.sh` from this folder
(the particle is `particle.data` / `particle.settings`, a copy of `../montmorillonite_edge/out_clay_montmorillonite/system.*`); from a clone of the repository the same steps are `python3 build_gouge.py -n 12 -o RUN_N12/gouge` and `pbs/submit_chain.sh RUN_N12 "5.6 56"`.

`pbs/submit_chain.sh RUN_N12 "5.6 56" "01 02"` runs only the compaction. Any wall time works, 1 h too: `submit_chain.sh` queues `N_JOBS` (10) identical jobs chained with
`qsub -W depend=afterany`; each job works on the first unfinished stage, continues it from its last restart, stops LAMMPS cleanly `TIME_MARGIN` (360 s) before the
end of its wall time (`timer timeout`, restart `*.restart.stop`), goes on with the next stage when one is finished, and the last job of the chain submits the next one
while work is left (if the compute nodes cannot `qsub`, run `./go.sh` again; it continues). A LAMMPS error stops the chain (`chain.failed`, the following jobs end at once).
`./status.sh RUN_N12` shows the speed and the computing time left of the current stage. What the jobs print goes at once to `RUN_N12/jobs.log` (the PBS
output files `ng_*.o*` appear only when a job has ended). The first job finds a LAMMPS with `fix press/langevin` that starts on the node (the centre's old
`lmp_icex` 3Mar20 has none: the notes about it in the log are normal) and checks that `mpirun` spreads the MPI ranks over all nodes (HPE MPT: else
`mpiexec_mpt`); both are kept in `RUN_N12/lmp.auto` for the next jobs (delete it to search again). Each stop at the wall time also writes the current
structure to `<stage>.latest.data` (e.g. `01_equil.latest.data`; opens in OVITO).
If the centre does not allow `depend`, submit the jobs by hand in the order of the table: `cd RUN_N12; qsub -v STAGE=01 job.sh`, etc. (`job.sh` reads `STAGE`, `SIGMA`).

Size (N = 12, m_w/m_c = 3): about 12 × 2,292 + 3 × 80,000 ≈ 268,000 atoms at the start (the builder takes seconds; S ≈ 0.13 of 300 configurations) and ≈ 40,000 after the dehydration to m_w/m_c = 0.2.
Defaults are small enough for a laboratory-scale first try; for a larger system use `-n 30` and the 24 × 24 × 40 nm box of the paper is approached.
Shear to γ = 3 at 1e8 /s is 30 ns (30 million steps of 1 fs); for a quick look use `-var rate 5.0e-7 -var gmax 1.0` through `EXTRA` (`qsub -v STAGE=04,SIGMA=5.6,EXTRA="-var rate 5.0e-7 -var gmax 1.0" job.sh`).

## Analysis (any computer; numpy and matplotlib; `analysis/` in the repository, next to `go.sh` in the pack)

```
python3 order_parameter.py RUN_N12/01_equil.data RUN_N12/02_dehydrate.data RUN_N12/03_consolidate_5.6.data RUN_N12/gouge.json   # S, m_w/m_c, dry bulk density
python3 shear_curve.py RUN_N12/shear_5.6.dat RUN_N12/shear_56.dat                                                              # mu(gamma), shear_curve.png
```

`RUN_N12/dehydrate.dat` lists the water number, m_w/m_c, dry bulk density and Lz of every dehydration step; `RUN_N12/shear_<σ>.dump` holds the clay atoms (unwrapped, every 0.1 ns = strain 0.01) for the pore-size and stacking analysis; it is the only trajectory
(150 k – 1.5 M atoms: a full frame is 7–70 MB). Structures to look at: `dehyd_snap_<k>.data` (every 5 dehydration cycles), the `.data` file at the end of each
stage, and `<stage>.latest.data` at each job stop.

## Notes and limits

- σn in the input is in MPa (1 MPa = 9.869 atm); the box in x and y is fixed during stages 03-05, z is the free direction (NPzT, as in the paper).
- 6 nm particles are smaller than natural grains and carry many edges per area: read the result as the limit in which edges matter most.
- A real gouge compacted at 5.6 or 56 MPa is far from the dry bulk density reached by a short dehydration; `final_mwmc` (default 0.2) and the water removal speed decide what is sheared, so report them with μ.
- Not done: edge protonation states other than the one of the library, particles of other minerals, Si-O bond breaking.
