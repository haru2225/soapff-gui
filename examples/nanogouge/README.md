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
| 02 | `lammps/in.02_dehydrate` | every `step_ps` delete 2 % of the initial water (1 % below m_w/m_c = 1; each molecule with that probability, any LAMMPS version) and run NPzT: compaction, particles stack (paper: 1 % every 1 ns, 134 ns); stops at `final_mwmc` |
| 03 | `lammps/in.03_consolidate` | hold at the normal stress σn (z only; x, y fixed) |
| 04 | `lammps/in.04_shear` | shear by the cell deformation (box tilt xz, rate 1e-7 /fs = 1e8 /s; Lees-Edwards boundaries, no wall) at constant σn (Nose-Hoover on z, x and y fixed), Langevin thermostat on the velocity relative to the linear profile (Lin's thesis, Osaka 2024, sect. 2.2; `FLOW=1`, default); then 4 velocity steps up to 1e10 /s; `shear_<σ>.dat` = strain, shear stress, normal stress, μ, Lz, ρb, T, Lx, Ly, rate; figures in `figs/` after each job; particle and velocity-profile files for the slip analysis; restart files every 10 ps and at each job stop |
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
Python (numpy and scipy for `build_gouge.py`, matplotlib for the figures) is `./sifpy`: python 3.11 in a Singularity / Apptainer image that `make_sif.sh` makes once (go.sh calls it when it is missing: pulls python:3.11-slim, no root; the image's pip installs numpy, scipy, matplotlib into `sifpy_lib/`; a few minutes; no pip on the centre's python). Without internet on the login node: run make_sif.sh on a computer with internet and Singularity, copy `sifpy.sif` and `sifpy_lib/` here, run it again. `PYTHON=/path/python3 ./go.sh` uses another python. `./go.sh 12 "5.6 56" "01 02"` = N particles, normal stresses in MPa, stages (here only the compaction).
`./go.sh 67 "5.6 56" all 4 RUN_N67_b` (5th argument = run folder) runs the same calculation again in another folder: the initial structure, `config.env` and `lmp.auto`
of `RUN_N67` are copied. `./go.sh 67 "5.6 56" all 4 RUN_N67_s2 2` (6th argument = seed, default 1) is an independent run: another random initial structure,
and LAMMPS uses the seed too (velocities, water removal, thermostat; `RUN_N67_s2/lammps_seed`). Each folder has its own chain and can run at the same time as the others.
Wall time of the jobs: `env WALLTIME=3 ./go.sh 67 "5.6 56" all 4` (hours, or `3:00:00`) replaces every `walltime` of the qsub options for this folder and keeps it in
`RUN_N67/walltime` (later jobs and `./go.sh` runs use it). The waiting jobs of this folder (old wall time) are deleted and new ones queued after the job working now;
the running job goes on. The queue must allow that wall time: `qstat -Qf <queue> | grep walltime`.
`./pbs/jobs.sh` lists your PBS jobs with the run folder each one works in; `./pbs/jobs.sh keep RUN_N67` keeps only the jobs of RUN_N67 with its wall time
(`RUN_N67/walltime`) and deletes all your other jobs, also RUN_N67 jobs with an old wall time, running ones too (it asks first).
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
output files `ng_*.o*` appear only when a job has ended). The first job finds a LAMMPS that starts on the node, preferably one with `fix press/langevin` (needed only by `FLOW=0`; the centre's old
`lmp_icex` 3Mar20 has none: the notes about it in the log are normal; the default shear runs on any version) and checks that `mpirun` spreads the MPI ranks over all nodes (HPE MPT: else
`mpiexec_mpt`); both are kept in `RUN_N12/lmp.auto` for the next jobs (delete it to search again). Each stop at the wall time also writes the current
structure to `<stage>.latest.data` (e.g. `01_equil.latest.data`; opens in OVITO).
Trajectories (OVITO; clay atoms and Na, unwrapped, `id mol type xu yu zu`; each job appends to the file of its stage): `01_equil.dump` every 5 ps, `02_dehydrate.dump` every 20 ps,
`03_consolidate_<σ>.dump` every 10 ps, `shear_<σ>.dump` every 0.01 of strain (`dump_dg`; 100 ps at 1e8 /s), `heal_<σ>.dump` every 25 ps (`dump_ps`; the steps follow the time step).
With the water (all atoms, ~10× larger files): `env DUMP=all ./go.sh ...`. Time step 0.5 fs instead of 1: `env DT=0.5 ./go.sh ...` (all times stay in ps, so there are twice as many steps).
Both are kept in `RUN_N<N>/config.env` (`EXTRA=... # go.sh: ...` lines) for every job of the folder; other LAMMPS variables the same way, e.g. `EXTRA="${EXTRA:-} -var dump_ps 0"` (no trajectories).
`env MWMC=0.5 ./go.sh ...` stops the dehydration (02) at m_w/m_c = 0.5 instead of 0.2 (the water content of the sheared gouge; one folder per water content: 5th argument).
Example, about 0.5 M atoms (N = 22: 490,599 atoms), 0.5 fs, water in the trajectories, 30 min jobs: `env WALLTIME=0:30:00 DT=0.5 DUMP=all ./go.sh 22 "5.6 56" all 4`. New `in.*` files reach a run folder with `./go.sh` (or `submit_chain.sh`);
they replace the old files, so a LAMMPS running now is not disturbed and the next job uses them.
Stage 01 is 100 ps NVT + 1000 ps NPzT (1.1 M steps); `EXTRA="-var t_npt_ps 200"` in `RUN_N67/config.env` shortens it (a run already past the new end finishes the stage at once).
Each job prints the step it reached and the speed in `jobs.log` (`grep -E "== stage|   step" RUN_N67/jobs.log | tail`); stage 01 is done at step `(t_nvt_ps + t_npt_ps) × 1000`, then the same job goes on with 02.
If the centre does not allow `depend`, submit the jobs by hand in the order of the table: `cd RUN_N12; qsub -v STAGE=01 job.sh`, etc. (`job.sh` reads `STAGE`, `SIGMA`).

Size (N = 12, m_w/m_c = 3): about 12 × 2,292 + 3 × 80,000 ≈ 268,000 atoms at the start (the builder takes seconds; S ≈ 0.13 of 300 configurations) and ≈ 40,000 after the dehydration to m_w/m_c = 0.2.
Defaults are small enough for a laboratory-scale first try; for a larger system use `-n 30` and the 24 × 24 × 40 nm box of the paper is approached.
Shear to γ = 3 at 1e8 /s is 30 ns (30 million steps of 1 fs); for a quick look use `-var rate 5.0e-7 -var gmax 1.0` through `EXTRA` (`qsub -v STAGE=04,SIGMA=5.6,EXTRA="-var rate 5.0e-7 -var gmax 1.0" job.sh`).

## Slip between the clay particles: "hot spots" (Zhou et al., Am. Mineral. 108, 87, 2023)

Zhou et al. pull one clay layer over another with a constant force and find that the water between the layers shears like a viscous,
shear-thinning fluid (ln τ = k ln v + ln τi, velocity strengthening): the interlayers would be the hot spots where the plastic deformation of a gouge
concentrates. Here the same is looked at in the whole assemblage: the gouge is sheared by the cell deformation only (no wall, the slip can concentrate
anywhere), and every particle is followed in time. `in.04_shear` writes, besides `shear_<σ>.dat`:
`particles_<σ>.<step>.dat` (every `part_dg` = 0.002 of strain: centre of mass, velocity, gyration tensor = normal and size, internal temperature of each
particle) and `prof_<σ>.<step>.dat` (averaged over `prof_dg` = 0.05 of strain, 40 layers across the gouge: vx of the water and of the clay, their numbers, T).

| shear (`env FLOW=... ./go.sh`) | drive | thermostat | normal stress |
|---|---|---|---|
| 1 (default) | cell deformation, Lees-Edwards (velocity ramp at the start) | Langevin on the deviation from the linear profile (`temp/deform`): Lin's thesis | Nose-Hoover on z (`fix nph z`), x, y fixed |
| 2 | same | Langevin on vy only (does not pull a localized profile back to the linear one) | same |
| 3 | same | none (the frictional heat stays) | same |
| 0 | SLLOD (imposed homogeneous flow) | Nose-Hoover on all atoms | `fix press/langevin` (LAMMPS 21Nov2023 or newer) |

**One chain gives all of it, with the usual command** (e.g. `env WALLTIME=0:30:00 ./go.sh 67 "5.6 56" all 4`): after the shear at 1e8 /s to γ = 3
(phase 0; its end is the start of the healing 05), `in.04_shear` steps the rate up 4 times by ×3.16 (3.2e8, 1e9, 3.2e9, 1e10 /s; 0.5 of strain each,
~8 % more time): the shear stress against the slip velocity (their Fig. 8, 10) from one run, like a velocity-step test of a friction experiment.
At the start and at the end of each phase it writes the whole system with the water (`shear_<σ>_all_<step>.dump`, for their Fig. 9; the trajectory
`shear_<σ>.dump` stays clay only). After every job of the shear, the job draws the figures into `RUN_N67/figs/` (`hot_spots.py`, copied into the run folder;
a trajectory frame is analysed once and kept in `figs/.cache`, so the jobs stay quick). They use `./sifpy` (see above), which every `RUN_*/config.env` gets as `PYTHON=...` (`./make_sif.sh` run later: a running chain draws the figures from its next job on). A shear begun
with an older `in.04_shear` (no `shear_<σ>.v2`) is done again: its files, and those of the healing after it, are moved to `old_04_<σ>/`.
Options: `env NVS=0` (no velocity steps), `NVS=2`, `VFAC=10` (×10 each), `RATE=1e9` (1/s), `GMAX=2`.

Everything in one input file, production lengths: `lammps/in.all` (dt 1 fs) runs 01 (minimisation, NVT 0.1 ns, NPzT 0.1 ns), 02 (2 ns: removals of
5 % of the initial water down to m_w/m_c = 0.2, twice the time per removal below m_w/m_c = 1, the form of Underwood & Bourg 2020, who took 134 ns;
none with `DRYNS=0`: the shear at the water content of the structure, `MWMC0`), then at each of the two normal stresses 03 (NPzT in pieces of 50 ps until Lz changes by less
than 0.1 % over a piece, 0.1-1 ns; `03_consolidate_<s>.pieces`), 04 (1e9 /s to the strain 3, then 3.2e9 and 1e10 /s, 0.5 of strain each) and 05
(hold 50 ps): about 11 ns; after 02 the system is small (N=22, m_w/m_c 0.2: about 80,000 atoms). With PBS:
`env STACK=3 MWMC0=1.5 WALLTIME=12 ./go.sh 22 "5.6 56" one 4 RUN_N22_fol` (a gouge of stacks built at m_w/m_c 1.5, compacted to 0.2 by 02, sheared;
each run ends with `all.restart` and `all.state`, so the next job goes on). Lengths: `NEQ`, `DRYNS` (02 in ns), `DFRAC`, `NCONS`, `CONSTOL`,
`NCONSMAX`, `RATE`, `GMAX`, `NVS`, `GSTEP`, `HOLD`; `MWMC` = m_w/m_c at the end of 02, `MWMC0` = of the new structure; `NFAC=2` doubles 02-05;
`NRUN=2000`: a test, every run 2000 steps. The speed: `grep Performance RUN_N22_fol/log.one | tail -1` (ns/day).

Dense gouges at a set water content, two fabrics (as Sakuma et al. 2023: nominally randomly oriented NRO and preferentially oriented PO):
`MWMC0` = 0, 0.1, 0.2 (dry, 10 %, 20 %) builds a loose structure (the box of m_w/m_c 1.5 for po, 2 for nro) with the water on the clay surfaces and in the
interlayers first; LAMMPS compacts it in 01 at the first normal stress (`ORIENT=nro`: particles at random in a cube, NPT iso; `ORIENT=po`:
the layers within 10 deg of the shear plane, NPzT), with no dehydration. N=54 gives a compacted gouge of about two particle diameters each
way: `env MWMC0=0 ORIENT=nro WALLTIME=12 ./go.sh 54 "5.6 56" one 4 RUN_N54_w0_nro` (dry 0.12 M atoms, 20 %: 0.20 M). `NOSUB=1` builds
only: `RUN_N54_w0_nro/gouge.xyz` for OVITO (`analysis/data2xyz.py` converts any data file).

Fabrics made by a dehydration, as Underwood & Bourg 2020, and the shear of snapshots taken on the way (two steps). 1) `FABRIC=dis`: the
particles at random in a cube (m_w/m_c 2), dehydrated the same in x, y, z (NPT iso): removals of 1 % of the initial water (`DFRAC`), twice
the time per removal below m_w/m_c 1, 100 ns in all (`DRYNS`; about their rate, 1 % per ns; they took 134 ns). The first state at each m_w/m_c of
`WSNAP="0.5 0.3 0.2 0.1 0"` is written, `02_w<value>.data` (listed in `snapshots.dat`), and the calculation ends at the last one; no shear:
`env FABRIC=dis WALLTIME=24 ./go.sh 54 "5.6" one 4 RUN_N54_fab_dis` (0.84 M atoms at the start). `FABRIC=ori`: the particles placed as
sedimented (layers within 30 deg of the horizontal, `TILT`) in a tall box (m_w/m_c 1.5), dehydrated along z only (NPzT with x, y fixed: a
uniaxial dehydration). `analysis/stacks.py` follows the stacking and the order along the dehydration (02_dehydrate.dump).
2) `FROMFAB=RUN_N54_fab_dis MWMC0=0.2`: the shear of the snapshot `02_w0.2.data` with its own water (it can start while the dehydration goes
on): minimisation, NVT 0.1 ns, the slow compaction at the first normal stress (dis: x, y, z; ori: along z), 0.3 ns NPT, 03, 04, 05 as for
the dense gouges: `env FROMFAB=RUN_N54_fab_dis MWMC0=0.2 WALLTIME=24 ./go.sh 54 "5.6 56" one 4 RUN_N54_dis_w0.2`. A water content without
a snapshot (or `REWATER=packmol`), once the dehydration is done: `rewater.py` takes all the water out of `02_dehydrate.data` and puts MWMC0
in with packmol (`./get_packmol.sh` makes `bin/packmol`; without it a grid, the sites nearest to the clay first).

GPU: `env GPU=1 ...  ./go.sh ...` runs the jobs of the folder on a GPU node (GPUs per node) with a LAMMPS of the GPU or KOKKOS package
(`pbs/job.sh` looks for one in /home/app/lammps*, or `GPU_LMP=/path/lmp`). The qsub options: `GPU_QSUB="-q <queue> -l
select=1:ncpus=..:mpiprocs=..:ngpus=1"`, else the line `bash pbs/gpu_check.sh` suggests (it lists the queues, the GPU nodes, the GPU builds of
LAMMPS and the CUDA modules); `GPU_MODULES="module load cuda"`, `GPU_PKG=gpu|kokkos`, `GPU_RPG=8` (MPI ranks per GPU, GPU package). The first
job checks with a tiny run that LAMMPS works on the GPU and whether its pppm takes a triclinic box (LAMMPS 3Mar2020: no; then the shear and
the hold run pppm on the CPU, the rest on the GPU); kept in `lmp_gpu.auto`. `GPU=0`: the CPU again.

The jobs of a folder form a chain: go.sh submits 10 jobs of WALLTIME, each starting when the one before has ended and continuing from its
restart (all.restart, all.state), and the last one submits more while work is left, until the calculation is done. `MAXJOBS=5` limits
it: from this submission on at most 5 jobs (a job running now not counted; its waiting jobs are replaced), then it stops; the same
command again gives 5 more (`MAXJOBS=0`: no limit). `MAXHOURS=72` does the same in hours.

MPI x OpenMP: `env OMP=2 ... ./go.sh ...` (stages "one") runs 2 OpenMP threads per MPI rank: the same cores, half the ranks (40 cores per
node: `select=..:ncpus=40:mpiprocs=20:ompthreads=2`), LAMMPS `-sf omp -pk omp 2` (its OPENMP package; go.sh checks the LAMMPS of the folder),
HPE MPT `omplace -nt 2` places the threads on cores of their own; the shear and the hold run plain pppm (pppm/omp may not take a triclinic
box). Fewer ranks: less communication (PPPM); faster only when that takes much of the time (`grep -A12 "MPI task timing" log.one`).
`OMP=1`: MPI only again. A change of the qsub options (nodes, OMP) replaces the waiting jobs of the folder.

Stacks: the particles start apart at random (U&B 2020) and come together face to face during the slow dehydration (U&B: S from 0.1 to 0.7-0.8 over
134 ns; the particles rotate in tens of ns). `env STACK=3 ./go.sh ...` builds them in stacks
(tactoids) of 3 at 15.5 A (2 water layers; `D001=12.5`: 1 layer), the stacks at random (structure in `RUN_N<N>_stack3`). Whether they are stacked:
`./figs.sh RUN_N22_prod` (also before the shear) runs `analysis/stacks.py`: stacked pairs, particles in stacks, film and order S at the end of 01, along 02,
at the end of 03 and along 04 (`figs/stacks.png`, `stacks.csv`), and the foliation forming in the shear against the strain: S, tilt of the
particles from the shear plane, share of particles within 20 deg of it, particles in stacks (`figs/foliation.png`).

### Shear at several water contents (porosities): `go_water.sh`, Lin & Hatano (2024)-type figures

`env WALLTIME=0:30:00 ./go_water.sh 22 "2 1 0.5 0.2" "5.6 56" 1` (N, m_w/m_c from wet to dry, normal stresses, nodes; folder `RUN_N22_water`).
`lammps/in.water` dehydrates once from the wet start (m_w/m_c = 3) and stops at each listed water content: the gouge of that moment goes to
`w<m_w/m_c>/` (`02_dehydrate.data`), is consolidated (03), sheared (04: 1e9 /s to the strain `GMAX` = 1, then `NVS` = 2 velocity steps ×3.16 of
`GSTEP` = 0.5 strain each, to 1e10 /s) and held at σn (05, `HOLD` = 20000 steps), at each normal stress; then the dehydration goes on.
`water_contents.dat` lists the reached m_w/m_c, water number, dry bulk density and Lz. Porosity of water-filled pores φ = (m_w/m_c) / (m_w/m_c + ρw/ρs):
m_w/m_c 2, 1, 0.5, 0.2 → φ 0.84, 0.73, 0.57, 0.35. `env BY=phi ./go_water.sh 22 "0.85 0.7 0.55 0.4"` takes porosities. A test of everything with every
run 1000 steps: `env NRUN=1000 WALLTIME=0:30:00 ./go_water.sh 22 "2 0.5" "5.6 56" 1` (folder `RUN_N22_water_test1000`).
The jobs draw (`analysis/lin_figs.py`; by hand `./sifpy analysis/lin_figs.py RUN_N22_water`) in `RUN_N22_water/figs`:
`lin_fig5_rate.png` (μ and φ against the strain rate per water content), `lin_fig6_profiles_<σ>.png` (velocity profile per rate and particle-normal
angle to z across the gouge), `lin_fig7_stacking.png` (parallel radial distribution function g∥(r) at the start / end of the shear, stacked neighbours
per particle against the strain), `lin_fig8_healing.png` (stacked neighbours against the rest time in the hold, a + b log10 t; b against φ),
`lin_fig9_a_b.png` ((a − b) = Δμss / ln(V2/V1) between the velocity steps against φ), `lin_fig10_band.png` (R² of the velocity profile against the
strain; share of the height that shears against φ), `lin_porosity.png` (μ, (a − b), R², stacking against φ), numbers in `lin_points.csv`, `lin_heal.csv`;
and the hot-spot figures of each water content in `w<m_w/m_c>/figs`.

A quick test of the whole chain (all stages, all outputs, the figures): `env NRUN=1000 WALLTIME=0:30:00 ./go.sh 22 "5.6 56" all 4` runs every
run of every stage for 1000 steps (02: one water removal; 04: each of its 5 phases; outputs every 100 steps) in its own folder `RUN_N22_test1000`
(the structure of `RUN_N22` is copied when it is there). The numbers mean nothing; it shows that everything runs and is written.

To draw them by hand from the results so far: `./figs.sh RUN_N22_one RUN_N22_water` (no argument: every `RUN_*` with a shear; `env JOB=1 ./figs.sh ...` as a short PBS job).
It finds or makes the python (`./sifpy` of `make_sif.sh`, else a python of the machine with numpy and matplotlib: PATH, modules, `/home/app`, other `.sif`
images, else the node's python3 + matplotlib from pip in `pyfigs_lib/`; what was tried: `figs_python.log`), writes it into the run folders' `config.env`
(their next jobs draw the figures themselves) and draws hot_spots.py (and lin_figs.py for a water series). Also several folders together:
`./sifpy analysis/hot_spots.py RUN_N67 --outdir figs` (any other computer: `python3 ...`).
Several folders at other rates from the same compacted gouge are possible too (each its own chain):
`env FROM=RUN_N22 RATE=1e9 GMAX=2 ./go.sh 22 "5.6" all 4 RUN_N22_r1e9` (after `RUN_N22/03_consolidate_5.6.done`; `FROM` copies its stages 01-03 and settings),
then `./sifpy analysis/hot_spots.py RUN_N22 RUN_N22_r1e9 --sigma 5.6`.

The slip velocity across the box is V = γ̇ Lz (1e8 /s with Lz = 10 nm: 1 m/s). `hot_spots.py` writes
`fig6_interface_slip.png` (relative in-plane displacement of stacked particle pairs in time over all phases, the fastest, median and slowest interface at
the first rate, and how much faster than the gouge the water films shear: (v/h)/γ̇), `fig7_velocity_profile.png` (vx / (γ̇ Lz) of water and clay against
z / Lz in the steady state of each rate: the gray diagonal is uniform shear, a step is a shear band; where the water is, T(z)), `fig8_tau_v.png` (τ against
the slip velocity of the gouge and of the fastest 10 % of the interfaces, one point per rate, log-log fit ln τ = k ln v + ln τi), `fig9_water_structure.png`
(hydrogen bonds per water molecule with their criteria, and g_OO(r) of the water between two particles, at one surface and in the pores, at each
snapshot) and `fig10_viscosity.png` (τ / shear rate against the shear rate), with the main numbers in `hot_spots_summary.txt` and the numbers behind each
figure in `fig6_interfaces.csv`, `fig7_profiles.csv`, `fig8_fig10_points.csv` and `fig9_water.csv` (to plot them again elsewhere). The steady state of
each rate is the second half of its strain (`--gmin` for phase 0).

## Analysis (on the supercomputer `./sifpy` in place of `python3`; any other computer: numpy and matplotlib; `analysis/` in the repository, next to `go.sh` in the pack)

```
python3 order_parameter.py RUN_N12/01_equil.data RUN_N12/02_dehydrate.data RUN_N12/03_consolidate_5.6.data RUN_N12/gouge.json   # S, m_w/m_c, dry bulk density
python3 shear_curve.py RUN_N12/shear_5.6.dat RUN_N12/shear_56.dat                                                              # mu(gamma), shear_curve.png
python3 hot_spots.py RUN_N12                                                                                                   # slip between the particles (see above)
```

`RUN_N12/dehydrate.dat` lists the water number, m_w/m_c, dry bulk density and Lz of every dehydration step; `RUN_N12/shear_<σ>.dump` holds the clay atoms (unwrapped, every 0.01 of strain, `dump_dg`) for the pore-size and stacking analysis; the other stages have their own `.dump` (see above)
(150 k – 1.5 M atoms: a frame of all atoms is 7–70 MB; of the clay atoms only, ~10× less). Structures to look at: `dehyd_snap_<k>.data` (every 5 dehydration cycles), the `.data` file at the end of each
stage, and `<stage>.latest.data` at each job stop.

## Notes and limits

- σn in the input is in MPa (1 MPa = 9.869 atm); the box in x and y is fixed during stages 03-05, z is the free direction (NPzT, as in the paper).
- 6 nm particles are smaller than natural grains and carry many edges per area: read the result as the limit in which edges matter most.
- A real gouge compacted at 5.6 or 56 MPa is far from the dry bulk density reached by a short dehydration; `final_mwmc` (default 0.2) and the water removal speed decide what is sheared, so report them with μ.
- Not done: edge protonation states other than the one of the library, particles of other minerals, Si-O bond breaking.
