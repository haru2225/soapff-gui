#!/bin/sh
# Diffusion at 300, 330 and 366 K (Greathouse et al. 2016 conditions) with the force field built by soapff-cui (make_start.sh),
# continued from start/final.data (same types, charges and bonds), new velocities per temperature:
#   NVT 20 ps + NPT along z (P = 0) 100 ps (equilibration, 0.5 fs), then 4 segments of NPT (z, P = 0) 200 ps each; every segment
#   defines a new compute msd (time origin at the start of the segment); D_xy of each segment from 50-200 ps.
# usage: make_T.sh T      (runs in T$T/; LMP=lmp MPIRUN=mpirun NP=4)
LMP=${LMP:-lmp}; MPIRUN=${MPIRUN:-mpirun}
T=$1
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
mkdir -p T$T/equil
cp start/final.data T$T/equil/start.data
sed 's/^read_data .*/read_data start.data/' start/system.settings > T$T/equil/system.settings
cp start/charges.csv start/report.json T$T/equil/
MSD=$(grep -E "^(group a_|compute msd_|fix msd_)" start/in.lammps)
DUMP=$(grep -E "^(dump d |dump_modify d )" start/in.lammps | sed 's/^dump d all custom [0-9]* /dump d all custom 20000 /')
cat > T$T/equil/in.lammps <<EOT
include system.settings
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
thermo_style custom step temp press pe ke etotal vol density lx ly lz
thermo 2000
reset_timestep 0
velocity all create $T.0 4928459 mom yes rot no dist gaussian
timestep 0.5
fix md all nvt temp $T.0 $T.0 100.0
run 40000
unfix md
fix md all npt temp $T.0 $T.0 100.0 z 0.0 0.0 1000.0
run 200000
write_data final.data
EOT
(cd T$T/equil && OMP_NUM_THREADS=1 $MPIRUN -np ${NP:-4} $LMP -in in.lammps -log log.lammps -screen none)
prev=T$T/equil/final.data
for s in 1 2 3 4; do
  d=T$T/seg$s; mkdir -p $d
  cp $prev $d/start.data
  sed 's/^read_data .*/read_data start.data/' start/system.settings > $d/system.settings
  cp start/charges.csv start/report.json $d/
  cat > $d/in.lammps <<EOT
include system.settings
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
thermo_style custom step temp press pe ke etotal vol density lx ly lz
thermo 2000
reset_timestep 0
timestep 0.5
fix md all npt temp $T.0 $T.0 100.0 z 0.0 0.0 1000.0
$DUMP
$MSD
run 400000
write_data final.data
EOT
  (cd $d && OMP_NUM_THREADS=1 $MPIRUN -np ${NP:-4} $LMP -in in.lammps -log log.lammps -screen none)
  prev=$d/final.data
done
