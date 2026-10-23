#!/bin/sh
# Na-montmorillonite 2W (../initial.xyz: 3x2x1 supercell of the 32-water state of examples/montmorillonite_swelling, 1,554 atoms)
# with the plain ClayFF + flexible SPC force field built by soapff-cui (bundled with soapff-gui) from the bundled library:
# minimisation, 10 ps relaxation NVT, 10 ps NVT, 100 ps NPT along z (P = 0), 300 K, 0.5 fs.
# usage: sh make_start.sh   [PY=python LMP=lmp MPIRUN=mpirun NP=3 MPIOPTS=...]   (RUN_ONLY=1: skip the force-field generation)
PY=${PY:-python3}; LMP=${LMP:-lmp}; MPIRUN=${MPIRUN:-mpirun}
HERE=$(cd "$(dirname "$0")" && pwd)
if [ -z "$RUN_ONLY" ]; then
  cd "$HERE/../../.."             # repository root
  $PY -W ignore -m soapff_cui make examples/montmorillonite_2W/initial.xyz -o examples/montmorillonite_2W/diffusion/start \
      --library src/soapff_gui/data/library_clayff_bulk.npz --relax-ps 10 --nvt-ps 10 --npt-ps 100 --baro z --dump-ps 10
fi
cd "$HERE/start" && OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl $MPIRUN $MPIOPTS -np ${NP:-3} $LMP -in in.lammps -log log.lammps -screen none
