#!/bin/sh
# Na-montmorillonite 2W (../initial.xyz) with the plain ClayFF + flexible SPC force field built by soapff-cui (bundled with
# soapff-gui) from the bundled library: minimisation, 10 ps relaxation NVT, 10 ps NVT, 100 ps NPT along z (P = 0), 300 K, 0.5 fs.
# (start/ was made with the earlier 5-structure library; the types and charges of all 2,502 atoms are the same with the
#  13-structure library bundled now.)   usage: sh make_start.sh   [PY=python LMP=lmp MPIRUN=mpirun NP=3]
PY=${PY:-python3}; LMP=${LMP:-lmp}; MPIRUN=${MPIRUN:-mpirun}
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE/../../.."             # repository root
$PY -W ignore -m soapff_cui make examples/montmorillonite_2W/initial.xyz -o examples/montmorillonite_2W/diffusion/start \
    --library src/soapff_gui/data/library_clayff_bulk.npz --relax-ps 10 --nvt-ps 10 --npt-ps 100 --baro z --dump-ps 10
cd examples/montmorillonite_2W/diffusion/start && OMP_NUM_THREADS=1 $MPIRUN -np ${NP:-3} $LMP -in in.lammps -log log.lammps -screen none
