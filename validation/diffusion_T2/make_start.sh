#!/bin/sh
# Na-montmorillonite 2W (JAEA initial structure) with the plain ClayFF + flexible SPC force field built by soapff-cui from
# library_clayff_bulk.npz: minimisation, 10 ps relaxation NVT, 10 ps NVT, 100 ps NPT along z (P = 0), 300 K, 0.5 fs.
PY=${PY:-/Users/harutokono/ScoreMD/.pixi/envs/default/bin/python}
BIN=${BIN:-/Users/harutokono/ScoreMD/.pixi/envs/default/bin}
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE/../.."
$PY -W ignore -m soapff_cui make examples/jaea/initial.xyz -o examples/diffusion_T2/start --library examples/minerals/library_clayff_bulk.npz \
    --relax-ps 10 --nvt-ps 10 --npt-ps 100 --baro z --dump-ps 10
cd examples/diffusion_T2/start && OMP_NUM_THREADS=1 FI_PROVIDER=tcp $BIN/mpirun -np ${NP:-3} $BIN/lmp -in in.lammps -log log.lammps -screen none
