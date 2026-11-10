#!/bin/bash
# One command: build the system and submit the whole PBS chain.
#   ./go.sh [N particles = 12] ["normal stresses in MPa" = "5.6 56"] [stages = all]
#   atoms ≈ 22,300 × N  (N=12: 0.27 M, N=67: 1.49 M).  NODES=<n> ./go.sh 67  uses n nodes.
# Before the first run, put the supercomputer's settings in config.env (copy them from /home/app/lammps/run.sh).
set -eu
cd "$(dirname "$0")"
N=${1:-12}; SIG=${2:-"5.6 56"}; ST=${3:-all}; RUN=RUN_N$N
PY=${PYTHON:-python3}
# random configurations to try (the most random one is kept): 300 for small systems; 30 for N > 20 (keeps the login node busy for minutes, not hours)
NCONF=${NCONF:-$([ "$N" -gt 20 ] && echo 30 || echo 300)}
$PY -c "import numpy, scipy" 2>/dev/null || { echo "installing numpy and scipy for $PY (pip --user)"; $PY -m pip install --user -q numpy scipy; }
[ -f $RUN/gouge.data ] || $PY build_gouge.py -n $N --mwmc 3.0 --nconf $NCONF --seed 1 -o $RUN/gouge
if [ -f submit_chain.sh ]; then S=./submit_chain.sh; else S=pbs/submit_chain.sh; fi
$S $RUN "$SIG" "$ST"
echo "submitted. watch with: qstat -u \$USER   results in $RUN/ (shear_<sigma>.dat; python3 shear_curve.py $RUN/shear_5.6.dat)"
