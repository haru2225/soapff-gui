#!/bin/bash
# One command: build the system and submit the whole PBS chain.
#   ./go.sh [N particles = 12] ["normal stresses in MPa" = "5.6 56"] [stages = all]
#   atoms ≈ 22,300 × N  (N=12: 0.27 M, N=22: 0.49 M, N=67: 1.49 M).  ./go.sh 67 "5.6 56" all 4   uses 4 nodes (4th argument).
#   ./go.sh 67 "5.6 56" all 4 RUN_N67_b      the same run in another folder (5th argument; the initial structure of RUN_N67 is copied)
#   ./go.sh 67 "5.6 56" all 4 RUN_N67_s2 2   an independent run: another random initial structure and LAMMPS seed (6th argument, default 1)
#   env DT=0.5 DUMP=all ./go.sh 22 "5.6 56" all 4   time step 0.5 fs (default 1) and the water in the trajectories (default: clay and Na only);
#   env MWMC=0.5 ./go.sh ...                 the dehydration (02) stops at water/clay = 0.5 (default 0.2): a wetter gouge
#                                            both are kept in RUN_N22/config.env for every job of the folder (WALLTIME=... too, in RUN_N22/walltime)
# Before the first run, put the supercomputer's settings in config.env (copy them from /home/app/lammps/run.sh).
set -eu
cd "$(dirname "$0")"
N=${1:-12}; SIG=${2:-"5.6 56"}; ST=${3:-all}; RUN=${5:-RUN_N$N}; SEED=${6:-1}; BASE=RUN_N$N
[ -n "${4:-}" ] && export NODES=$4      # 4th argument = number of nodes (works in any shell, no "NODES=4 ./go.sh" needed)
PY=${PYTHON:-python3}
# random configurations to try (the most random one is kept): 300 for small systems; 30 for N > 20 (keeps the login node busy for minutes, not hours)
NCONF=${NCONF:-$([ "$N" -gt 20 ] && echo 30 || echo 300)}
$PY -c "import numpy, scipy" 2>/dev/null || { echo "installing numpy and scipy for $PY (pip --user)"; $PY -m pip install --user -q numpy scipy; }
if [ ! -f "$RUN/gouge.data" ]; then
  mkdir -p "$RUN"
  if [ "$SEED" = 1 ] && [ "$RUN" != "$BASE" ] && [ -f "$BASE/gouge.data" ]; then
    # seed 1 builds the structure of RUN_N<N> again: copy it (minutes saved)
    cp "$BASE/gouge.data" "$BASE/gouge.settings" "$RUN/"; [ -f "$BASE/gouge.json" ] && cp "$BASE/gouge.json" "$RUN/"
    echo "$RUN: the same initial structure as $BASE (copied)"
  else
    $PY build_gouge.py -n $N --mwmc 3.0 --nconf $NCONF --seed $SEED -o "$RUN/gouge"
  fi
  # another seed: LAMMPS uses it too (velocities, water removal, thermostat); job.sh reads lammps_seed
  [ "$SEED" != 1 ] && echo "$SEED" > "$RUN/lammps_seed"
fi
# a new folder takes the settings the jobs of RUN_N<N> found and checked (your config.env, the LAMMPS build and MPI launcher in lmp.auto)
if [ "$RUN" != "$BASE" ]; then for f in config.env lmp.auto; do [ -f "$BASE/$f" ] && [ ! -f "$RUN/$f" ] && cp "$BASE/$f" "$RUN/"; done; fi
if [ -f submit_chain.sh ]; then S=./submit_chain.sh; CE=./config.env; else S=pbs/submit_chain.sh; CE=pbs/config.env; fi
# DT / DUMP: one line each in $RUN/config.env (job.sh passes EXTRA to LAMMPS); given again, the line is replaced
set_opt() {
  local f=$RUN/config.env; [ -f "$f" ] || cp "$CE" "$f"
  { grep -v "# go.sh: $1\$" "$f" || true; echo "EXTRA=\"\${EXTRA:-} -var $1 $2\"   # go.sh: $1"; } > "$f.tmp" && mv -f "$f.tmp" "$f"
  echo "$RUN: LAMMPS -var $1 $2 for every job of this folder ($f)"
}
if [ -n "${DT:-}" ]; then set_opt dt "$DT"; fi
if [ -n "${DUMP:-}" ]; then set_opt dump_group "$DUMP"; fi
if [ -n "${MWMC:-}" ]; then set_opt final_mwmc "$MWMC"; fi   # water / clay mass ratio at which the dehydration (02) stops (default 0.2)
$S $RUN "$SIG" "$ST"
echo "submitted. watch with: qstat -u \$USER;  bash $RUN/status.sh $RUN   results in $RUN/ (shear_<sigma>.dat; python3 shear_curve.py $RUN/shear_5.6.dat)"
