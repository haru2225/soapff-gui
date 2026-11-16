#!/bin/bash
# One command: build the system and submit the whole PBS chain.
#   ./go.sh [N particles = 12] ["normal stresses in MPa" = "5.6 56"] [stages = all]
#   atoms ≈ 22,300 × N  (N=12: 0.27 M, N=22: 0.49 M, N=67: 1.49 M).  ./go.sh 67 "5.6 56" all 4   uses 4 nodes (4th argument).
#   ./go.sh 67 "5.6 56" all 4 RUN_N67_b      the same run in another folder (5th argument; the initial structure of RUN_N67 is copied)
#   ./go.sh 67 "5.6 56" all 4 RUN_N67_s2 2   an independent run: another random initial structure and LAMMPS seed (6th argument, default 1)
# Options (env VAR=... ./go.sh ...; kept in RUN_<..>/config.env for every job of the folder, WALLTIME in RUN_<..>/walltime):
#   DT=0.5        time step in fs (default 1)
#   DUMP=all      the water in the trajectories too (default: clay and Na only)
#   MWMC=0.5      the dehydration (02) stops at water/clay = 0.5 (default 0.2)
#   FLOW=1        how the shear (04) is driven: 1 (default) Lees-Edwards + Langevin on the deviation from the linear profile + barostat on z
#                 (Lin's thesis); 2 the same with the Langevin thermostat on vy only; 3 no thermostat; 0 SLLOD
#   RATE=1e9      shear rate of 04 in 1/s (default 1e8)
#   GMAX=2        final shear strain of 04 (default 3)
#   NVS=4         velocity steps after GMAX (default 4: the rate x3.16 each, 0.5 of strain each; 1e8 -> 1e10 /s; 0 = none), VFAC=3.16227766 their factor
#   The jobs draw the figures (analysis/hot_spots.py) in RUN_<..>/figs after each job of the shear; a shear begun with an older input is done again.
#   FROM=RUN_N22  start from the compacted gouge of another folder (its stages 01-03 are copied; only 03 of a missing sigma, 04, 05 run):
#                 a series of shear rates on the same gouge, e.g.
#                 env FROM=RUN_N22 RATE=1e9 GMAX=2 ./go.sh 22 "5.6" all 4 RUN_N22_r1e9
# Before the first run, put the supercomputer's settings in config.env (copy them from /home/app/lammps/run.sh).
set -eu
cd "$(dirname "$0")"
N=${1:-12}; SIG=${2:-"5.6 56"}; ST=${3:-all}; RUN=${5:-RUN_N$N}; SEED=${6:-1}; BASE=RUN_N$N
[ -n "${4:-}" ] && export NODES=$4      # 4th argument = number of nodes (works in any shell, no "NODES=4 ./go.sh" needed)
PY=${PYTHON:-python3}
# FROM: the new folder starts where the other one finished the compaction
if [ -n "${FROM:-}" ]; then
  FROM=${FROM%/}
  [ -d "$FROM" ] || { echo "FROM=$FROM: no such folder here"; exit 1; }
  [ "$FROM" != "$RUN" ] || { echo "FROM=$FROM is the run folder itself: give another 5th argument (e.g. ${FROM}_r1e9)"; exit 1; }
  [ -f "$FROM/02_dehydrate.done" ] || { echo "FROM=$FROM has not finished the dehydration (02) yet: wait for $FROM/02_dehydrate.done"; exit 1; }
  mkdir -p "$RUN"
  for f in gouge.data gouge.settings gouge.json lammps_seed config.env lmp.auto walltime 01_equil.done 01_equil.data 02_dehydrate.done 02_dehydrate.data dehydrate.dat dehyd_nw0.txt; do
    [ -f "$FROM/$f" ] && [ ! -f "$RUN/$f" ] && cp "$FROM/$f" "$RUN/"
  done
  for s in $SIG; do for f in 03_consolidate_$s.done 03_consolidate_$s.data 03_consolidate_$s.dat; do [ -f "$FROM/$f" ] && [ ! -f "$RUN/$f" ] && cp "$FROM/$f" "$RUN/"; done; done
  echo "$RUN: starts from the compacted gouge of $FROM (finished there: $(cd "$FROM" && ls *.done 2>/dev/null | tr '\n' ' '))"
fi
# random configurations to try (the most random one is kept): 300 for small systems; 30 for N > 20 (keeps the login node busy for minutes, not hours)
NCONF=${NCONF:-$([ "$N" -gt 20 ] && echo 30 || echo 300)}
$PY -c "import numpy, scipy" 2>/dev/null || { echo "installing numpy and scipy for $PY (pip --user)"; $PY -m pip install --user -q numpy scipy; }
# the figures of the shear (hot_spots.py, drawn by the jobs in RUN/figs): matplotlib; the runs do not need it
$PY -c "import matplotlib" 2>/dev/null || { echo "installing matplotlib for $PY (pip --user; for the figures in $RUN/figs)"; $PY -m pip install --user -q matplotlib || echo "note: matplotlib could not be installed: the jobs run, but draw no figures"; }
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
# options: one line each in $RUN/config.env (job.sh passes EXTRA to LAMMPS); given again, the line is replaced
set_opt() {
  local f=$RUN/config.env; [ -f "$f" ] || cp "$CE" "$f"
  { grep -v "# go.sh: $1\$" "$f" || true; echo "EXTRA=\"\${EXTRA:-} -var $1 $2\"   # go.sh: $1"; } > "$f.tmp" && mv -f "$f.tmp" "$f"
  echo "$RUN: LAMMPS -var $1 $2 for every job of this folder ($f)"
}
if [ -n "${DT:-}" ]; then set_opt dt "$DT"; fi
if [ -n "${DUMP:-}" ]; then set_opt dump_group "$DUMP"; fi
if [ -n "${MWMC:-}" ]; then set_opt final_mwmc "$MWMC"; fi
if [ -n "${FLOW:-}" ]; then set_opt flow "$FLOW"; fi
if [ -n "${RATE:-}" ]; then set_opt rate "$(awk -v r="$RATE" 'BEGIN{printf "%.6g", r/1e15}')"; fi   # 1/s -> 1/fs
if [ -n "${GMAX:-}" ]; then set_opt gmax "$GMAX"; fi
if [ -n "${NVS:-}" ]; then set_opt nvs "$NVS"; fi
if [ -n "${VFAC:-}" ]; then set_opt vfac "$VFAC"; fi
$S $RUN "$SIG" "$ST"
echo "submitted. watch with: qstat -u \$USER;  bash $RUN/status.sh $RUN   results in $RUN/ (shear_<sigma>.dat; python3 analysis/hot_spots.py $RUN)"
