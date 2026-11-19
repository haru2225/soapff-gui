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
#   NRUN=1000     a test of the whole chain: every run of every stage is 1000 steps (02: one water removal; 04: each of its 5 phases), outputs
#                 every 100 steps; in its own folder RUN_N<N>_test1000 (the structure of RUN_N<N> is copied if there).  Checks that all runs and
#                 all outputs work; the numbers mean nothing.   env NRUN=1000 WALLTIME=0:30:00 ./go.sh 22 "5.6 56" all 4
#   ./go.sh 22 "5.6 56" one 4   everything in ONE input file, lammps/in.all, production lengths (dt 1 fs; a job stopped by the wall time is
#                 continued by the next one, all.state), in its own folder RUN_N22_one (5th argument: another folder):
#                 01 minimisation, NVT 0.1 ns, NPzT 0.1 ns;  02 10 water removals, 0.2 ns NPzT after each (2 ns);  then at each normal stress:
#                 03 NPzT in pieces of 50 ps until Lz changes by < 0.1 % over a piece (0.1-1 ns);  04 1e9 /s to the strain 1 (1 ns), then
#                 3.2e9 and 1e10 /s, 0.5 of strain each;  05 hold 50 ps.  About 7 ns in all.  NRUN=2000: a test, every run 2000 steps.
#                 Lengths (one): NEQ=100000 (01), NDWELL=200000 (02, each removal), NCONS=50000 (03 piece), CONSTOL=0.001, NCONSMAX=20,
#                 RATE=1e9, GMAX=1, NVS=2, GSTEP=0.5, HOLD=50000 (05)
#   STACK=3       the particles in stacks (tactoids) of 3, face to face at 15.5 A (2 water layers; D001=12.5: 1 layer), the stacks at random
#                 (default 1: single particles at random, U&B 2020; within 2 ns of dehydration they hardly come together face to face);
#                 its own structure in RUN_N<N>_stack3 (copied to the run folder)
#   NFAC=2        "one": 02, 03, 04, 05 twice as long (04: twice the strain); also in a running folder: from its next run on
#   NVS=4         velocity steps after GMAX (default 4: the rate x3.16 each, 0.5 of strain each; 1e8 -> 1e10 /s; 0 = none), VFAC=3.16227766 their factor
#   The jobs draw the figures (analysis/hot_spots.py) in RUN_<..>/figs after each job of the shear; a shear begun with an older input is done again.
#   Python (build_gouge.py: numpy, scipy; the figures: matplotlib) = ./sifpy, python 3.11 in a Singularity image that make_sif.sh makes once
#   (go.sh calls it; no pip).  PYTHON=/path/python3 uses another one.
#   FROM=RUN_N22  start from the compacted gouge of another folder (its stages 01-03 are copied; only 03 of a missing sigma, 04, 05 run):
#                 a series of shear rates on the same gouge, e.g.
#                 env FROM=RUN_N22 RATE=1e9 GMAX=2 ./go.sh 22 "5.6" all 4 RUN_N22_r1e9
# Before the first run, put the supercomputer's settings in config.env (copy them from /home/app/lammps/run.sh).
set -eu
cd "$(dirname "$0")"
N=${1:-12}; SIG=${2:-"5.6 56"}; ST=${3:-all}; RUN=${5:-RUN_N$N}; SEED=${6:-1}; BASE=RUN_N$N
# the folder of the initial structure (seed 1): RUN_N<N>, or RUN_N<N>_stack<k> for stacks (STACK=k)
STACK=${STACK:-1}; case "$STACK" in [1-9]) ;; *) echo "STACK=$STACK: particles per stack, 1 to 9"; exit 1;; esac
GB=$BASE; [ "$STACK" -gt 1 ] && GB=${BASE}_stack$STACK
# NRUN (a short test of the whole chain) never goes into a running folder: its own folder RUN_N<N>_test<NRUN> unless a 5th argument is given
if [ -n "${NRUN:-}" ] && [ -z "${5:-}" ]; then RUN=RUN_N${N}_test$NRUN; fi
# stages "one": everything in one input file (lammps/in.all, every run 5000 steps), in its own folder RUN_N<N>_one
if [ "$ST" = one ] && [ -z "${5:-}" ]; then RUN=RUN_N${N}_one; [ "$STACK" -gt 1 ] && RUN=RUN_N${N}_one_stack$STACK; fi
# stages "water": a series of water contents (lammps/in.water; go_water.sh), in its own folder RUN_N<N>_water
if [ "$ST" = water ] && [ -z "${5:-}" ]; then RUN=RUN_N${N}_water; fi
[ -n "${4:-}" ] && export NODES=$4      # 4th argument = number of nodes (works in any shell, no "NODES=4 ./go.sh" needed)
# ---- one python for everything: ./sifpy (Singularity image of make_sif.sh: python 3.11, numpy, scipy, matplotlib), made once when missing;
#      PYTHON=... overrides it.  No pip.  A failed make_sif.sh (no internet?) is not tried again by go.sh (.sifpy_failed): run ./make_sif.sh by hand.
tmo() { local t=$1; shift; if command -v timeout >/dev/null 2>&1; then timeout "$t" "$@"; else "$@"; fi; }
py_ok() { tmo 120 "$@" -c "import numpy, scipy, matplotlib" >/dev/null 2>&1 </dev/null; }
PY=""
get_py() {
  [ -n "$PY" ] && return 0
  if [ -n "${PYTHON:-}" ]; then PY=$PYTHON; return 0; fi
  if ! { [ -x sifpy ] && py_ok ./sifpy; } && [ ! -f .sifpy_failed ]; then
    echo "making ./sifpy, the python of this folder (numpy, scipy, matplotlib in a Singularity image; once, a few minutes): ./make_sif.sh"
    ./make_sif.sh || { touch .sifpy_failed; echo "make_sif.sh failed (see above); go.sh does not try again: ./make_sif.sh by hand"; }
  fi
  if [ -x sifpy ] && py_ok ./sifpy; then PY=$PWD/sifpy; return 0; fi
  return 1
}
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
if [ ! -f "$RUN/gouge.data" ]; then
  mkdir -p "$RUN"
  if [ "$SEED" = 1 ] && [ "$RUN" != "$GB" ] && [ -f "$GB/gouge.data" ]; then
    # seed 1 builds the structure of RUN_N<N> again: copy it (minutes saved)
    cp "$GB/gouge.data" "$GB/gouge.settings" "$RUN/"; [ -f "$GB/gouge.json" ] && cp "$GB/gouge.json" "$RUN/"
    echo "$RUN: the same initial structure as $GB (copied)"
  else
    # numpy and scipy are needed only here (building a new structure)
    get_py || { [ -z "${PYTHON:-}" ] && python3 -c "import numpy, scipy" 2>/dev/null </dev/null && PY=python3 && echo "(using python3 of this node, which has numpy and scipy)"; } \
      || { echo "no python with numpy and scipy to build the structure: ./sifpy could not be made (see above; ./make_sif.sh needs internet once)"; exit 1; }
    # seed 1 in another folder (a test, a rate series): built in RUN_N<N> and copied, so all of them start from the same structure
    OUT=$RUN; if [ "$SEED" = 1 ] && [ "$RUN" != "$GB" ]; then OUT=$GB; mkdir -p "$GB"; fi
    $PY build_gouge.py -n $N --mwmc 3.0 --nconf $NCONF --seed $SEED --stack $STACK ${D001:+--d001 $D001} -o "$OUT/gouge" || { echo "build_gouge.py failed"; exit 1; }
    if [ "$OUT" != "$RUN" ]; then cp "$GB/gouge.data" "$GB/gouge.settings" "$RUN/"; [ -f "$GB/gouge.json" ] && cp "$GB/gouge.json" "$RUN/"; echo "$RUN: structure built in $GB (for its real run too) and copied"; fi
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
if [ -n "${NRUN:-}" ]; then set_opt nrun "$NRUN"; fi
if [ -n "${NFAC:-}" ]; then set_opt nfac "$NFAC"; fi
if [ -n "${VFAC:-}" ]; then set_opt vfac "$VFAC"; fi
if [ -n "${HEAL:-}" ]; then set_opt heal "$HEAL"; fi
if [ -n "${GSTEP:-}" ]; then set_opt gstep "$GSTEP"; fi
if [ -n "${HOLD:-}" ]; then set_opt nhold "$HOLD"; fi
if [ -n "${DFRAC:-}" ]; then set_opt dfrac "$DFRAC"; fi
if [ -n "${NEQ:-}" ]; then set_opt neq "$NEQ"; fi
if [ -n "${NDWELL:-}" ]; then set_opt ndwell "$NDWELL"; fi
if [ -n "${NCONS:-}" ]; then set_opt ncons "$NCONS"; fi
if [ -n "${CONSTOL:-}" ]; then set_opt constol "$CONSTOL"; fi
if [ -n "${NCONSMAX:-}" ]; then set_opt nconsmax "$NCONSMAX"; fi
# WLIST (go_water.sh): the water contents m_w/m_c of in.water, wet to dry, at most 8
if [ -n "${WLIST:-}" ]; then
  i=0; for w in $WLIST; do i=$((i+1)); set_opt w$i "$w" >/dev/null; done
  set_opt nw $i
fi
$S $RUN "$SIG" "$ST"
echo "submitted. watch with: qstat -u \$USER;  bash $RUN/status.sh $RUN   results in $RUN/ (shear_<sigma>.dat; figures in $RUN/figs)"

# ---- the figures (hot_spots.py, drawn by the jobs in $RUN/figs) need numpy and matplotlib: ./sifpy (made now if missing; the jobs are already
#      submitted), kept in $RUN/config.env (PYTHON=...) for the jobs.  None: the jobs run without figures.
if get_py && tmo 120 $PY -c "import numpy, matplotlib" >/dev/null 2>&1 </dev/null; then
  { grep -v "# go.sh: python" "$RUN/config.env" || true; echo "PYTHON=\"$PY\"   # go.sh: python"; } > "$RUN/config.env.tmp" && mv -f "$RUN/config.env.tmp" "$RUN/config.env"
  echo "figures: drawn by the jobs with $PY in $RUN/figs"
else
  echo "note: no python with matplotlib (./sifpy could not be made, see above): the jobs run without figures.  ./figs.sh $RUN (finds or makes a python, draws them; then the jobs draw them too)"
fi
