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
#                 01 minimisation, NVT 0.1 ns, NPzT 0.1 ns;  02 2 ns (removals of 5 % of the initial water down to m_w/m_c 0.2; below
#                 m_w/m_c = 1 twice the time per removal: the form of Underwood & Bourg 2020, who took 134 ns);  then at each normal stress:
#                 03 NPzT in pieces of 50 ps until Lz changes by < 0.1 % over a piece (0.1-1 ns);  04 1e9 /s to the strain 3 (3 ns), then
#                 3.2e9 and 1e10 /s, 0.5 of strain each;  05 hold 50 ps.  About 11 ns in all.  NRUN=2000: a test, every run 2000 steps.
#                 The shear of a gouge of stacks (the foliation forming):  env STACK=3 MWMC0=1.5 WALLTIME=12 ./go.sh 22 "5.6 56" one 4 RUN_N22_fol
#                 DRYNS=0: no dehydration, the shear at the water content of the structure (MWMC0):
#                 env STACK=3 MWMC0=1.5 DRYNS=0 WALLTIME=12 ./go.sh 22 "5.6 56" one 4 RUN_N22_fol
#                 Lengths (one): NEQ=100000 (01), DRYNS=2 (02 in ns), DFRAC=0.05 (02 per removal), NCONS=50000 (03 piece), CONSTOL=0.001, NCONSMAX=20,
#                 RATE=1e9, GMAX=1, NVS=2, GSTEP=0.5, HOLD=50000 (05)
#   STACK=3       the particles in stacks (tactoids) of 3, face to face at 15.5 A (2 water layers; D001=12.5: 1 layer), the stacks at random
#                 (default 1: single particles at random, U&B 2020: they stack only during a slow dehydration, ~100 ns);
#                 its own structure in RUN_N<N>_stack3 (copied to the run folder)
#   MWMC0=1.5     water content of the new structure (default 3, U&B 2020): less water, fewer atoms (N=22: 0.27 M instead of 0.49 M), a short 02
#   Dense gouges at a set water content, two fabrics (as Sakuma et al. 2023: nominally randomly oriented NRO, preferentially oriented PO):
#   MWMC0=0, 0.1, 0.2 (dry, 10 %, 20 % = m_w/m_c), below 1.5: a loose structure (the box of m_w/m_c 1.5 (po) or 2 (nro), BOXW) with the water on the clay
#                 surfaces and in the interlayers first, compacted slowly in 01 (the box shrinks at CRATE = 1e9 /s, fix deform, in pieces of 10 ps, until the pressure reaches the first normal
#                 stress; then 0.3 ns NPT at it; 03 in pieces of 20 ps), no dehydration
#   ORIENT=nro    the particles (single, or stacks with STACK) at random orientation in a cube, compacted the same in x, y, z (01 NPT iso):
#                 nominally random (default)
#   ORIENT=po     the particles with their layers within 10 deg of the shear plane (TILT=10) in a tall box Lx = Ly = LX (default ~ the
#                 compacted height), compacted along z (01 NPzT): preferentially oriented.  E.g. (N=54 single particles, the compacted gouge
#                 about 2 particle diameters each way; dry 0.12 M atoms, 20 %: 0.20 M; NOSUB=1 first to look at RUN/gouge.xyz):
#                 env MWMC0=0   ORIENT=nro WALLTIME=12 ./go.sh 54 "5.6 56" one 4 RUN_N54_w0_nro
#                 env MWMC0=0.1 ORIENT=po  WALLTIME=12 ./go.sh 54 "5.6 56" one 4 RUN_N54_w0.1_po
#   Fabrics made by a dehydration, then the water of choice (packmol) and the shear (two steps):
#   FABRIC=dis    1) a fabric: the particles at random (MWMC0=2, a cube), dehydrated the same in x y z (NPT iso) as Underwood & Bourg
#                 2020: removals of 1 % of the initial water (DFRAC=0.01), twice the time per removal below m_w/m_c 1 (where the particles
#                 meet and stack), DRYNS=100 ns in all (about their rate, 1 % per ns; they took 134 ns); on the way the first state at each m_w/m_c of
#                 WSNAP="0.5 0.3 0.2 0.1 0" is written: 02_w<value>.data (snapshots.dat lists them); it ends at the last one
#                 (02_dehydrate.data, the trajectory 02_dehydrate.dump: analysis/stacks.py follows the stacking and the order)
#                 env FABRIC=dis WALLTIME=24 ./go.sh 54 "5.6" one 4 RUN_N54_fab_dis
#   FABRIC=ori    the same from particles placed as sedimented (layers within TILT=30 deg of the horizontal, a tall box, MWMC0=1.5),
#                 dehydrated along z only (NPzT, x y fixed: a uniaxial dehydration, the particles align and stack)
#                 env FABRIC=ori WALLTIME=24 ./go.sh 54 "5.6" one 4 RUN_N54_fab_ori
#   FROMFAB=RUN_N54_fab_dis MWMC0=0.2    2) the shear of the snapshot 02_w0.2.data (its own water); it can start while the dehydration
#                 goes on.  MWMC0 without a snapshot (or REWATER=packmol / grid), when the dehydration is done: rewater.py takes out all the
#                 water of 02_dehydrate.data and puts MWMC0 into the space between the particles (packmol: ./bin/packmol made by ./get_packmol.sh, called by go.sh once;
#                 no packmol: a grid, the sites nearest to the clay first; REWATER=grid chooses that); then in.all: minimisation, NVT
#                 0.1 ns, the compaction at the first normal stress (ori: along z, dis: x y z) and the shear as for the dense gouges
#                 env FROMFAB=RUN_N54_fab_dis MWMC0=0.2 WALLTIME=24 ./go.sh 54 "5.6 56" one 4 RUN_N54_dis_w0.2
#   GPU=1         the jobs on a GPU node (GPUs per node), with a LAMMPS of the GPU or KOKKOS package (found by pbs/job.sh, or GPU_LMP=/path/lmp);
#                 GPU_QSUB="-q <queue> -l select=1:ncpus=..:mpiprocs=..:ngpus=1" (default: the line of pbs/gpu_check.sh), GPU_MODULES="module
#                 load cuda", GPU_PKG=gpu|kokkos, GPU_RPG=8 (ranks per GPU).  The 4th argument (nodes) does not count then; GPU=0: the CPU again.
#                 env GPU=1 FABRIC=dis WALLTIME=24 ./go.sh 54 "5.6" one 1 RUN_N54_fab_dis
#   MAXJOBS=5     at most 5 jobs of WALLTIME from this submission on (a job running now not counted), then the chain stops; the same command
#                 again: 5 more (MAXJOBS=0: no limit).  MAXHOURS=72: the same in hours (WALLTIME=24: 3 jobs).  Without them the chain goes on
#                 until the calculation is done.
#   NOSUB=1       build the structure only (RUN/gouge.data and RUN/gouge.xyz to look at in OVITO), submit nothing; the same command
#                 without NOSUB then submits it
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
MWMC0_GIVEN=${MWMC0:-}
FABRIC=${FABRIC:-}; FROMFAB=${FROMFAB:-}; FFAB=""
case "$FABRIC" in
  "") ;;
  ori) ORIENT=po; TILT=${TILT:-30}; MWMC0=${MWMC0:-1.5};;
  dis) ORIENT=nro; TILT=${TILT:-90}; MWMC0=${MWMC0:-2};;
  *) echo "FABRIC=$FABRIC: ori (aligned and stacked by a uniaxial dehydration) or dis (dispersed: an isotropic one)"; exit 1;;
esac
if [ -n "$FABRIC$FROMFAB" ] && [ "$ST" != one ]; then echo "FABRIC / FROMFAB: with the stages \"one\" (3rd argument)"; exit 1; fi
if [ -n "$FROMFAB" ]; then
  FROMFAB=${FROMFAB%/}
  [ -z "$FABRIC" ] || { echo "FABRIC and FROMFAB: one or the other"; exit 1; }
  [ -n "$MWMC0_GIVEN" ] || { echo "FROMFAB=$FROMFAB: give the water content of the new structure, e.g. MWMC0=0, 0.1 or 0.2 (m_w/m_c)"; exit 1; }
  FFAB=$(sed -n 's/^FABRIC=//p' "$FROMFAB/fabric.env" 2>/dev/null | head -1)
  case "$FFAB" in ori) ORIENT=po;; dis) ORIENT=nro;; *) FFAB="";; esac
fi
MWMC0=${MWMC0:-3}
awk -v w="$MWMC0" 'BEGIN{ exit !(w + 0 >= 0 && w + 0 <= 5) }' || { echo "MWMC0=$MWMC0: the water content of the new structure, 0 <= m_w/m_c <= 5"; exit 1; }
ORIENT=${ORIENT:-nro}; case "$ORIENT" in nro) TILT=${TILT:-90};; po) TILT=${TILT:-10};; *) echo "ORIENT=$ORIENT: nro or po"; exit 1;; esac
# below m_w/m_c 1.5 the particles do not fit at random: a loose structure (the box of BOXW: po 1.5, nro 2 = N=54 stacks at random in a cube),
# compacted by LAMMPS in 01, no dehydration
BOXW=${BOXW:-$([ "$ORIENT" = po ] && echo 1.5 || echo 2)}; LOOSE=$(awk -v w="$MWMC0" -v b="$BOXW" 'BEGIN{print (w + 0 < b + 0) ? 1 : 0}')
# a fabric is built full of water (no loose box); the shear of a fabric (FROMFAB) is compacted and sheared as a dense gouge
[ -n "$FABRIC" ] && LOOSE=0
[ -n "$FROMFAB" ] && LOOSE=1
[ -z "${D001:-}" ] && [ "$(awk -v w="$MWMC0" 'BEGIN{print (w + 0 < 0.15)}')" = 1 ] && D001=12.5     # dry, 10 %: 1 water layer (or none)
SHAPE=""
if [ "$LOOSE" = 1 ] || [ -n "$FABRIC" ]; then
  # nro: a cube (compacted iso, stays a cube); po: Lx = Ly about the compacted height (clay 24,630 A^3 per particle, porosity ~0.4)
  if [ "$ORIENT" = po ]; then SHAPE="--lx ${LX:-$(awk -v n="$N" 'BEGIN{printf "%.1f", (n * 41000) ^ (1 / 3)}')}"; else SHAPE="--aspect 1"; fi
fi
GB=$BASE; [ "$STACK" -gt 1 ] && GB=${GB}_stack$STACK; [ "$(awk -v w="$MWMC0" 'BEGIN{print (w + 0 == 3)}')" = 1 ] || GB=${GB}_w$MWMC0
[ "$ORIENT" = po ] && GB=${GB}_po
[ "$TILT" != "$([ "$ORIENT" = po ] && echo 10 || echo 90)" ] && GB=${GB}_t$TILT
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
py_ok() { tmo 120 "$@" -c "import numpy, matplotlib; from scipy.spatial import cKDTree" >/dev/null 2>&1 </dev/null; }
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
# a python with numpy to build a structure: ./sifpy first; else python3 of this node (numpy is enough: without a working scipy.spatial
# build_gouge.py uses its numpy cell list, the same structure, about 10 times slower)
pick_py() {
  get_py && return 0
  [ -n "${PYTHON:-}" ] && { echo "PYTHON=$PYTHON does not work"; return 1; }
  local c
  for c in "env PYTHONNOUSERSITE=1 python3" python3; do
    if tmo 120 $c -c "import numpy" >/dev/null 2>&1 </dev/null; then
      PY=$c
      if tmo 120 $c -c "from scipy.spatial import cKDTree" >/dev/null 2>&1 </dev/null; then echo "(no ./sifpy: using '$c' of this node, numpy and scipy)"
      else echo "(no ./sifpy: using '$c' of this node, numpy only: no scipy.spatial, build_gouge.py uses its numpy cell list, slower)"; fi
      return 0
    fi
  done
  echo "no python with numpy to build the structure: ./sifpy could not be made (see above; ./make_sif.sh needs internet once)"; return 1
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
# FROMFAB: the snapshot of the dehydration of a FABRIC run at m_w/m_c MWMC0 (02_w<MWMC0>.data, with its own water); none (or REWATER
# given): the clay of the fabric (02_dehydrate.data) with MWMC0 of water put in by rewater.py (packmol, else a grid)
WS0=$(awk -v w="$MWMC0" 'BEGIN{printf "%g", w}')
if [ -n "$FROMFAB" ] && [ ! -f "$RUN/gouge.data" ] && [ -z "${REWATER:-}" ] && [ -f "$FROMFAB/02_w$WS0.data" ]; then
  [ "$FROMFAB" != "$RUN" ] || { echo "FROMFAB=$FROMFAB is the run folder itself: give another 5th argument"; exit 1; }
  mkdir -p "$RUN"
  cp "$FROMFAB/02_w$WS0.data" "$RUN/gouge.data"; cp "$FROMFAB/gouge.settings" "$RUN/"
  { echo "FROMFAB=$FROMFAB"; echo "FABRIC=${FFAB:-unknown}"; echo "MWMC0=$MWMC0"; echo "SNAPSHOT=02_w$WS0.data"; } > "$RUN/fabric.env"
  echo "$RUN: the snapshot of the dehydration of $FROMFAB at m_w/m_c $WS0 (02_w$WS0.data, its own water: $(grep "^$WS0 " "$FROMFAB/snapshots.dat" 2>/dev/null | awk '{print "m_w/m_c " $2 ", " $3 " molecules"}'))"
fi
if [ -n "$FROMFAB" ] && [ ! -f "$RUN/gouge.data" ]; then
  [ "$FROMFAB" != "$RUN" ] || { echo "FROMFAB=$FROMFAB is the run folder itself: give another 5th argument"; exit 1; }
  if [ -z "${REWATER:-}" ] && [ ! -f "$FROMFAB/02_dehydrate.done" ]; then
    echo "FROMFAB=$FROMFAB: no snapshot at m_w/m_c $WS0 yet ($FROMFAB/02_w$WS0.data); there: $(cd "$FROMFAB" 2>/dev/null && ls 02_w*.data 2>/dev/null | tr '\n' ' ')"
    echo "  wait for it (the dehydration is going on), or REWATER=packmol to put the water into the end of the dehydration when it is done"; exit 1
  fi
  [ -f "$FROMFAB/02_dehydrate.done" ] && [ -f "$FROMFAB/02_dehydrate.data" ] || { echo "FROMFAB=$FROMFAB: the fabric is not ready (no $FROMFAB/02_dehydrate.done yet)"; exit 1; }
  [ -n "$FFAB" ] || echo "(no $FROMFAB/fabric.env: compacted as ORIENT=$ORIENT, $([ "$ORIENT" = po ] && echo "along z" || echo "the same in x y z"))"
  mkdir -p "$RUN"
  pick_py || exit 1
  if [ "${REWATER:-auto}" != grid ] && [ -z "${PACKMOL:-}" ] && [ ! -x bin/packmol ] && ! command -v packmol >/dev/null 2>&1 && [ ! -f .packmol_failed ]; then
    echo "making ./bin/packmol (once): ./get_packmol.sh"; ./get_packmol.sh || { touch .packmol_failed; echo "(get_packmol.sh failed, not tried again: the water goes on a grid; ./get_packmol.sh by hand)"; }
  fi
  $PY rewater.py "$FROMFAB/02_dehydrate.data" --mwmc "$MWMC0" --seed "$SEED" --method "${REWATER:-auto}" -o "$RUN/gouge.data" \
    || { rm -f "$RUN/gouge.data"; echo "rewater.py failed"; exit 1; }
  cp "$FROMFAB/gouge.settings" "$RUN/"
  { echo "FROMFAB=$FROMFAB"; echo "FABRIC=${FFAB:-unknown}"; echo "MWMC0=$MWMC0"; } > "$RUN/fabric.env"
  echo "$RUN: the fabric of $FROMFAB (${FFAB:-?}) with m_w/m_c $MWMC0 of water"
fi
# random configurations to try (the most random one is kept): 300 for small systems; 30 for N > 20 (keeps the login node busy for minutes, not hours)
NCONF=${NCONF:-$(if [ -n "$FABRIC" ]; then echo 2; elif [ "$STACK" -gt 1 ] || [ "$LOOSE" = 1 ]; then echo 10; elif [ "$N" -gt 20 ]; then echo 30; else echo 300; fi)}   # stacks: 10 (N=22: about 20 s each)
if [ ! -f "$RUN/gouge.data" ]; then
  mkdir -p "$RUN"
  if [ "$SEED" = 1 ] && [ "$RUN" != "$GB" ] && [ -f "$GB/gouge.data" ]; then
    # seed 1 builds the structure of RUN_N<N> again: copy it (minutes saved)
    cp "$GB/gouge.data" "$GB/gouge.settings" "$RUN/"; [ -f "$GB/gouge.json" ] && cp "$GB/gouge.json" "$RUN/"
    echo "$RUN: the same initial structure as $GB (copied)"
  else
    # numpy and scipy are needed only here (building a new structure)
    # ./sifpy (Singularity) first; else python3 of this node: build_gouge.py needs only numpy (without a working scipy.spatial it uses a
    # numpy cell list: the same structure, about 10 times slower)
    pick_py || exit 1
    echo "build_gouge.py with: $PY"
    # seed 1 in another folder (a test, a rate series): built in RUN_N<N> and copied, so all of them start from the same structure
    OUT=$RUN; if [ "$SEED" = 1 ] && [ "$RUN" != "$GB" ]; then OUT=$GB; mkdir -p "$GB"; fi
    $PY build_gouge.py -n $N --mwmc $MWMC0 --nconf $NCONF --seed $SEED --stack $STACK ${D001:+--d001 $D001} --tilt $TILT \
        $([ "$LOOSE" = 1 ] && echo "--box-mwmc $BOXW") $SHAPE -o "$OUT/gouge" || { echo "build_gouge.py failed"; exit 1; }
    if [ "$OUT" != "$RUN" ]; then cp "$GB/gouge.data" "$GB/gouge.settings" "$RUN/"; [ -f "$GB/gouge.json" ] && cp "$GB/gouge.json" "$RUN/"; echo "$RUN: structure built in $GB (for its real run too) and copied"; fi
  fi
  # another seed: LAMMPS uses it too (velocities, water removal, thermostat); job.sh reads lammps_seed
  [ "$SEED" != 1 ] && echo "$SEED" > "$RUN/lammps_seed"
fi
# the initial structure for OVITO: $RUN/gouge.xyz (extended XYZ; analysis/data2xyz.py, no numpy needed)
D2X=analysis/data2xyz.py; [ -f "$D2X" ] || D2X=data2xyz.py        # the pack: flat
if [ -f "$RUN/gouge.data" ] && [ ! -f "$RUN/gouge.xyz" ]; then
  { python3 $D2X "$RUN/gouge.data" 2>/dev/null </dev/null || ${PY:-false} $D2X "$RUN/gouge.data" </dev/null; } || echo "(no $RUN/gouge.xyz: python3 $D2X $RUN/gouge.data)"
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
if [ -n "${DRYNS:-}" ]; then set_opt tdry "$DRYNS"; fi
if [ -n "${NCONS:-}" ]; then set_opt ncons "$NCONS"; fi
if [ -n "${CONSTOL:-}" ]; then set_opt constol "$CONSTOL"; fi
if [ -n "${NCONSMAX:-}" ]; then set_opt nconsmax "$NCONSMAX"; fi
# a loose structure: compacted in 01 at the first normal stress (nro: iso, po: along z), the box damped (drag); no dehydration
if [ "$LOOSE" = 1 ] && [ "$ST" = one ]; then
  set_opt c01 "$([ "$ORIENT" = po ] && echo z || echo iso)"
  set_opt p01 "$(awk -v s="${SIG%% *}" 'BEGIN{printf "%.4g", s * 9.86923}')"
  set_opt pdrag 1.0
  # first a slow compaction at a set rate (fix deform, CRATE 1/s, default 1e9: nro ~0.35 ns, po ~1 ns) until the pressure reaches the
  # first normal stress, then 0.3 ns NPT at it
  [ -z "${NRUN:-}" ] && set_opt crate "$(awk -v r="${CRATE:-1e9}" 'BEGIN{printf "%.6g", r / 1e15}')"
  [ -z "${NRUN:-}" ] && set_opt npt01 300000
  # 03 in pieces of 20 ps (each a LAMMPS start: PPPM set for the box of the moment while it still shrinks), at most 1 ns
  [ -z "${NRUN:-}" ] && [ -z "${NCONS:-}" ] && set_opt ncons 20000 && set_opt nconsmax 50
  [ -z "${NRUN:-}" ] && [ -z "${CONSTOL:-}" ] && set_opt constol 0.0005
  [ -z "${DRYNS:-}" ] && set_opt tdry 0
fi
# a fabric: 01 at 1 atm (ori: NPzT, dis: NPT iso), 02 dehydration (ori: along z, dis: x y z) with snapshots, then the end
if [ -n "$FABRIC" ]; then
  set_opt c01 "$([ "$FABRIC" = ori ] && echo z || echo iso)"
  set_opt c02 "$([ "$FABRIC" = ori ] && echo z || echo iso)"
  # as U&B 2020: removals of 1 % of the initial water (DFRAC), twice the time per removal below m_w/m_c 1, DRYNS ns in all (default 100:
  # about their rate, 1 % per ns, 134 ns), down to the driest snapshot; snapshots WSNAP on the way (wet to dry, at most 8): 02_w<value>.data,
  # the initial structures of the shear runs (FROMFAB=<this folder> MWMC0=<value>)
  WS=$(for w in ${WSNAP:-0.5 0.3 0.2 0.1 0}; do echo "$w"; done | awk -v w0="$MWMC0" '$1 + 0 >= 0 && $1 + 0 < w0 + 0 {printf "%g\n", $1}' | sort -g -r -u | head -8 | tr '\n' ' ')
  [ -n "$WS" ] || { echo "WSNAP=${WSNAP:-}: no water content between 0 and MWMC0=$MWMC0"; exit 1; }
  i=0; for w in $WS; do i=$((i+1)); set_opt ws$i "$w" >/dev/null; done; set_opt nsnap $i
  echo "$RUN: snapshots at m_w/m_c $WS(02_w<value>.data)"
  [ -z "${MWMC:-}" ] && set_opt final_mwmc "$(echo $WS | awk '{print $NF}')"
  [ -z "${DFRAC:-}" ] && set_opt dfrac 0.01
  [ -z "${DRYNS:-}" ] && set_opt tdry 100
  set_opt stop02 1
  echo "FABRIC=$FABRIC" > "$RUN/fabric.env"
fi
# GPU=<GPUs per node>: the jobs on a GPU node with a LAMMPS of the GPU or KOKKOS package (pbs/job.sh finds it, or GPU_LMP=/path/lmp);
# GPU_QSUB: their qsub options (default: the line pbs/gpu_check.sh suggests), GPU_MODULES: modules for them (CUDA ...), GPU_PKG=gpu|kokkos,
# GPU_RPG: MPI ranks per GPU (GPU package, default 8).  Kept in $RUN/config.env; GPU=0: back to the CPU.
set_cfg() {
  { grep -v "# go.sh: $1\$" "$RUN/config.env" 2>/dev/null || true; if [ -n "$2" ]; then echo "$1=\"$2\"   # go.sh: $1"; fi; } > "$RUN/config.env.tmp"
  mv -f "$RUN/config.env.tmp" "$RUN/config.env"
}
if [ -n "${GPU:-}" ]; then
  case "$GPU" in
    0) for k in GPU GPU_QSUB GPU_LMP GPU_MODULES GPU_PKG GPU_RPG; do set_cfg $k ""; done; echo "$RUN: the jobs run on the CPU again";;
    [1-9]|[1-9][0-9])
      set_cfg GPU "$GPU"
      for k in GPU_LMP GPU_MODULES GPU_PKG GPU_RPG GPU_QSUB; do [ -n "${!k:-}" ] && set_cfg $k "${!k}"; done
      if ! grep -q "^GPU_QSUB=" "$RUN/config.env"; then
        GC=pbs/gpu_check.sh; [ -f "$GC" ] || GC=./gpu_check.sh
        GCO=$(bash "$GC" 2>&1); SUG=$(echo "$GCO" | sed -n 's/^SUGGEST: //p' | tail -1)
        if [ "${SUG#-q }" != "$SUG" ]; then
          set_cfg GPU_QSUB "$SUG"; echo "$RUN: GPU jobs with qsub $SUG (from $GC; another line: GPU_QSUB=\"...\")"
        else
          echo "$GCO"; echo "GPU=$GPU: the qsub options of the GPU jobs are not known; give them, e.g."
          echo "  env GPU=$GPU GPU_QSUB=\"-q <GPU queue> -l select=1:ncpus=<cores>:mpiprocs=<ranks>:ngpus=<GPUs>\" ... ./go.sh ..."; exit 1
        fi
      fi
      echo "$RUN: GPU jobs, $GPU GPU(s) per node ($RUN/config.env: $(grep '^GPU' "$RUN/config.env" | sed 's/   # go.sh.*//' | tr '\n' ' '))";;
    *) echo "GPU=$GPU: the number of GPUs per node (0: the CPU)"; exit 1;;
  esac
fi
# MAXJOBS=5 / MAXHOURS=72: at most that many jobs / hours of jobs from this submission on (not counting a job running now), then the chain
# stops; the same command again: as many more.  0: no limit again (default: none, the chain goes on until the calculation is done)
if [ -n "${MAXHOURS:-}" ]; then
  case "$MAXHOURS" in 0) set_cfg MAX_HOURS ""; echo "$RUN: no limit on the hours of the jobs";; [1-9]*) set_cfg MAX_HOURS "$MAXHOURS"; set_cfg MAX_JOBS "";; *) echo "MAXHOURS=$MAXHOURS: hours"; exit 1;; esac
fi
if [ -n "${MAXJOBS:-}" ]; then
  case "$MAXJOBS" in 0) set_cfg MAX_JOBS ""; echo "$RUN: no limit on the number of jobs";; [1-9]*) set_cfg MAX_JOBS "$MAXJOBS"; set_cfg MAX_HOURS "";; *) echo "MAXJOBS=$MAXJOBS: a number of jobs"; exit 1;; esac
fi
# WLIST (go_water.sh): the water contents m_w/m_c of in.water, wet to dry, at most 8
if [ -n "${WLIST:-}" ]; then
  i=0; for w in $WLIST; do i=$((i+1)); set_opt w$i "$w" >/dev/null; done
  set_opt nw $i
fi
if [ -n "${NOSUB:-}" ]; then
  echo "NOSUB: $RUN is ready (gouge.data, gouge.xyz for OVITO), nothing submitted.  The same command without NOSUB submits it."; exit 0
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
