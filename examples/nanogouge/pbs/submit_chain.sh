#!/bin/bash
# Usage:  pbs/submit_chain.sh RUNDIR "5.6 56" [stages]
#   RUNDIR  run folder; must already hold gouge.data and gouge.settings (from build_gouge.py)
#   "5.6 56" normal stresses (MPa) to shear at, one after the other after the dehydration
#   stages  default "all" = 01 02, then 03 04 05 for each sigma.  Or e.g. "01 02" to stop earlier.
# Submits a chain of N_JOBS identical jobs (qsub -W depend=afterany).  Each job works on the first unfinished stage, continues it from
# its last restart and stops LAMMPS cleanly before the wall time (any wall time, 1 h too); the last job of the chain submits the next one
# while work is left.  Running this again is safe: finished stages are skipped, and new jobs wait for a chain that is still queued.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$(dirname "$HERE")
RUNDIR=$1; SIGMAS=${2:-"5.6 56"}; STAGES=${3:-all}
mkdir -p "$RUNDIR"
if ls "$HERE"/in.* >/dev/null 2>&1; then cp "$HERE"/in.* "$RUNDIR"/; else cp "$TOP"/lammps/in.* "$RUNDIR"/; fi
cp "$HERE"/job.sh "$HERE"/run.sh "$RUNDIR"/
[ -f "$HERE/status.sh" ] && cp "$HERE/status.sh" "$RUNDIR"/
[ -f "$RUNDIR/config.env" ] || cp "$HERE/config.env" "$RUNDIR/"
[ -f "$HERE/autoconfig.sh" ] && bash "$HERE/autoconfig.sh" "$RUNDIR"
cd "$RUNDIR"; N_JOBS=10; QSUB_OPTS=""
[ -f ./config.auto ] && source ./config.auto
source ./config.env
for f in gouge.data gouge.settings; do [ -f $f ] || { echo "missing $RUNDIR/$f (run build_gouge.py -o $RUNDIR/gouge)"; exit 1; }; done
printf 'SIGMAS="%s"\nSTAGES="%s"\n' "$SIGMAS" "$STAGES" > chain.env
[ -f chain.failed ] && { echo "removing chain.failed: $(cat chain.failed)"; rm -f chain.failed; }
# a chain of this folder that is still queued or running: the new jobs start after it
J=""
if [ -f chain.tail ]; then T=$(cat chain.tail); qstat "$T" >/dev/null 2>&1 && { J=$T; echo "the chain of this folder is still queued ($T): the new jobs follow it"; }; fi
for i in $(seq 1 ${N_JOBS}); do
  dep=""; [ -n "$J" ] && dep="-W depend=afterany:$J"
  J=$(qsub -N ng_$i ${QSUB_OPTS} $dep -v "SUBDIR_PHYS=$(pwd -P),SUBDIR_LOGICAL=$PWD" job.sh)
  echo "$J" > chain.tail; echo "job $i: $J"
done
echo "stages: $STAGES   sigmas: $SIGMAS   ($N_JOBS jobs; the last one submits more while work is left)"
