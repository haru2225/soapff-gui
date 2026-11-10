#!/bin/bash
# Usage:  pbs/submit_chain.sh RUNDIR "5.6 56" [stages]
#   RUNDIR  run folder; must already hold gouge.data and gouge.settings (from build_gouge.py)
#   "5.6 56" normal stresses (MPa) to shear at, each as its own branch after the dehydration
#   stages  default "all" = 01 02 then for each sigma: 03, 04 (N_SHEAR_JOBS times), 05.  Or e.g. "01 02" to stop earlier.
# Each job starts when the previous one has ended (afterany: a job stopped by the wall time does not block the chain; job.sh skips finished stages and
# continues stages 02 and 04 from their last restart).
set -eu
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$(dirname "$HERE")
RUNDIR=$1; SIGMAS=${2:-"5.6 56"}; STAGES=${3:-all}
mkdir -p "$RUNDIR"
if ls "$HERE"/in.* >/dev/null 2>&1; then cp "$HERE"/in.* "$RUNDIR"/; else cp "$TOP"/lammps/in.* "$RUNDIR"/; fi
cp "$HERE"/job.sh "$HERE"/run.sh "$RUNDIR"/
[ -f "$RUNDIR/config.env" ] || cp "$HERE/config.env" "$RUNDIR/"
cd "$RUNDIR"; source ./config.env
for f in gouge.data gouge.settings; do [ -f $f ] || { echo "missing $RUNDIR/$f (run build_gouge.py -o $RUNDIR/gouge)"; exit 1; }; done
Q=""; [ -n "${QUEUE}" ] && Q="-q ${QUEUE}"
sub() {   # sub NAME WALLTIME DEPENDENCY VARS
  local dep=""; [ -n "$3" ] && dep="-W depend=afterany:$3"
  qsub -N "$1" $Q ${RES_FLAG} "${SELECT}" ${RES_FLAG} "walltime=$2" $dep -v "$4" job.sh
}
want() { [ "$STAGES" = all ] || [[ " $STAGES " == *" $1 "* ]]; }
J=""
if want 01; then J=$(sub ng01 "$WALLTIME_SHORT" "" "STAGE=01"); echo "01: $J"; fi
if want 02; then for i in $(seq 1 ${N_DEHYD_JOBS}); do J=$(sub ng02_$i "$WALLTIME_LONG" "$J" "STAGE=02"); echo "02 part $i: $J"; done; fi
BASE=$J
for S in $SIGMAS; do
  J=$BASE
  if want 03; then J=$(sub ng03_$S "$WALLTIME_SHORT" "$J" "STAGE=03,SIGMA=$S"); echo "03 sigma=$S: $J"; fi
  if want 04; then for i in $(seq 1 ${N_SHEAR_JOBS}); do J=$(sub ng04_${S}_$i "$WALLTIME_LONG" "$J" "STAGE=04,SIGMA=$S"); echo "04 sigma=$S part $i: $J"; done; fi
  if want 05; then J=$(sub ng05_$S "$WALLTIME_SHORT" "$J" "STAGE=05,SIGMA=$S"); echo "05 sigma=$S: $J"; fi
done
