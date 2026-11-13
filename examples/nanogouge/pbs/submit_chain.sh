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
# new inputs replace the old files (mv), not overwrite them: a LAMMPS running now keeps reading its own copy
if ls "$HERE"/in.* >/dev/null 2>&1; then INS="$HERE"; else INS="$TOP/lammps"; fi
for f in "$INS"/in.*; do cp "$f" "$RUNDIR/.new.$$" && mv -f "$RUNDIR/.new.$$" "$RUNDIR/$(basename "$f")"; done
cp "$HERE"/job.sh "$HERE"/run.sh "$RUNDIR"/
[ -f "$HERE/status.sh" ] && cp "$HERE/status.sh" "$RUNDIR"/
[ -f "$RUNDIR/config.env" ] || cp "$HERE/config.env" "$RUNDIR/"
[ -f "$HERE/autoconfig.sh" ] && bash "$HERE/autoconfig.sh" "$RUNDIR"
cd "$RUNDIR"; N_JOBS=10; QSUB_OPTS=""
[ -f ./config.auto ] && source ./config.auto
source ./config.env
# the wall time of the jobs (file "walltime" of the run folder, written by WALLTIME=3 ./go.sh): every walltime in the qsub options is replaced by it
set_walltime() { echo "$1" | sed -E 's/walltime=[0-9:]+,//g; s/,walltime=[0-9:]+//g; s/-l[[:space:]]*walltime=[0-9:]+//g' | awk -v w="$2" '{ if (NF) { $1 = $1; print $0 " -l walltime=" w } else print "-l walltime=" w }'; }
# WALLTIME=3 (hours) or 3:00:00: kept in this folder for every later job; the waiting jobs of this folder (old wall time) are replaced
NEWWALL=""
if [ -n "${WALLTIME:-}" ]; then
  W=$WALLTIME; case $W in *:*:*) ;; *:*) W="$W:00";; *) W="$W:00:00";; esac
  echo "$W" | grep -qE '^[0-9]+:[0-9]+:[0-9]+$' || { echo "WALLTIME=$WALLTIME: give hours (WALLTIME=3) or h:mm:ss (WALLTIME=3:00:00)"; exit 1; }
  OLDWALL=$(cat walltime 2>/dev/null || true); echo "$W" > walltime; NEWWALL=1; echo "wall time of the jobs of this folder: $W (kept in $RUNDIR/walltime)"
fi
[ -f walltime ] && QSUB_OPTS=$(set_walltime "$QSUB_OPTS" "$(cat walltime)")
# jobs of this user that work in this folder, in the given states (qselect -s: R running, QH waiting)
folder_jobs() {
  command -v qselect >/dev/null 2>&1 || return 0
  local j f
  for j in $(qselect -u "$USER" -s "$1" 2>/dev/null); do
    f=$(qstat -f "$j" 2>/dev/null | tr -d ' \t\n')
    case $f in *"SUBDIR_PHYS=$(pwd -P),"*|*"PBS_O_WORKDIR=$PWD,"*|*"PBS_O_WORKDIR=$(pwd -P),"*) echo "$j";; esac
  done
}
for f in gouge.data gouge.settings; do [ -f $f ] || { echo "missing $RUNDIR/$f (run build_gouge.py -o $RUNDIR/gouge)"; exit 1; }; done
printf 'SIGMAS="%s"\nSTAGES="%s"\n' "$SIGMAS" "$STAGES" > chain.env
[ -f chain.failed ] && { echo "removing chain.failed: $(cat chain.failed)"; rm -f chain.failed; }
# a chain of this folder that is still queued or running: the new jobs start after it
J=""; OLD=""
if [ -n "$NEWWALL" ]; then
  # new wall time: the waiting jobs of this folder still have the old one.  The new jobs start after the job working now; the waiting ones are deleted below.
  OLD=$(folder_jobs QH | tr '\n' ' ' | sed 's/ *$//')
  command -v qselect >/dev/null 2>&1 || echo "WARNING: no qselect: delete the waiting jobs of this folder by hand (qstat -u $USER; qdel <id>), they still have the old wall time"
elif [ -f chain.tail ]; then T=$(cat chain.tail); qstat "$T" >/dev/null 2>&1 && { J=$T; echo "the chain of this folder is still queued ($T): the new jobs follow it"; }; fi
# a job working in this folder now: the new jobs start after it (never two jobs in one folder)
if [ -z "$J" ]; then
  J=$(folder_jobs R | tail -1)
  [ -z "$J" ] && [ -f .chain.lock ] && H=$(cat .chain.lock) && [ -n "$H" ] && qstat "$H" 2>/dev/null | grep -q " R " && J=$H
  [ -n "$J" ] && echo "job $J is working in this folder: the new jobs start after it"
fi
for i in $(seq 1 ${N_JOBS}); do
  dep=""; [ -n "$J" ] && dep="-W depend=afterany:$J"
  if ! J=$(qsub -N ng_$i ${QSUB_OPTS} $dep -v "SUBDIR_PHYS=$(pwd -P),SUBDIR_LOGICAL=$PWD" job.sh 2>&1); then
    echo "ERROR: qsub refused job $i: $J"
    if [ -n "$NEWWALL" ] && [ "$i" = 1 ]; then
      # e.g. "Job exceeds queue resource limits": nothing was changed; the old wall time stays for this folder
      if [ -n "$OLDWALL" ]; then echo "$OLDWALL" > walltime; else rm -f walltime; fi
      echo "the wall time is not changed, the waiting jobs stay.  The longest wall time of the queues: qstat -Qf | grep -E 'Queue:|walltime'"
    fi
    exit 1
  fi
  echo "$J" > chain.tail; echo "job $i: $J"
done
[ -n "$OLD" ] && { qdel $OLD && echo "deleted the waiting jobs with the old wall time: $OLD"; }
echo "stages: $STAGES   sigmas: $SIGMAS   ($N_JOBS jobs; the last one submits more while work is left)   qsub options: ${QSUB_OPTS}"
