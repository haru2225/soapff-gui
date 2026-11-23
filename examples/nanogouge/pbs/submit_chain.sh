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
# the analysis of the shear (job.sh draws the figures with it after each job): analysis/ of the repository, or next to this script in the pack
for f in "$HERE/hot_spots.py" "$TOP/analysis/hot_spots.py"; do if [ -f "$f" ]; then cp "$f" "$RUNDIR/"; break; fi; done
for f in "$HERE/lin_figs.py" "$TOP/analysis/lin_figs.py"; do if [ -f "$f" ]; then cp "$f" "$RUNDIR/"; break; fi; done
[ -f "$HERE/status.sh" ] && cp "$HERE/status.sh" "$RUNDIR"/
[ -f "$RUNDIR/config.env" ] || cp "$HERE/config.env" "$RUNDIR/"
[ -f "$HERE/autoconfig.sh" ] && bash "$HERE/autoconfig.sh" "$RUNDIR"
cd "$RUNDIR"; N_JOBS=10; QSUB_OPTS=""
[ -f ./config.auto ] && source ./config.auto
source ./config.env
# GPU jobs (go.sh GPU=...): their own qsub options
GPUMODE=$(sed -n 's/^GPU="\{0,1\}\([1-9][0-9]*\).*/\1/p' ./config.env | tail -1)
[ -n "$GPUMODE" ] && { [ -n "${GPU_QSUB:-}" ] && QSUB_OPTS=$GPU_QSUB; echo "GPU jobs ($GPUMODE GPU(s) per node): qsub $QSUB_OPTS"; }
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
# OMP (go.sh OMP=2): ncpus of a chunk as asked, mpiprocs = ncpus / OMP, ompthreads = OMP
OMPN=$(sed -n 's/^OMP="\{0,1\}\([0-9][0-9]*\).*/\1/p' ./config.env | tail -1)
if [ -n "$OMPN" ] && [ "$OMPN" -gt 1 ] && [ -z "$GPUMODE" ]; then
  CH=$(echo "$QSUB_OPTS" | grep -o 'select=[0-9][^ ,]*' | head -1)
  NCPU=$(echo "$CH" | grep -o 'ncpus=[0-9]*' | head -1 | cut -d= -f2)
  [ -n "$NCPU" ] || { echo "OMP=$OMPN: no ncpus= in the qsub options ($QSUB_OPTS)"; exit 1; }
  [ $((NCPU % OMPN)) -eq 0 ] || { echo "OMP=$OMPN does not divide ncpus=$NCPU (cores per node): OMP=2, 4, 5, 8, 10 or 20 for 40"; exit 1; }
  NCH=$(echo "$CH" | sed -E 's/:mpiprocs=[0-9]+//; s/:ompthreads=[0-9]+//'):mpiprocs=$((NCPU / OMPN)):ompthreads=$OMPN
  QSUB_OPTS=$(echo "$QSUB_OPTS" | sed "s|$CH|$NCH|")
  echo "OMP=$OMPN: $((NCPU / OMPN)) MPI ranks x $OMPN OpenMP threads per node ($NCH)"
fi
# other qsub options than at the last submission (nodes, OMP, ...): the waiting jobs of this folder are replaced by new ones
NEWQ=""
if [ -f qsub.opts ]; then [ "$(cat qsub.opts)" != "$QSUB_OPTS" ] && NEWQ=1; elif [ -n "$OMPN" ] && [ "$OMPN" -gt 1 ]; then NEWQ=1; fi
# MAX_HOURS (go.sh MAXHOURS=72): from this submission on, at most that many hours of jobs (72 h of 24 h jobs: 3 jobs, no more after them;
# chain.left counts the jobs left; ./go.sh again gives MAX_HOURS more)
MAXH=$(sed -n 's/^MAX_HOURS="\{0,1\}\([0-9][0-9]*\).*/\1/p' ./config.env | tail -1)
MAXJ=$(sed -n 's/^MAX_JOBS="\{0,1\}\([0-9][0-9]*\).*/\1/p' ./config.env | tail -1)
if [ -n "$MAXJ" ] && [ "$MAXJ" -gt 0 ]; then
  N_JOBS=$MAXJ; echo "$N_JOBS" > chain.left
  echo "MAX_JOBS=$MAXJ: $N_JOBS job(s), then the chain stops (./go.sh again: $MAXJ more)"
elif [ -n "$MAXH" ] && [ "$MAXH" -gt 0 ]; then
  WJ=$(echo "$QSUB_OPTS" | grep -o 'walltime=[0-9:]*' | head -1 | cut -d= -f2); WJ=${WJ:-24:00:00}
  N_JOBS=$(echo "$WJ" | awk -F: -v m="$MAXH" '{h = $1 + $2 / 60 + $3 / 3600; n = int(m / h + 0.999); if (n < 1) n = 1; print n}')
  echo "$N_JOBS" > chain.left
  echo "MAX_HOURS=$MAXH: $N_JOBS job(s) of $WJ, then the chain stops (./go.sh again: $MAXH h more)"
else
  rm -f chain.left
fi
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
# a shear (04) begun with an older in.04_shear (no shear_<s>.v2: no velocity steps, no particle outputs) is done again with the new one;
# its files and those of the healing (05) after it go to old_04_<s>/ (nothing is deleted)
# a job of this folder (its SUBDIR_PHYS): .chain.lock / chain.tail copied from another folder (cp -r) name jobs of that one
job_here() { qstat -f "$1" 2>/dev/null | tr -d ' \t\n' | grep -q "SUBDIR_PHYS=$(pwd -P),"; }
for f in .chain.lock chain.tail; do
  if [ -f $f ]; then H=$(cat $f); if [ -n "$H" ] && qstat "$H" >/dev/null 2>&1 && ! job_here "$H"; then echo "$f names job $H of another folder (a copy?): removed"; rm -f $f; fi; fi
done
RUNNING=""
if [ -f .chain.lock ]; then H=$(cat .chain.lock); if [ -n "$H" ] && qstat "$H" 2>/dev/null | grep -q " R "; then RUNNING=$H; fi; fi
for S in $SIGMAS; do
  if [ -f "shear_$S.v2" ] || ! ls "shear_$S.done" "shear_$S.dat" shear_$S.restart.* >/dev/null 2>&1; then continue; fi
  if [ -n "$RUNNING" ]; then echo "note: the shear at $S MPa was begun with an older input, and job $RUNNING works in this folder now: run ./go.sh again when it has ended"; continue; fi
  D=old_04_$S; mkdir -p "$D"
  for f in shear_$S.* shear_${S}_* log.04_shear_$S heal_$S.* heal_${S}_* log.05_heal_$S particles_$S.* prof_$S.*; do if [ -e "$f" ]; then mv -f "$f" "$D/"; fi; done
  echo "the shear at $S MPa was begun with an older in.04_shear: its files (and those of the healing after it) are now in $RUNDIR/$D; it is done again with the new input"
done
[ -f chain.failed ] && { echo "removing chain.failed: $(cat chain.failed)"; rm -f chain.failed; }
# a chain of this folder that is still queued or running: the new jobs start after it
J=""; OLD=""
if [ -n "$NEWWALL" ] || [ -n "$NEWQ" ] || [ -f chain.left ]; then
  # new wall time (the waiting jobs of this folder still have the old one), or a limit of jobs (MAX_JOBS / MAX_HOURS: the waiting jobs are
  # replaced by the counted ones).  The new jobs start after the job working now; the waiting ones are deleted below.
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
[ -n "$OLD" ] && { qdel $OLD && echo "deleted the waiting jobs of this folder (replaced by the new ones): $OLD"; }
# the qsub options of these jobs: the jobs that the chain submits itself (job.sh) use the same
echo "$QSUB_OPTS" > qsub.opts
echo "stages: $STAGES   sigmas: $SIGMAS   ($N_JOBS jobs; $([ -f chain.left ] && echo "no more after them: MAX_JOBS / MAX_HOURS" || echo "the last one submits more while work is left"))   qsub options: ${QSUB_OPTS}"
