#!/bin/bash
# Stop every job of a run folder now (the running one and the waiting ones), e.g. to go on with other options:
#   bash pbs/stop.sh RUN_N54_fab_dis          then:  env OMP=4 ... ./go.sh ... RUN_N54_fab_dis   (the same command with the new option)
#   bash pbs/stop.sh RUN_N54_fab_dis -n       only list them
# The calculation goes on from its last save (all.restart + all.state: written at the end of every run and when a job's time is up);
# what the running job computed after that save is done again.
set -u
D=${1:?usage: bash pbs/stop.sh RUN_folder [-n]}; LIST=${2:-}
cd "$D" || exit 1
P=$(pwd -P)
command -v qselect >/dev/null 2>&1 || { echo "no qselect here (run it on the login node)"; exit 1; }
J=""
for j in $(qselect -u "$USER" 2>/dev/null); do
  f=$(qstat -f "$j" 2>/dev/null | tr -d ' \t\n')
  case $f in *"SUBDIR_PHYS=$P,"*|*"PBS_O_WORKDIR=$P,"*|*"PBS_O_WORKDIR=$PWD,"*)
    s=$(echo "$f" | grep -o 'job_state=[A-Z]' | cut -d= -f2); J="$J $j"; echo "job $j  state ${s:-?}";;
  esac
done
[ -n "$J" ] || { echo "no job of $D in the queue"; exit 0; }
if [ -f all.restart ]; then
  echo "last save (all.restart): $(date -r all.restart '+%Y-%m-%d %H:%M' 2>/dev/null || ls -l all.restart | awk '{print $6, $7, $8}')   now: $(date '+%Y-%m-%d %H:%M')"
  echo "  (what the running job computed after the last save is done again)"
fi
[ "$LIST" = -n ] && exit 0
qdel $J && echo "deleted:$J"
# wait until they are gone (a running job takes some seconds to end)
for i in $(seq 1 30); do
  left=""; for j in $J; do qstat "$j" >/dev/null 2>&1 && qstat "$j" 2>/dev/null | grep -qE " [RQHE] " && left="$left $j"; done
  [ -z "$left" ] && break; sleep 2
done
[ -n "${left:-}" ] && echo "still ending:$left (check with qstat -u $USER)"
rm -f .chain.lock
echo "stopped.  Go on: the same ./go.sh command (with the new options); it goes on from the last save"
