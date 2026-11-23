#!/bin/bash
# The time the jobs used, from the job history of PBS (qstat -fx: resources_used.walltime, the number of nodes).
#   bash pbs/usage.sh RUN_N54_fab_dis          the jobs of these run folders (their IDs from jobs.log), one line each and the sum
#   bash pbs/usage.sh RUN_*                    several folders
#   bash pbs/usage.sh --all                    every job of this user that PBS still has in its history
# PBS keeps finished jobs for a while only (job_history_duration, often 1-2 weeks): older jobs show as "not in the history".
# Node hours = wall time x nodes.  The time charged by the centre (its own command, see /etc/motd) can count differently.
set -u
hms() { echo "$1" | awk -F: '{ if (NF == 3) print $1 + $2 / 60 + $3 / 3600; else if (NF == 2) print $1 / 60 + $2 / 3600; else print 0 }'; }
one() {   # $1 = job ID, $2 = folder: prints  id  state  wall(h)  nodes  node-h  name
  local f st w n nm
  f=$(qstat -fx "$1" 2>/dev/null) || { printf '%-22s %-5s %9s %5s %9s  %s\n' "$1" - - - - "not in the history ($2)"; return; }
  st=$(echo "$f" | awk -F' = ' '/^ *job_state =/ {print $2; exit}')
  w=$(echo "$f" | awk -F' = ' '/^ *resources_used.walltime =/ {print $2; exit}')
  n=$(echo "$f" | awk -F' = ' '/^ *Resource_List.nodect =/ {print $2; exit}')
  [ -z "$n" ] && n=$(echo "$f" | tr -d '\n\t' | grep -o 'Resource_List.select = [0-9]*' | grep -o '[0-9]*$')
  nm=$(echo "$f" | awk -F' = ' '/^ *Job_Name =/ {print $2; exit}')
  w=$(hms "${w:-0}"); n=${n:-1}
  printf '%-22s %-5s %9.2f %5s %9.2f  %s\n' "$1" "${st:-?}" "$w" "$n" "$(awk -v a="$w" -v b="$n" 'BEGIN{print a * b}')" "${nm:-} $2"
}
[ $# -gt 0 ] || { sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
printf '%-22s %-5s %9s %5s %9s  %s\n' job state wall_h nodes node_h "name folder"
if [ "${1:-}" = --all ]; then
  command -v qselect >/dev/null 2>&1 || { echo "no qselect here: give the run folders instead"; exit 1; }
  for j in $(qselect -x -u "$USER" 2>/dev/null); do one "$j" ""; done
else
  for d in "$@"; do
    [ -f "$d/jobs.log" ] || { echo "$d: no jobs.log"; continue; }
    for j in $(sed -n 's/^===== job \([^ ]*\) on .*/\1/p' "$d/jobs.log" | awk '!seen[$0]++'); do one "$j" "$d"; done
  done
fi | awk '{ print } $3 != "-" { w += $3; nh += $5; n++ } $3 == "-" { m++ }
  END { printf "\nsum of %d job(s): %.1f h of wall time, %.1f node hours", n, w, nh; if (m) printf "  (%d job(s) no longer in the history: not counted)", m; print "" }'
