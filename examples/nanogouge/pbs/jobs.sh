#!/bin/bash
# Usage (in examples/nanogouge):
#   ./pbs/jobs.sh                     your PBS jobs: id, name, state (R running, Q/H waiting), wall time and the run folder each one works in
#   ./pbs/jobs.sh keep RUN_N67        delete every job of yours that does NOT work in RUN_N67 (other run folders, other jobs too); asks first
U=${USER:-$(whoami)}
MODE=${1:-list}; KEEP=${2:-}
if [ "$MODE" = keep ]; then
  [ -n "$KEEP" ] && [ -d "$KEEP" ] || { echo "usage: ./pbs/jobs.sh keep RUN_FOLDER   (no folder '$KEEP' here: $(ls -d RUN_* 2>/dev/null | tr '\n' ' '))"; exit 1; }
  KL=$(cd "$KEEP" && pwd); KP=$(cd "$KEEP" && pwd -P)
fi
command -v qselect >/dev/null 2>&1 || { echo "qselect not found (run this on the login node)"; exit 1; }
DEL=""; NK=0
printf '%-22s %-10s %-2s %-9s %s\n' "job" "name" "S" "walltime" "run folder"
for j in $(qselect -u "$U" 2>/dev/null); do
  # qstat -f folds long values onto lines that start with a tab: join them first
  I=$(qstat -f "$j" 2>/dev/null | awk '/^\t/{sub(/^\t/,""); last=last $0; next} {if (NR>1) print last; last=$0} END{print last}')
  [ -z "$I" ] && continue
  get() { echo "$I" | awk -F' = ' -v k="$1" '{sub(/^ +/,"",$1)} $1==k{print $2; exit}'; }
  NAME=$(get Job_Name); ST=$(get job_state); WT=$(get Resource_List.walltime)
  VL=$(get Variable_List | tr ',' '\n')
  var() { echo "$VL" | awk -F= -v k="$1" '$1==k{sub(/^[^=]*=/,""); print; exit}'; }
  DL=$(var SUBDIR_LOGICAL); DP=$(var SUBDIR_PHYS); DW=$(var PBS_O_WORKDIR)
  D=${DL:-$DW}
  printf '%-22s %-10s %-2s %-9s %s\n' "$j" "$NAME" "$ST" "$WT" "${D:-?}"
  if [ "$MODE" = keep ]; then
    case "|$DL|$DP|$DW|" in
      *"|$KL|"*|*"|$KP|"*) NK=$((NK+1)) ;;
      *) DEL="$DEL $j" ;;
    esac
  fi
done
[ "$MODE" = keep ] || exit 0
echo
echo "kept: $NK job(s) of $KL"
[ -z "$DEL" ] && { echo "nothing else to delete"; exit 0; }
echo "to delete (all your other jobs, running ones too):$DEL"
printf 'delete them? [y/N] '; read -r A
case "$A" in y|Y|yes) qdel $DEL && echo "deleted:$DEL" ;; *) echo "nothing deleted" ;; esac
