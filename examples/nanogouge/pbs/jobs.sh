#!/bin/bash
# Usage (in examples/nanogouge):
#   ./pbs/jobs.sh                     your PBS jobs: id, name, state (R running, Q/H waiting), wall time and the run folder each one works in
#   ./pbs/jobs.sh keep RUN_N67        keep only the jobs of RUN_N67 with its wall time (RUN_N67/walltime, set by WALLTIME=... ./go.sh);
#                                     every other job of yours is deleted: other folders, other jobs, and RUN_N67 jobs with an old wall time
#                                     (running ones too).  It shows the list and asks first.
U=${USER:-$(whoami)}
MODE=${1:-list}; KEEP=${2:-}
tosec() { echo "$1" | awk -F: '{s=0; for(i=1;i<=NF;i++) s=s*60+$i; print int(s)}'; }
if [ "$MODE" = keep ]; then
  [ -n "$KEEP" ] && [ -d "$KEEP" ] || { echo "usage: ./pbs/jobs.sh keep RUN_FOLDER   (no folder '$KEEP' here: $(ls -d RUN_* 2>/dev/null | tr '\n' ' '))"; exit 1; }
  KL=$(cd "$KEEP" && pwd); KP=$(cd "$KEEP" && pwd -P)
  KW=""; [ -f "$KEEP/walltime" ] && KW=$(cat "$KEEP/walltime")
  echo "keeping the jobs of $KL${KW:+ with wall time $KW}"
fi
# job ids: qselect, or the numbers of the list of qstat -u (enough for qstat / qdel)
IDS=$(qselect -u "$U" 2>/dev/null)
[ -z "$IDS" ] && IDS=$(qstat -u "$U" 2>/dev/null | awk '$1 ~ /^[0-9]/{id=$1; sub(/\..*/,"",id); sub(/\*$/,"",id); print id}')
[ -z "$IDS" ] && { echo "no jobs of $U (qstat -u $U shows nothing)"; exit 0; }
DEL=""; NK=0
printf '%-24s %-10s %-2s %-9s %-10s %s\n' "job" "name" "S" "walltime" "" "run folder"
for j in $IDS; do
  # qstat -f folds long values onto lines that start with a tab: join them first
  I=$(qstat -f "$j" 2>/dev/null | awk '/^\t/{sub(/^\t/,""); last=last $0; next} {if (NR>1) print last; last=$0} END{print last}')
  [ -z "$I" ] && continue
  get() { echo "$I" | awk -F' = ' -v k="$1" '{sub(/^ +/,"",$1)} $1==k{print $2; exit}'; }
  NAME=$(get Job_Name); ST=$(get job_state); WT=$(get Resource_List.walltime)
  VL=$(get Variable_List | tr ',' '\n')
  var() { echo "$VL" | awk -F= -v k="$1" '$1==k{sub(/^[^=]*=/,""); print; exit}'; }
  DL=$(var SUBDIR_LOGICAL); DP=$(var SUBDIR_PHYS); DW=$(var PBS_O_WORKDIR)
  D=${DL:-$DW}
  W=""
  if [ "$MODE" = keep ]; then
    case "|$DL|$DP|$DW|" in
      *"|$KL|"*|*"|$KP|"*)
        if [ -n "$KW" ] && [ "$(tosec "$WT")" != "$(tosec "$KW")" ]; then W="DELETE(old)"; DEL="$DEL $j"; else W="keep"; NK=$((NK+1)); fi ;;
      *) W="DELETE"; DEL="$DEL $j" ;;
    esac
  fi
  printf '%-24s %-10s %-2s %-9s %-10s %s\n' "$j" "$NAME" "$ST" "$WT" "$W" "${D:-?}"
done
[ "$MODE" = keep ] || exit 0
echo
echo "keep: $NK job(s)"
[ -z "$DEL" ] && { echo "nothing to delete"; exit 0; }
echo "delete:$DEL"
printf 'delete them? [y/N] '; read -r A
case "$A" in
  y|Y|yes)
    qdel $DEL; echo "deleted:$DEL"
    echo "now:"; sleep 3; qstat -u "$U" 2>/dev/null | tail -n +4
    echo "(the kept jobs that waited for a deleted one start by themselves; if they stay in H with nothing running, run the same ./go.sh again)" ;;
  *) echo "nothing deleted" ;;
esac
