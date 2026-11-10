#!/bin/bash
# Usage: ./status.sh [RUNDIR]   what is done, what is running, the last lines of the current log (step, ps done, speed)
cd "${1:-$(ls -d RUN_N* 2>/dev/null | head -1)}" || exit 1
echo "== $PWD"; echo "-- qstat"; qstat -u "$USER" 2>/dev/null | grep -E "ng0|Job id|---" | head -20
echo "-- finished stages:"; ls *.done 2>/dev/null | tr '\n' ' '; echo
for L in $(ls -t log.0* 2>/dev/null | head -1); do
  echo "-- latest log: $L ($(( $(date +%s) - $(stat -c %Y "$L" 2>/dev/null || stat -f %m "$L") )) s since the last write)"
  grep -E "^ +[0-9]+ +[0-9.-]+ " "$L" | tail -2 | cut -c1-110
  grep -E "ERROR|Lost atoms|Bond atoms missing|Loop time" "$L" | tail -3
done
echo "-- job outputs (last line each):"; for f in $(ls -t ng0*.o* 2>/dev/null | head -3); do echo "$f: $(tail -1 $f | cut -c1-150)"; done
[ -f dehydrate.dat ] && { echo "-- dehydration (k n_water mw/mc rho_b Lz):"; tail -2 dehydrate.dat; }
for f in shear_*.dat; do [ -f "$f" ] && { echo "-- $f (step gamma tau sigma_n mu):"; tail -1 "$f" | cut -c1-100; }; done
