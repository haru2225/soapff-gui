#!/bin/bash
# Usage: ./status.sh [RUNDIR]   what is done, what is running, speed and the computing time left of the current stage, MPI ranks per node
cd "${1:-$(ls -d RUN_N* 2>/dev/null | head -1)}" || exit 1
echo "== $PWD"; echo "-- qstat"; qstat -u "$USER" 2>/dev/null | grep -E "ng|Job id|---" | head -20
echo "-- finished stages: $(ls *.done 2>/dev/null | tr '\n' ' ')"
[ -f chain.env ] && echo "-- chain: $(tr '\n' ' ' < chain.env)  last job: $(cat chain.tail 2>/dev/null)"
[ -f chain.failed ] && echo "-- FAILED: $(cat chain.failed)   (fix it, then ./go.sh again)"
L=$(ls -t log.0* 2>/dev/null | head -1)
if [ -n "$L" ]; then
  echo "-- latest log: $L ($(( $(date +%s) - $(stat -c %Y "$L" 2>/dev/null || stat -f %m "$L") )) s since the last write)"
  TH=$(grep -E "^ +[0-9]+ +[0-9.-]+ " "$L" | tail -1)
  echo "$TH" | cut -c1-110
  grep -E "ERROR|Lost atoms|Bond atoms missing" "$L" | tail -3
  grep -E "Loop time|Performance:|CPU use" "$L" | tail -3 | cut -c1-110      # CPU use ~100% per MPI task = every rank has its own core
  # speed (last thermo column = time steps per second) and the steps left of this stage (default lengths; EXTRA changes them)
  echo "$TH" | awk -v L="$L" -v D="$(tail -1 dehydrate.dat 2>/dev/null)" '{
    step=$1; sp=$NF; if (sp+0 <= 0) exit
    tgt=0
    if (L ~ /01_equil/) tgt=1100000
    else if (L ~ /03_consolidate/) tgt=2000000
    else if (L ~ /04_shear/) tgt=30000000
    else if (L ~ /05_heal/) tgt=5000000
    else if (L ~ /02_dehyd/) { split(D, d, " "); r=d[3]+0; if (r <= 0) r=3.0
      cyc = (r > 1 ? (r-1)/(3*0.02) + 0.8/(3*0.01) : (r-0.2)/(3*0.01)); tgt = step + cyc*200000 }
    left = tgt - step; if (left < 0) left = 0
    printf "-- speed %.2f steps/s = %.3f ns/day (1 fs).  This stage: step %d of %d, about %.1f h of computing left (= jobs of 1 h minus ~8 min each)\n", sp, sp*86400/1e6, step, tgt, left/sp/3600 }'
fi
echo "-- job outputs (last line each):"; for f in $(ls -t ng*.o* 2>/dev/null | head -3); do echo "$f: $(tail -1 $f | cut -c1-150)"; done
[ -f dehydrate.dat ] && { echo "-- dehydration (k n_water mw/mc rho_b Lz):"; tail -2 dehydrate.dat; }
for f in shear_*.dat; do [ -f "$f" ] && { echo "-- $f (step gamma tau sigma_n mu):"; tail -1 "$f" | cut -c1-100; }; done
if [ -f procs.txt ]; then
  echo "-- MPI ranks per node (procs.txt of the latest job; the same number on every node = all nodes are used):"
  awk -F': ' '$1 ~ /^[0-9]/{print $3}' procs.txt | sort | uniq -c | awk '{printf "  %s: %s ranks\n", $2, $1}'
fi
