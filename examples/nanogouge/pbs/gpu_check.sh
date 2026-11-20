#!/bin/bash
# What GPUs this PBS system has: the queues, the nodes with GPUs, the LAMMPS builds with the GPU or KOKKOS package, the CUDA modules.
# The last line, "SUGGEST: <qsub options>", is a resource line for one GPU node (go.sh GPU=1 takes it when GPU_QSUB is not given).
#   bash pbs/gpu_check.sh        (on the login node)
APP=${APP_ROOT:-/home/app}
echo "-- queues (qstat -Q)"; qstat -Q 2>/dev/null | head -30
echo "-- nodes with GPUs (pbsnodes -a): node, ngpus, ncpus, queue"
NODES=$(pbsnodes -a 2>/dev/null | awk '
  /^[^ \t]/ { if (g + 0 > 0) print n, g, c, q; n = $1; g = 0; c = 0; q = "" }
  $1 == "resources_available.ngpus" { g = $3 } $1 == "resources_available.ncpus" { c = $3 }
  $1 == "queue" || $1 == "resources_available.Qlist" { q = $3 }
  END { if (g + 0 > 0) print n, g, c, q }')
if [ -n "$NODES" ]; then echo "$NODES" | sort -k2,2nr -k3,3nr | awk '{c[$2" "$3" "$4]++} END{for (k in c) print "  " c[k] " nodes: ngpus ncpus queue = " k}'
else echo "  none listed (no pbsnodes access, or no ngpus resource)"; fi
echo "-- LAMMPS builds with the GPU or KOKKOS package"
for f in $(ls -d "$APP"/lammps*/bin/lmp* "$APP"/lammps*/*/bin/lmp* "$APP"/lammps*/*/*/bin/lmp* "$APP"/lammps*/*/*/*/bin/lmp* "$APP"/lammps*/*/*/*/*/bin/lmp* 2>/dev/null | grep -v '\.so'); do
  p=$(LC_ALL=C grep -aoE 'pppm/(gpu|kk)' "$f" 2>/dev/null | sort -u | tr '\n' ' ')
  [ -z "$p" ] && so=$(ldd "$f" 2>/dev/null | awk '/liblammps/ && $3 ~ /^\//{print $3; exit}') && [ -n "$so" ] && p=$(LC_ALL=C grep -aoE 'pppm/(gpu|kk)' "$so" 2>/dev/null | sort -u | tr '\n' ' ')
  [ -n "$p" ] && echo "  $f: $p"
done
echo "-- modules (cuda, lammps)"; { type module >/dev/null 2>&1 || { [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh; }; module avail cuda lammps 2>&1; } | grep -iE "cuda|lammps" | head -10
# a resource line: the queue of the GPU nodes (else the first queue with "gpu" in its name), one node, all its cores, its GPUs
G=$(echo "$NODES" | awk 'NF {print $2, $3, $4; exit}')
NG=$(echo "$G" | awk '{print $1}'); NC=$(echo "$G" | awk '{print $2}'); Q=$(echo "$G" | awk '{print $3}')
[ -z "$Q" ] && Q=$(qstat -Q 2>/dev/null | awk 'NR > 2 && tolower($1) ~ /gpu/ {print $1; exit}')
if [ -n "$NG" ] && [ -n "$Q" ]; then
  MP=$((NG * ${GPU_RPG:-8})); [ "$MP" -gt "${NC:-$MP}" ] && MP=$NC
  echo "SUGGEST: -q $Q -l select=1:ncpus=${NC}:mpiprocs=${MP}:ngpus=${NG}"
else
  echo "SUGGEST: (none: give GPU_QSUB yourself, e.g. GPU_QSUB=\"-q <GPU queue> -l select=1:ncpus=<cores>:mpiprocs=<ranks>:ngpus=<GPUs>\")"
fi
