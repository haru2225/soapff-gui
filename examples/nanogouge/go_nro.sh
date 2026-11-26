#!/bin/bash
# NRO: the initial random state of a fabric run, sheared as Lin & Hatano (2024, JGR Solid Earth 129, e2024JB028850), at several water
# contents at once (one folder and one PBS chain each, so the jobs run side by side), and their figures for all of them together.
#   ./go_nro.sh <fabric folder> ["m_w/m_c list" = "0 0.05"] ["normal stresses MPa" = "5.6"] [nodes per job = 16]
#   env WALLTIME=24 ./go_nro.sh RUN_N54_fab_dis "0 0.05" "5.6" 16      -> RUN_N54_nro_w0 (dry) and RUN_N54_nro_w0.05 (5 % = m_w/m_c 0.05)
# Each folder (go.sh FROMFAB=<fabric> NRO=1 MWMC0=<w>): the particles of the fabric at random as they were before the dehydration
# (<fabric>/01_equil.data), all their water out, m_w/m_c <w> of water on the clay surfaces first (rewater.py, a grid); then in.all:
#   01  minimisation, NVT NEQ steps at the volume of the random state; compaction the same in x y z at CRATE until the pressure is the
#       first normal stress; NPT NPT01 steps at it, x y z each to the pressure (C01=aniso: the flexible box of Lin & Hatano)
#   03  NPzT at sigma in pieces of 20 ps until Lz changes by < 0.05 % over a piece, at most NCONSMAX pieces
#   04  shear: Lees-Edwards (the box tilts), Langevin on the deviation from the linear profile, barostat on z (as Lin & Hatano); RATE to
#       the strain GMAX, then NVS velocity steps x VFAC of GSTEP strain each
#   05  hold HOLD steps at sigma
# Defaults (about 4 ns per folder: 01 + 03 ~0.6 ns, 04 3.3 ns), any of them given as env VAR=... instead:
#   NEQ=50000 CRATE=2e9 NPT01=150000 C01=aniso NCONSMAX=10  RATE=3.16e9 GMAX=10 NVS=2 VFAC=3.16227766 GSTEP=1 (3.16e9 /s to the strain 10,
#   then 1e10 and 3.16e10 /s, one unit of strain each) HOLD=50000;  WALLTIME=24 MAXJOBS=2 (jobs per folder)
#   NAME=nro: the folders RUN_N<N>_<NAME>_w<w>
# Figures: each job draws those of its folder (figs/: hot_spots.py) and Lin & Hatano's of all the folders (figs/lin_*.png, lin_figs.py;
# lin_group.txt lists them).  By hand at any time:  ./figs.sh RUN_N54_nro_w0   (or ./sifpy analysis/lin_figs.py RUN_N54_nro_w0 RUN_N54_nro_w0.05)
# The same command again: the folders that exist go on from where they are (with the defaults above or the env VAR=... given).
set -eu
cd "$(dirname "$0")"
FAB=${1:-}; WL=${2:-"0 0.05"}; SIG=${3:-"5.6"}; NOD=${4:-16}
[ -n "$FAB" ] || { echo "usage: ./go_nro.sh <fabric folder, e.g. RUN_N54_fab_dis> [\"0 0.05\"] [\"5.6\"] [nodes = 16]"; exit 1; }
FAB=${FAB%/}
[ -d "$FAB" ] || { echo "$FAB: no such folder here"; exit 1; }
[ -f "$FAB/01_equil.data" ] || [ -f "$FAB/gouge.data" ] || { echo "$FAB: no 01_equil.data or gouge.data (the initial random state)"; exit 1; }
N=${N:-$(echo "$FAB" | sed -n 's/^RUN_N\([0-9][0-9]*\).*/\1/p')}
[ -n "$N" ] || { echo "the number of particles: N=<n> ./go_nro.sh ... (not in the folder name $FAB)"; exit 1; }
NAME=${NAME:-nro}
for w in $WL; do
  awk -v w="$w" 'BEGIN{ exit !(w ~ /^[0-9.]+$/ && w + 0 >= 0 && w + 0 <= 2) }' || { echo "$w: a water content 0 <= m_w/m_c <= 2"; exit 1; }
done
WL=$(for w in $WL; do awk -v w="$w" 'BEGIN{printf "%g\n", w}'; done | sort -g -u | tr '\n' ' ')
RUNS=$(for w in $WL; do echo "RUN_N${N}_${NAME}_w$w"; done)
# the protocol (as Lin & Hatano 2024, shortened for an all-atom gouge)
LIN="NEQ=50000 CRATE=2e9 NPT01=150000 C01=aniso NCONSMAX=10 RATE=3.16e9 GMAX=10 NVS=2 VFAC=3.16227766 GSTEP=1 HOLD=50000"
export WALLTIME=${WALLTIME:-24} MAXJOBS=${MAXJOBS:-2}
for w in $WL; do
  RUN=RUN_N${N}_${NAME}_w$w
  # each time (go.sh sets the options of a loose structure again when it is given again): the defaults, or the env VAR given
  ENVS=""
  for kv in $LIN; do k=${kv%%=*}; ENVS="$ENVS $k=${!k:-${kv#*=}}"; done
  echo "==== $RUN: m_w/m_c $w ($(awk -v w=$w 'BEGIN{ if (w == 0) print "dry"; else printf "%g %% water", 100 * w }'))  $ENVS"
  env $ENVS FROMFAB="$FAB" NRO=1 MWMC0="$w" ./go.sh "$N" "$SIG" one "$NOD" "$RUN"
  echo "$RUNS" > "$RUN/lin_group.txt"
done
echo
echo "two chains (one per folder), $NOD nodes each: qstat -u \$USER"
for r in $RUNS; do echo "  progress: tail -2 $r/log.one;  grep 'in.all: 0' $r/log.one | tail -2"; done
echo "figures (each job draws them; by hand at any time): ./figs.sh $(echo $RUNS | awk '{print $1}')   -> <folder>/figs/lin_*.png"
