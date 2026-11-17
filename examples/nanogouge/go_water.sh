#!/bin/bash
# One command: shear at several water contents (porosities), Lin & Hatano (2024)-type figures against the porosity.
#   ./go_water.sh [N = 22] ["m_w/m_c list, wet to dry" = "2 1 0.5 0.2"] ["normal stresses MPa" = "5.6 56"] [nodes = 1] [run folder] [seed]
#   env WALLTIME=0:30:00 ./go_water.sh 22 "2 1 0.5 0.2" "5.6 56" 1        -> RUN_N22_water/
# One dehydration (lammps/in.water) from the wet start (m_w/m_c = 3, U&B 2020) goes down the list; at each water content the gouge of that
# moment is consolidated (03), sheared (04: RATE to the strain GMAX, then NVS velocity steps x3.16 of GSTEP strain each) and held (05: HOLD
# steps) at each normal stress, in RUN_N<N>_water/w<m_w/m_c>/; then the dehydration goes on.  The PBS chain runs it all (job.sh, as go.sh).
# Porosity (water-filled pores, rho_s 2.7, rho_w 1.0):  phi = (m_w/m_c) / (m_w/m_c + rho_w/rho_s)
#   m_w/m_c   3     2     1     0.5   0.3   0.2
#   phi       0.89  0.84  0.73  0.57  0.45  0.35
# Options (env VAR=... ./go_water.sh ..., as go.sh; kept in the folder's config.env):
#   BY=phi        the list is porosities, e.g.  env BY=phi ./go_water.sh 22 "0.85 0.7 0.55 0.4"
#   RATE=1e9      first shear rate in 1/s (default 1e9);  GMAX=1  its strain (default 1);  NVS=2  velocity steps (default 2: 1e9 -> 1e10 /s);
#   GSTEP=0.5     strain of each velocity step;  HOLD=20000  steps of the hold after the shear (default 20000 = 20 ps; HEAL=0: no hold)
#   NRUN=1000     a quick test of everything: every run NRUN steps (also the shear and the hold), folder RUN_N<N>_water_test<NRUN>
#   DFRAC=0.1     each water removal takes this fraction of the initial water (default 0.1)
#   DT, DUMP, FLOW, VFAC: as go.sh
# Figures (the jobs draw them after each job; by hand:  ./sifpy analysis/lin_figs.py RUN_N22_water):  RUN_N22_water/figs/lin_*.png
set -eu
cd "$(dirname "$0")"
N=${1:-22}; WL=${2:-"2 1 0.5 0.2"}; SIG=${3:-"5.6 56"}; NOD=${4:-1}; RUNARG=${5:-}; SEED=${6:-1}
if [ "${BY:-mwmc}" = phi ]; then
  WL=$(for p in $WL; do awk -v p="$p" 'BEGIN{ if (p <= 0 || p >= 1) { print "bad" } else printf "%.4g\n", p * 1.0 / ((1 - p) * 2.7) }'; done | tr '\n' ' ')
  echo "porosities ${2:-} -> m_w/m_c $WL"
fi
# wet to dry, each below the start (3) and above 0, at most 8
for w in $WL; do
  awk -v w="$w" 'BEGIN{ exit !(w ~ /^[0-9.eE+-]+$/ && w + 0 > 0 && w + 0 < 3) }' || { echo "$w: every water content must be 0 < m_w/m_c < 3 (porosity 0 < phi < 0.89)"; exit 1; }
done
WL=$(for w in $WL; do awk -v w="$w" 'BEGIN{print w + 0}'; done | sort -g -r -u | tr '\n' ' ')
NWC=$(echo $WL | wc -w); [ "$NWC" -ge 1 ] && [ "$NWC" -le 8 ] || { echo "1 to 8 water contents, not $NWC"; exit 1; }
echo "water contents m_w/m_c (wet to dry): $WL"
echo "porosities:                          $(for w in $WL; do awk -v w=$w 'BEGIN{printf "%.2f ", w / (w + 1.0 / 2.7)}'; done)"
if [ -n "${NRUN:-}" ]; then
  RUN=${RUNARG:-RUN_N${N}_water_test$NRUN}
  export GMAX=${GMAX:-0} HOLD=${HOLD:-0}
else
  RUN=${RUNARG:-RUN_N${N}_water}
  export RATE=${RATE:-1e9} GMAX=${GMAX:-1} NVS=${NVS:-2} GSTEP=${GSTEP:-0.5} HOLD=${HOLD:-20000}
fi
export WLIST="$WL"
exec ./go.sh "$N" "$SIG" water "$NOD" "$RUN" "$SEED"
