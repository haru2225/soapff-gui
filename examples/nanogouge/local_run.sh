#!/bin/bash
# A run folder of go.sh (stages "one", lammps/in.all) continued on a PC without PBS: one run per LAMMPS start with a time limit (a save
# every TSAVE s: all.restart + all.state), started again until all.done.  Run it in tmux, so a lost connection does not stop it:
#   tmux new -s ng1 'bash local_run.sh RUN_N54_nro_w0.05 8; bash'          (detach: Ctrl-b d;  back: tmux attach -t ng1)
#   bash local_run.sh <folder> [MPI ranks = all hardware threads] [OpenMP threads per rank = 1]
# Options (env): LMPSF="-sf gpu -pk gpu 1" (a LAMMPS with the GPU package: the pair and PPPM on the GPU), LMP=/usr/bin/lmp (the LAMMPS: MOLECULE, KSPACE, RIGID), MPIRUN=mpirun, TSAVE=3600, MAXFAIL=3 (failed starts in a row, then
# it stops: local.failed).  The options of the folder are those of its config.env (EXTRA lines of go.sh) and chain.env (SIGMAS); a line
# changed there counts from the next LAMMPS start.  Log: <folder>/local.log (also the errors of all ranks); rank 0: <folder>/lmp.screen.local
set -u
D=${1:?usage: bash local_run.sh <run folder> [MPI ranks] [OpenMP threads per rank]}
NP=${2:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}; OMP=${3:-1}
LMP=${LMP:-$(command -v lmp || command -v lmp_mpi)}; MPIRUN=${MPIRUN:-mpirun}
TSAVE=${TSAVE:-3600}; MAXFAIL=${MAXFAIL:-3}
cd "$D" || exit 1
[ -f in.all ] || { echo "$D: no in.all"; exit 1; }
rm -f local.failed
MPIOPT="--bind-to none"; $MPIRUN --version 2>&1 | grep -qi "open mpi" && MPIOPT="$MPIOPT --use-hwthread-cpus --oversubscribe"
fails=0
while [ ! -f all.done ]; do
  EXTRA=""; eval "$(grep '^EXTRA=' config.env 2>/dev/null)"
  [ -f lammps_seed ] && EXTRA="-var seed $(cat lammps_seed) $EXTRA"
  SIGMAS=5.6; eval "$(grep '^SIGMAS=' chain.env 2>/dev/null)"; set -- $SIGMAS
  RES=0; [ -f all.state ] && RES=1
  SF=${LMPSF:-}; [ "$OMP" -gt 1 ] && SF="$SF -sf omp -pk omp $OMP"
  echo "== $(date '+%F %T')  next: $(head -1 all.state 2>/dev/null | awk '{print $4}')  $NP ranks x $OMP threads, $LMP" | tee -a local.log
  OMP_NUM_THREADS=$OMP $MPIRUN $MPIOPT -np $NP $LMP $SF -in in.all -var data gouge.data -var settings gouge.settings \
      -var sig1 $1 -var sig2 ${2:-0} -var resume $RES -var tlimit $TSAVE $EXTRA -screen lmp.screen.local >> local.log 2>&1
  rc=$?
  L=$(grep -E "^ +[0-9]+ +[0-9.-]+ " log.one | tail -1 | awk '{print "step " $1 ", " $NF " steps/s"}')
  if [ $rc -ne 0 ]; then
    fails=$((fails + 1))
    echo "LAMMPS failed (status $rc; $fails of $MAXFAIL in a row; last: $L): $(grep -h "ERROR" lmp.screen.local local.log 2>/dev/null | tail -2 | tr '\n' ' ')" | tee -a local.log
    [ $fails -ge $MAXFAIL ] && { echo "stopped: $MAXFAIL failures in a row (local.failed).  Again: bash local_run.sh $D $NP $OMP" | tee -a local.log; date > local.failed; exit 1; }
    sleep 30
  else
    fails=0; echo "   ok ($L)" | tee -a local.log
  fi
done
echo "== $(date '+%F %T')  all done (all.done)" | tee -a local.log
