#!/bin/bash
#PBS -j oe
# One job of the chain (submit_chain.sh): it works on the first unfinished stage (01, 02, then 03 04 05 for each sigma of chain.env),
# continues it from its last restart, stops LAMMPS cleanly before the wall time, and the next job goes on.  Any wall time works (1 h too).
# STAGE=01..05 SIGMA=<MPa> (run.sh, qsub -v) runs one given stage only.  EXTRA = more LAMMPS -var options.
# It runs in the submit directory (the run folder holding the in.* files, gouge.data and gouge.settings).
set -u
# go to the run folder.  The login node's path may not exist on the compute node (symlink / automount): try the physical path too, retry, and stop with a diagnosis rather than running in the wrong place.
for d in "${PBS_O_WORKDIR:-}" "${SUBDIR_PHYS:-}" "${SUBDIR_LOGICAL:-}"; do
  [ -n "$d" ] || continue
  for t in 1 2 3; do cd "$d" 2>/dev/null && [ -f ./in.01_equil ] && break 2; sleep 5; done
done
if [ ! -f ./in.01_equil ]; then
  echo "ERROR: cannot enter the run folder on $(hostname).  tried: PBS_O_WORKDIR=${PBS_O_WORKDIR:-} SUBDIR_PHYS=${SUBDIR_PHYS:-}"
  echo "pwd=$PWD  HOME=$HOME"; for d in "${PBS_O_WORKDIR:-}" "${SUBDIR_PHYS:-}"; do [ -n "$d" ] && ls -ld "$d" "$(dirname "$d")" 2>&1; done
  ls -ld /home /home/* 2>&1 | head -5; df -h "$HOME" 2>&1 | tail -2; exit 1
fi
# everything this job prints also goes at once to jobs.log in the run folder (PBS writes its output file ng_*.o* only when the job has ended)
exec > >(tee -a jobs.log) 2>&1
echo "===== job ${PBS_JOBID:-$$} on $(hostname)   start $(date)"
QSUB_OPTS=""; MODULES=""; SIF=""; SING_OPTS=""; LMP="lmp"; LMP_IN_IMAGE="lmp"; MPIRUN="mpirun"; MPI_OPTS=""; HOSTFILE_FMT="none"; CORES_PER_NODE="${CORES_PER_NODE:-}"
LMP_LIBDIRS=""; MPIRUN_OK=""
[ -f ./config.auto ] && source ./config.auto
# lmp.auto: the LAMMPS build and the MPI launcher an earlier job of this folder found and checked (no search again; delete lmp.auto to search again)
[ -f ./lmp.auto ] && source ./lmp.auto
[ -f ./config.env ] && source ./config.env
# a command with a time limit (a search through a large shared file system must not eat the job's wall time)
tmo() { local t=$1; shift; if command -v timeout >/dev/null 2>&1; then timeout "$t" "$@"; else "$@"; fi; }
SEARCH_ROOTS=${SEARCH_ROOTS:-"/home/app /opt /usr/local $HOME /work"}   # where other LAMMPS builds and libraries are looked for
# safety: MPI_OPTS must hold launcher options only; cut at the first word that is the program name or a redirect (old config.auto files)
if [ -n "${MPI_OPTS}" ]; then
  MPI_OPTS=$(echo "${MPI_OPTS}" | awk '{o=""; for(i=1;i<=NF;i++){ if ($i ~ /(^|\/)(lmp[_a-z0-9]*|lammps[_a-z0-9]*)$/ || $i ~ /^[<>]/ || $i ~ /^\$/ ) break; o=o" "$i } sub(/^ /,"",o); print o}')
fi
if [ -n "${MODULES}" ]; then set +u; eval "${MODULES}"; set -u; fi
[ -n "${LMP_LIBDIRS}" ] && export LD_LIBRARY_PATH="${LMP_LIBDIRS}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1
HOSTOPT=""
NPROC=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN)
if [ -n "${PBS_NODEFILE:-}" ]; then
  # cores per node = the largest of: slots listed per host in $PBS_NODEFILE, the cores this job may use (nproc), the cores requested (NCPUS)
  SLOTS=$(sort "$PBS_NODEFILE" | uniq -c | awk '{print $1}' | sort -rn | head -1)
  NREQ=${NCPUS:-${PBS_NCPUS:-0}}
  C=${CORES_PER_NODE:-$(echo "$SLOTS $NPROC $NREQ" | tr ' ' '\n' | sort -rn | head -1)}
  NODES=$(sort -u "$PBS_NODEFILE" | wc -l | tr -d " ")
  echo "cores: slots/host in PBS_NODEFILE=${SLOTS}  nproc=${NPROC}  NCPUS=${NREQ}  -> ${C} per node x ${NODES} nodes"
  [ "${C}" -le 1 ] && echo "WARNING: only 1 core per node. Ask for more cores in the resource line (e.g. -l select=${NODES}:ncpus=<cores>:mpiprocs=<cores>); see QSUB_OPTS in config.auto" 
  # default: no host option (a PBS-aware mpirun reads the allocation itself). Only if config sets HOSTFILE_FMT=openmpi|mpich a host file is passed.
  HF=hostfile.${PBS_JOBID:-$$}
  if [ "${HOSTFILE_FMT}" = mpich ]; then sort -u "$PBS_NODEFILE" | awk -v c="$C" '{print $1":"c}' > "$HF"; HOSTOPT="-machinefile $HF"
  elif [ "${HOSTFILE_FMT}" = openmpi ]; then sort -u "$PBS_NODEFILE" | awk -v c="$C" '{print $1" slots="c}' > "$HF"; HOSTOPT="--hostfile $HF"; fi
  NP=${NP:-$((NODES*C))}
else
  NP=${NP:-$NPROC}
fi
if [ -n "${SIF}" ]; then LMPCMD="singularity exec ${SING_OPTS} --bind $PWD ${SIF} ${LMP_IN_IMAGE}"; else LMPCMD="${LMP}"; fi
if [ -z "${SIF}" ]; then
  # the program is given to mpirun with its absolute path (the other nodes may not have the module's PATH)
  LB=${LMP%% *}; LMPABS=$(command -v "$LB" 2>/dev/null)
  [ -z "$LMPABS" ] && [ -x "$LB" ] && LMPABS=$(readlink -f "$LB" 2>/dev/null || echo "$LB")
  if [ -z "$LMPABS" ]; then
    echo "note: '$LB' (from the centre's sample script) is not in PATH on this node; searching for it ..."
    for root in $SEARCH_ROOTS; do
      [ -d "$root" ] || continue
      LMPABS=$(tmo 120 find "$root" -maxdepth 6 -type f -name "$(basename "$LB")" 2>/dev/null | head -1); [ -n "$LMPABS" ] && break
    done
  fi
  if [ -z "$LMPABS" ]; then
    echo "ERROR: LAMMPS executable '$LB' not found. Run on the login node:  which $LB;  module avail 2>&1 | grep -i lammps;  cat /home/app/lammps/run.sh"
    echo "then put the module/PATH lines in config.env (MODULES=...) or the full path (LMP=/full/path/$LB)"; exit 1
  fi
  echo "LAMMPS: $LMPABS"; LMPCMD="$LMPABS"
  ADDED_LIBDIRS=""
  # shared libraries a binary needs but cannot find (e.g. libfftw3.so.3): look for them in the software area and add their directories
  fixlibs() {
    local BIN=$1 MISS LMPDIR l f root
    MISS=$(ldd "$BIN" 2>/dev/null | awk '/not found/{print $1}'); [ -z "$MISS" ] && return 0
    echo "missing shared libraries of $BIN: $(echo $MISS)"
    LMPDIR=$(dirname "$(dirname "$BIN")")
    for l in $MISS; do
      f=""
      for root in "$LMPDIR" "$(dirname "$LMPDIR")" /home/app /opt /usr/lib64 /usr/local "$HOME" /work; do
        [ -d "$root" ] && f=$(tmo 60 find "$root" -maxdepth 8 -name "${l}*" 2>/dev/null | head -1) && [ -n "$f" ] && break
      done
      if [ -n "$f" ]; then
        export LD_LIBRARY_PATH="$(dirname "$f")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"; ADDED_LIBDIRS="$(dirname "$f")${ADDED_LIBDIRS:+:$ADDED_LIBDIRS}"
        echo "  $l -> $(dirname "$f")"
      else echo "  $l: not found anywhere; try 'module avail 2>&1 | grep -i ${l#lib}' and add the module to config.env (MODULES=...)"; fi
    done
    ldd "$BIN" 2>/dev/null | grep "not found" && echo "WARNING: still missing libraries (see above)"
    return 0
  }
  fixlibs "$LMPABS"
  # the shear stage needs "fix press/langevin" (LAMMPS 21Nov2023 or newer).  The style name is stored in the program (or in its liblammps.so), so it is
  # looked for there; running "lmp -h" printed nothing when a library was missing or MPI did not start without mpirun, and a good build looked old.
  has_pl() {
    LC_ALL=C grep -aq 'press/langevin' "$1" 2>/dev/null && return 0
    local so; so=$(ldd "$1" 2>/dev/null | awk '/liblammps/ && $3 ~ /^\//{print $3; exit}')
    [ -n "$so" ] && LC_ALL=C grep -aq 'press/langevin' "$so" 2>/dev/null
  }
  # can it start here?  A library not found, or a build made with another MPI version (lmp_icex 3Mar20 with HPE MPT 2.33:
  # "symbol lookup error: undefined symbol: mpi_sgi_status_ignore"), stops before LAMMPS begins.  ldd -r resolves every symbol without running it.
  load_errors() { command -v ldd >/dev/null 2>&1 || return 0; ldd -r "$1" 2>&1 | grep -E 'not found|undefined symbol' | sort -u | head -5; }
  E=$(load_errors "$LMPABS")
  if ! has_pl "$LMPABS" || [ -n "$E" ]; then
    # normal on the first job with the centre's lmp_icex (3Mar20): another build is chosen here and kept in lmp.auto for the next jobs
    has_pl "$LMPABS" || echo "note: $LMPABS has no 'fix press/langevin' (LAMMPS older than 21Nov2023)"
    [ -n "$E" ] && { echo "note: $LMPABS cannot start in this environment:"; echo "$E" | sed 's/^/    /'; }
    echo "looking for another LAMMPS build (with press/langevin, starting here) ..."
    for root in $SEARCH_ROOTS; do
      [ -d "$root" ] || continue
      for c in $(tmo 120 find "$root" -maxdepth 7 -type f \( -name 'lmp' -o -name 'lmp_*' \) -perm -u+x 2>/dev/null | grep -v '\.so' | sort -r | head -30); do
        [ "$c" = "$LMPABS" ] && continue
        has_pl "$c" || continue
        OLDLD=${LD_LIBRARY_PATH:-}; OLDADD=$ADDED_LIBDIRS; fixlibs "$c" >/dev/null 2>&1
        CE=$(load_errors "$c")
        if [ -z "$CE" ]; then LMPABS=$c; LMPCMD=$c; echo "using $c (has press/langevin; the notes above are no problem)"; break 2; fi
        echo "  skipped $c:"; echo "$CE" | sed 's/^/    /'; export LD_LIBRARY_PATH="$OLDLD"; ADDED_LIBDIRS=$OLDADD
      done
    done
    has_pl "$LMPABS" || echo "WARNING: no LAMMPS with press/langevin found; stages 01-03 run, stage 04 (shear) will fail. Put LMP=/path/to/newer/lmp in config.env"
  fi
  E=$(load_errors "$LMPABS")
  if [ -n "$E" ]; then
    echo "WARNING: $LMPABS will probably not start:"; echo "$E" | sed 's/^/    /'
    echo "$E" | grep -q 'undefined symbol: mpi' && echo "  it was built with another MPI version than the one loaded ($(command -v "${MPIRUN%% *}" 2>/dev/null)): use another LAMMPS (LMP=... in config.env) or load the MPI it was built with (MODULES=... in config.env)"
    echo "$E" | grep -q 'not found' && echo "  a library is missing: find it (e.g. find /home/app -name 'libfftw3.so*') and add its folder in config.env:  MODULES='export LD_LIBRARY_PATH=/that/folder:\$LD_LIBRARY_PATH'"
  fi
  echo "MPI library of LAMMPS: $(ldd "$LMPABS" 2>/dev/null | awk '/libmpi/{print $3; exit}')   launcher: $(command -v "${MPIRUN%% *}" 2>/dev/null)"

  # ---- are the MPI ranks spread over all nodes of the job?  HPE MPT's plain "mpirun -np N prog" may start every rank on the first node
  #      (4 nodes asked, 1 used, 4x slower); its PBS launcher mpiexec_mpt uses them all.  A tiny LAMMPS run (a few s) writes the node of each rank.
  nodes_used() {   # $@ = launcher with its options: the number of different nodes the ranks ran on (0 = it did not start)
    rm -f probe_procs.txt
    printf 'processors * * * file probe_procs.txt\nregion b block 0 10 0 10 0 10\ncreate_box 1 b\n' > probe.in
    tmo 180 "$@" ${LMPCMD} -in probe.in -log none -screen none >/dev/null 2>&1 </dev/null
    if [ -f probe_procs.txt ]; then awk -F': ' '$1 ~ /^[0-9]/{print $3}' probe_procs.txt | sort -u | wc -l | tr -d ' '; else echo 0; fi
    rm -f probe.in probe_procs.txt
  }
  NODES=${NODES:-1}; MPIRUN_FROM_NOW=$MPIRUN
  if [ "$NODES" -gt 1 ] && [ "${MPI_CHECK:-1}" = 1 ]; then
    if [ -n "${MPIRUN_OK}" ] && [ "${MPIRUN}" = "${MPIRUN_FROM:-}" ]; then
      [ "$MPIRUN_OK" != "$MPIRUN" ] && { MPIRUN=$MPIRUN_OK; MPI_OPTS=""; HOSTOPT=""; }
      echo "MPI launcher: ${MPIRUN} (spreads the ranks over all nodes; checked by an earlier job, lmp.auto)"
    else
      MPIRUN_OK=""
      U=$(nodes_used ${MPIRUN} -np ${NP} ${MPI_OPTS} ${HOSTOPT})
      echo "MPI launcher check: '${MPIRUN} -np ${NP}' puts the ranks on ${U} of ${NODES} nodes"
      if [ "$U" -ge "$NODES" ]; then MPIRUN_OK=$MPIRUN
      elif command -v mpiexec_mpt >/dev/null 2>&1; then
        U2=$(nodes_used mpiexec_mpt -np ${NP})
        echo "MPI launcher check: 'mpiexec_mpt -np ${NP}' (HPE MPT's launcher for PBS) puts them on ${U2} of ${NODES} nodes"
        if [ "$U2" -gt "$U" ]; then MPIRUN=mpiexec_mpt; MPI_OPTS=""; HOSTOPT=""; echo "-> using mpiexec_mpt"; fi
        [ "$U2" -ge "$NODES" ] && MPIRUN_OK=mpiexec_mpt
      fi
      [ -z "$MPIRUN_OK" ] && echo "WARNING: the ranks do not reach every node of the job (see the lines above).  Set MPIRUN (and HOSTFILE_FMT) in config.env."
    fi
  fi
  # keep what was found for the next jobs (no search, no check)
  if [ -z "$(load_errors "$LMPABS")" ] && has_pl "$LMPABS"; then
    {
      echo "# written by job.sh $(date): the LAMMPS build and MPI launcher found and checked in this folder.  Delete this file to search again."
      echo "LMP=\"$LMPABS\""
      echo "LMP_LIBDIRS=\"${ADDED_LIBDIRS}${LMP_LIBDIRS:+${ADDED_LIBDIRS:+:}$LMP_LIBDIRS}\""
      [ -n "${MPIRUN_OK}" ] && echo "MPIRUN_FROM=\"${MPIRUN_FROM_NOW}\"; MPIRUN_OK=\"${MPIRUN_OK}\""
    } > lmp.auto
  fi
fi
RUN="${MPIRUN} -np ${NP} ${MPI_OPTS} ${HOSTOPT} ${LMPCMD} -screen none"
SIGMA=${SIGMA:-5.6}; EXTRA=${EXTRA:-}
[ -f ./chain.env ] && source ./chain.env          # SIGMAS and STAGES of the chain (written by submit_chain.sh)
SIGMAS=${SIGMAS:-$SIGMA}; STAGES=${STAGES:-all}
echo "MPI ${NP} (${C:-?} per node)  $(date)"
echo "${RUN}"

# ---- wall time.  LAMMPS gets the time left minus TIME_MARGIN (s), stops cleanly, writes a restart; the next job continues from it.
WALL=""
[ -n "${PBS_JOBID:-}" ] && WALL=$(qstat -f "$PBS_JOBID" 2>/dev/null | awk -F'= ' '/Resource_List.walltime/{print $2; exit}')
[ -z "$WALL" ] && WALL=$(echo "${QSUB_OPTS:-}" | grep -o 'walltime=[0-9:]*' | head -1 | cut -d= -f2)
tosec() { echo "$1" | awk -F: '{s=0; for(i=1;i<=NF;i++) s=s*60+$i; print int(s)}'; }
WALLSEC=${WALLSEC:-$([ -n "$WALL" ] && tosec "$WALL" || echo 0)}
[ "$WALLSEC" -le 0 ] && [ -n "${PBS_WALLTIME:-}" ] && WALLSEC=$(tosec "$PBS_WALLTIME")
if [ "$WALLSEC" -le 0 ] && [ -n "${PBS_JOBID:-}" ]; then
  # not known: assume a short one rather than letting LAMMPS run into the kill (a stop too early only costs one more job)
  WALLSEC=${DEFAULT_WALLSEC:-3600}; WALL="unknown, assumed"; echo "WARNING: the wall time of this job was not found (qstat -f); assuming ${WALLSEC} s (DEFAULT_WALLSEC in config.env)"
fi
MARGIN=${TIME_MARGIN:-360}
echo "wall time of this job: ${WALL:-unknown} (${WALLSEC} s); LAMMPS stops ${MARGIN} s before its end"
tleft() { if [ "$WALLSEC" -gt 0 ]; then echo $(( WALLSEC - SECONDS - MARGIN )); else echo off; fi; }

donefile() { case $1 in 01) echo 01_equil.done;; 02) echo 02_dehydrate.done;; 03) echo 03_consolidate_$2.done;; 04) echo shear_$2.done;; 05) echo heal_$2.done;; esac; }
want() { [ "$STAGES" = all ] || [[ " $STAGES " == *" $1 "* ]]; }
next_work() {   # "STAGE SIGMA" of the first unfinished stage: 01, 02, then 03 04 05 for each sigma in turn; nothing when all is done
  local st S
  for st in 01 02; do want $st && [ ! -f "$(donefile $st)" ] && { echo "$st -"; return; }; done
  for S in $SIGMAS; do for st in 03 04 05; do want $st && [ ! -f "$(donefile $st $S)" ] && { echo "$st $S"; return; }; done; done
}
newest() { local b="" f; for f in "$@"; do [ -f "$f" ] && { [ -z "$b" ] || [ "$f" -nt "$b" ]; } && b=$f; done; echo "$b"; }   # the first one wins a tie (.stop)

# the input of a stage is missing: say which earlier stage did not finish and why, then stop
missing_input() {
  local st d
  echo "ERROR: stage ${STAGE} cannot start: $1 is missing."
  echo "  finished stages: $(ls *.done 2>/dev/null | tr '\n' ' ')"
  for st in 01 02 03 04; do
    d=$(donefile $st "$SIGMA"); [ -f "$d" ] && continue
    echo "  the first unfinished stage is ${st}. The end of its LAMMPS log and of the last job outputs:"
    for f in $(ls log.${st}* 2>/dev/null) $(ls -t ng*.o* 2>/dev/null | head -2); do
      echo "  --- $f"; grep -E "ERROR|walltime|Killed|Lost atoms|Out of range|Bond atoms missing|cannot|not found" "$f" | tail -4 | cut -c1-200 | grep . || tail -3 "$f" | cut -c1-200
    done
    break
  done
  exit 1
}

# run (or continue) one stage until it is done or the time is up.  A restart of the stage is continued; otherwise it starts from the data file of the stage before.
run_stage() {
  local st=$1 S=$2 R TL rc
  TL=$(tleft)
  if [ "$TL" != off ] && [ "$TL" -lt "${MIN_LEFT:-180}" ]; then echo "only ${TL} s left: stage $st is left to the next job"; return 1; fi
  STAGE=$st; SIGMA=$S
  echo "== stage $st  sigma $S  $(date)  LAMMPS time limit ${TL} s"
  local V="-var tlimit $TL ${EXTRA}"
  case $st in
  01) R=$(newest 01_equil.restart.stop 01_equil.restart.1 01_equil.restart.2 01_equil.restart.min)
      if [ -z "$R" ]; then ${RUN} -in in.01_equil -var data gouge.data -var settings gouge.settings $V
      elif [ "$R" = 01_equil.restart.min ]; then ${RUN} -in in.01_equil -var resume 1 -var rfile $R -var domin 1 $V
      else ${RUN} -in in.01_equil -var resume 1 -var rfile $R -var domin 0 $V; fi ;;
  02) R=$(newest dehyd.restart.stop dehyd.restart.1 dehyd.restart.2)
      if [ -n "$R" ]; then ${RUN} -in in.02_dehydrate -var resume 1 -var rfile $R -var nw0 "$(cat dehyd_nw0.txt)" $V
      else [ -f 01_equil.data ] || missing_input 01_equil.data
           ${RUN} -in in.02_dehydrate -var data 01_equil.data -var settings gouge.settings $V; fi ;;
  03) R=$(newest 03_consolidate_$S.restart.stop 03_consolidate_$S.restart.1 03_consolidate_$S.restart.2)
      if [ -n "$R" ]; then ${RUN} -in in.03_consolidate -var resume 1 -var rfile $R -var sigma $S $V
      else [ -f 02_dehydrate.data ] || missing_input 02_dehydrate.data
           ${RUN} -in in.03_consolidate -var data 02_dehydrate.data -var settings gouge.settings -var sigma $S $V; fi ;;
  04) R=$(newest shear_$S.restart.stop shear_$S.restart.1 shear_$S.restart.2)
      if [ -n "$R" ]; then ${RUN} -in in.04_shear -var resume 1 -var rfile $R -var sigma $S $V
      else [ -f 03_consolidate_$S.data ] || missing_input 03_consolidate_$S.data
           ${RUN} -in in.04_shear -var data 03_consolidate_$S.data -var settings gouge.settings -var sigma $S $V; fi ;;
  05) R=$(newest heal_$S.restart.stop heal_$S.restart.1 heal_$S.restart.2)
      if [ -n "$R" ]; then ${RUN} -in in.05_heal -var resume 1 -var rfile $R -var sigma $S $V
      else [ -f shear_${S}_final.data ] || missing_input shear_${S}_final.data
           ${RUN} -in in.05_heal -var data shear_${S}_final.data -var settings gouge.settings -var sigma $S $V; fi ;;
  *)  echo "unknown STAGE $st"; exit 1 ;;
  esac
  rc=$?
  [ -f procs.txt ] && echo "MPI ranks per node: $(awk -F': ' '$1 ~ /^[0-9]/{print $3}' procs.txt | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"
  # how far this job got: the last thermo line of the stage (step ... last column = time steps per second)
  local L; L=$(ls -t log.${st}* 2>/dev/null | head -1)
  [ -n "$L" ] && echo "LAMMPS ended $(date) with status $rc.  last thermo line of $L (step temp pe ... steps/s):" && grep -E "^ +[0-9]+ +[0-9.-]+ " "$L" | tail -1 | awk '{print "   step", $1, "  T", $2, "  steps/s", $NF}'
  if [ $rc -ne 0 ]; then
    # LAMMPS stopped with an error (a stop by the time limit ends with status 0): do not let the following jobs repeat it
    echo "stage $st sigma $S job ${PBS_JOBID:-} $(date): LAMMPS ended with status $rc" > chain.failed
    echo "ERROR: LAMMPS ended with status $rc.  The chain stops here (chain.failed); the following jobs end at once."
    for f in $(ls -t log.${st}* 2>/dev/null | head -1); do echo "--- $f"; grep -E "ERROR|WARNING" "$f" | tail -5; tail -3 "$f"; done
    echo "Fix the cause, then run ./go.sh again (it deletes chain.failed and continues from the last restart)."
    exit 1
  fi
  return 0
}

# the last job of the chain submits the next one while work is left (needs qsub on the compute node; otherwise ./go.sh again).
# $1 = the job it waits for (default: this one)
extend_chain() {
  [ -n "${PBS_JOBID:-}" ] && [ -f chain.tail ] && [ "$(cat chain.tail)" = "$PBS_JOBID" ] || return 0
  local J dep=${1:-$PBS_JOBID}
  if J=$(qsub -N ng_next ${QSUB_OPTS} -W depend=afterany:${dep} -v "SUBDIR_PHYS=$(pwd -P),SUBDIR_LOGICAL=$PWD" job.sh 2>&1); then
    echo "$J" > chain.tail; echo "submitted the next job of the chain: $J (starts after ${dep})"
  else
    echo "WARNING: could not submit the next job from the compute node ($J).  The chain ends after the jobs already queued: run ./go.sh again to continue."
  fi
}

# one job in this folder at a time (two chains would write the same files)
LOCK=.chain.lock
if [ -f $LOCK ]; then
  H=$(cat $LOCK)
  if [ -n "$H" ] && [ "$H" != "${PBS_JOBID:-$$}" ] && qstat "$H" 2>/dev/null | grep -q " R "; then
    echo "job $H is working in this folder; this job ends"; extend_chain "$H"; exit 0
  fi
fi
echo "setup of this job: ${SECONDS} s"
echo "${PBS_JOBID:-$$}" > $LOCK; trap '[ "$(cat $LOCK 2>/dev/null)" = "${PBS_JOBID:-$$}" ] && rm -f $LOCK' EXIT
[ -f chain.failed ] && { echo "an earlier job failed: $(cat chain.failed).  Fix it and run ./go.sh again (it deletes chain.failed)."; exit 1; }

if [ -n "${STAGE:-}" ] && [ "${STAGE}" != auto ]; then
  # one given stage (run.sh, qsub -v STAGE=..): continue it until it is done or the time is up
  [ -f "$(donefile $STAGE $SIGMA)" ] && { echo "stage $STAGE (sigma $SIGMA) is already done"; exit 0; }
  run_stage $STAGE $SIGMA
else
  # chain job: work on the first unfinished stage; when it is done and time is left, go on with the next one
  while W=$(next_work) && [ -n "$W" ]; do
    set -- $W; S=$2; [ "$S" = - ] && S=$SIGMA
    run_stage $1 $S || break
    [ -f "$(donefile $1 $S)" ] || break          # stopped by the time limit: the next job continues
  done
  W=$(next_work)
  if [ -z "$W" ]; then echo "ALL STAGES DONE"; else echo "next: stage $W"; extend_chain; fi
fi
echo "end $(date)"
