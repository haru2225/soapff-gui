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
# a job queued before job.sh was updated (PBS keeps the script of the time of qsub): run the current job.sh of the run folder instead
if [ -f ./in.01_equil ] && [ -f ./job.sh ] && [ -z "${JOB_SH_CURRENT:-}" ] && ! cmp -s "$0" ./job.sh 2>/dev/null; then
  export JOB_SH_CURRENT=1; echo "running the current job.sh of the run folder (this job was queued with an older one)"; exec bash ./job.sh "$@"
fi
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
# GPU (go.sh GPU=<GPUs per node>, in config.env): the jobs run a LAMMPS with the GPU or KOKKOS package; its own cache lmp_gpu.auto
GPUMODE=""; [ -f ./config.env ] && GPUMODE=$(sed -n 's/^GPU="\{0,1\}\([1-9][0-9]*\).*/\1/p' ./config.env | tail -1)
AUTO=lmp.auto; [ -n "$GPUMODE" ] && AUTO=lmp_gpu.auto
# lmp.auto: the LAMMPS build and the MPI launcher an earlier job of this folder found and checked (no search again; delete lmp.auto to search again)
[ -f ./$AUTO ] && source ./$AUTO
[ -f ./config.env ] && source ./config.env
if [ -n "$GPUMODE" ]; then
  # GPU_QSUB: the qsub options of the GPU jobs (also for the jobs this one submits); GPU_MODULES: modules for them (CUDA ...);
  # GPU_LMP: a LAMMPS with the GPU or KOKKOS package (else one is looked for, as for press/langevin below)
  [ -n "${GPU_QSUB:-}" ] && QSUB_OPTS=$GPU_QSUB
  [ -n "${GPU_MODULES:-}" ] && MODULES="${MODULES:+$MODULES; }$GPU_MODULES"
  SIF=""; [ -n "${GPU_LMP:-}" ] && LMP=$GPU_LMP
fi
# the wall time of the jobs (file "walltime" of the run folder, written by WALLTIME=3 ./go.sh): every walltime in the qsub options is replaced by it
set_walltime() { echo "$1" | sed -E 's/walltime=[0-9:]+,//g; s/,walltime=[0-9:]+//g; s/-l[[:space:]]*walltime=[0-9:]+//g' | awk -v w="$2" '{ if (NF) { $1 = $1; print $0 " -l walltime=" w } else print "-l walltime=" w }'; }
[ -f ./walltime ] && QSUB_OPTS=$(set_walltime "${QSUB_OPTS:-}" "$(cat ./walltime)")
# the qsub options of the last ./go.sh (submit_chain.sh: nodes, OMP, GPU, wall time applied): the next job of the chain gets the same
[ -f ./qsub.opts ] && QSUB_OPTS=$(cat ./qsub.opts)
# OpenMP (go.sh OMP=<threads per MPI rank>, config.env; not with the GPU): cores per node / OMP ranks per node, OMP threads each
OMPN=""; [ -z "$GPUMODE" ] && [ -f ./config.env ] && OMPN=$(sed -n 's/^OMP="\{0,1\}\([0-9][0-9]*\).*/\1/p' ./config.env | tail -1)
[ -n "$OMPN" ] && [ "$OMPN" -le 1 ] && OMPN=""
# a command with a time limit (a search through a large shared file system must not eat the job's wall time)
tmo() { local t=$1; shift; if command -v timeout >/dev/null 2>&1; then timeout "$t" "$@"; else "$@"; fi; }
APP=${APP_ROOT:-/home/app}                     # the centre's software area
SEARCH_ROOTS=${SEARCH_ROOTS:-"$APP /opt /usr/local $HOME /work"}   # where other LAMMPS builds and libraries are looked for
# safety: MPI_OPTS must hold launcher options only; cut at the first word that is the program name or a redirect (old config.auto files)
if [ -n "${MPI_OPTS}" ]; then
  MPI_OPTS=$(echo "${MPI_OPTS}" | awk '{o=""; for(i=1;i<=NF;i++){ if ($i ~ /(^|\/)(lmp[_a-z0-9]*|lammps[_a-z0-9]*)$/ || $i ~ /^[<>]/ || $i ~ /^\$/ ) break; o=o" "$i } sub(/^ /,"",o); print o}')
fi
if [ -n "${MODULES}" ]; then set +u; eval "${MODULES}"; set -u; fi
[ -n "${LMP_LIBDIRS}" ] && export LD_LIBRARY_PATH="${LMP_LIBDIRS}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=${OMPN:-1}
HOSTOPT=""
NPROC=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN)
if [ -n "${PBS_NODEFILE:-}" ]; then
  # cores per node = the largest of: slots listed per host in $PBS_NODEFILE, the cores this job may use (nproc), the cores requested (NCPUS)
  SLOTS=$(sort "$PBS_NODEFILE" | uniq -c | awk '{print $1}' | sort -rn | head -1)
  NREQ=${NCPUS:-${PBS_NCPUS:-0}}
  C=${CORES_PER_NODE:-$(echo "$SLOTS $NPROC $NREQ" | tr ' ' '\n' | sort -rn | head -1)}
  NODES=$(sort -u "$PBS_NODEFILE" | wc -l | tr -d " ")
  echo "cores: slots/host in PBS_NODEFILE=${SLOTS}  nproc=${NPROC}  NCPUS=${NREQ}  -> ${C} per node x ${NODES} nodes"
  # MPI ranks per node: the cores, or with OMP the cores / OMP (at most the slots PBS gives: mpiprocs)
  # with OMP: the ranks per node PBS gives (mpiprocs = ncpus / OMP, submit_chain.sh).  A job submitted before OMP was set (no ompthreads in
  # its select): MPI only (fewer ranks would crowd the first nodes).  PBS sets NCPUS to ompthreads, so C may not be the cores then.
  RPN=$C
  if [ -n "$OMPN" ]; then
    JSEL=""; [ -n "${PBS_JOBID:-}" ] && JSEL=$(qstat -f "$PBS_JOBID" 2>/dev/null | tr -d '\n\t' | grep -o 'Resource_List.select = [^ ]*' | head -1)
    case "$JSEL" in
      *ompthreads=$OMPN*) RPN=${SLOTS:-0}; [ "$RPN" -gt 0 ] || RPN=$((C / OMPN));;
      "") RPN=$((C / OMPN));;
      *) echo "note: OMP=$OMPN, but this job was submitted for MPI only (${JSEL#Resource_List.}): MPI only in this job, the next ones use OMP"; OMPN="";;
    esac
    [ "$RPN" -lt 1 ] && RPN=1
  fi
  export OMP_NUM_THREADS=${OMPN:-1}
  [ "${C}" -le 1 ] && echo "WARNING: only 1 core per node. Ask for more cores in the resource line (e.g. -l select=${NODES}:ncpus=<cores>:mpiprocs=<cores>); see QSUB_OPTS in config.auto" 
  # default: no host option (a PBS-aware mpirun reads the allocation itself). Only if config sets HOSTFILE_FMT=openmpi|mpich a host file is passed.
  HF=hostfile.${PBS_JOBID:-$$}
  if [ "${HOSTFILE_FMT}" = mpich ]; then sort -u "$PBS_NODEFILE" | awk -v c="$RPN" '{print $1":"c}' > "$HF"; HOSTOPT="-machinefile $HF"
  elif [ "${HOSTFILE_FMT}" = openmpi ]; then sort -u "$PBS_NODEFILE" | awk -v c="$RPN" '{print $1" slots="c}' > "$HF"; HOSTOPT="--hostfile $HF"; fi
  NP=${NP:-$((NODES*RPN))}
else
  NP=${NP:-$((NPROC / ${OMPN:-1}))}; [ "$NP" -lt 1 ] && NP=1
fi
if [ -n "${SIF}" ]; then LMPCMD="singularity exec ${SING_OPTS} --bind $PWD ${SIF} ${LMP_IN_IMAGE}"; else LMPCMD="${LMP}"; fi
if [ -z "${SIF}" ]; then
  # the program is given to mpirun with its absolute path (the other nodes may not have the module's PATH)
  LB=${LMP%% *}; LMPABS=$(command -v "$LB" 2>/dev/null)
  [ -z "$LMPABS" ] && [ -x "$LB" ] && LMPABS=$(readlink -f "$LB" 2>/dev/null || echo "$LB")
  if [ -z "$LMPABS" ]; then
    echo "note: '$LB' (from the centre's sample script) is not in PATH on this node; searching for it ..."
    BN=$(basename "$LB")
    # the usual places first (a find through the shared file system can be slow and was cut by its time limit), then find
    for f in "$APP"/lammps/*/bin/$BN "$APP"/*/bin/$BN "$APP"/*/*/bin/$BN /opt/*/bin/$BN /usr/local/*/bin/$BN "$HOME"/*/bin/$BN "$HOME"/bin/$BN; do
      [ -f "$f" ] && [ -x "$f" ] && { LMPABS=$f; break; }
    done
    if [ -z "$LMPABS" ]; then
      for root in $SEARCH_ROOTS; do
        [ -d "$root" ] || continue
        LMPABS=$(tmo 300 find "$root" -maxdepth 6 -type f -name "$BN" 2>/dev/null | head -1); [ -n "$LMPABS" ] && break
      done
    fi
  fi
  if [ -z "$LMPABS" ]; then
    echo "ERROR: LAMMPS executable '$LB' not found on $(hostname) (lmp.auto: $([ -f lmp.auto ] && echo there || echo missing)). Run on the login node:  which $LB;  ls /home/app/lammps/*/bin;  cat /home/app/lammps/run.sh"
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
      # the usual places first (seconds), then find (a find through the shared file system can be slow)
      for g in "$LMPDIR"/lib*/"$l"* "$(dirname "$LMPDIR")"/*/lib*/"$l"* "$APP"/*/lib*/"$l"* "$APP"/*/*/lib*/"$l"* "$APP"/*/*/*/lib*/"$l"* "$APP"/*/*/*/*/lib*/"$l"*; do
        [ -e "$g" ] && { f=$g; break; }
      done
      [ -z "$f" ] && for root in "$LMPDIR" "$(dirname "$LMPDIR")" /home/app /opt /usr/lib64 /usr/local "$HOME" /work; do
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
  # GPU: a build with the GPU or KOKKOS package (the style names are in the program or its liblammps.so)
  in_lmp() {
    LC_ALL=C grep -aqE "$2" "$1" 2>/dev/null && return 0
    local so; so=$(ldd "$1" 2>/dev/null | awk '/liblammps/ && $3 ~ /^\//{print $3; exit}')
    [ -n "$so" ] && LC_ALL=C grep -aqE "$2" "$so" 2>/dev/null
  }
  has_gpu() { in_lmp "$1" 'pppm/(gpu|kk)'; }
  ok_build() { if [ -n "$GPUMODE" ]; then has_gpu "$1"; else has_pl "$1"; fi; }
  # can it start here?  A library not found, or a build made with another MPI version (lmp_icex 3Mar20 with HPE MPT 2.33:
  # "symbol lookup error: undefined symbol: mpi_sgi_status_ignore"), stops before LAMMPS begins.  ldd -r resolves every symbol without running it.
  load_errors() { command -v ldd >/dev/null 2>&1 || return 0; ldd -r "$1" 2>&1 | grep -E 'not found|undefined symbol' | sort -u | head -5; }
  # the first program of the list with press/langevin that can start here becomes LAMMPS (missing libraries are looked for)
  pick() {
    local c CE OLDLD OLDADD
    for c in "$@"; do
      [ -f "$c" ] && [ -x "$c" ] && [ "$c" != "$LMPABS" ] || continue
      ok_build "$c" || continue
      OLDLD=${LD_LIBRARY_PATH:-}; OLDADD=$ADDED_LIBDIRS; fixlibs "$c" >/dev/null 2>&1
      CE=$(load_errors "$c")
      if [ -z "$CE" ]; then
        LMPABS=$c; LMPCMD=$c; echo "using $c (has $([ -n "$GPUMODE" ] && echo "the GPU or KOKKOS package" || echo press/langevin); the notes above are no problem)"
        [ -n "$ADDED_LIBDIRS" ] && echo "  libraries from: $ADDED_LIBDIRS"; return 0
      fi
      echo "  skipped $c:"; echo "$CE" | sed 's/^/    /'; export LD_LIBRARY_PATH="$OLDLD"; ADDED_LIBDIRS=$OLDADD
    done
    return 1
  }
  E=$(load_errors "$LMPABS")
  if ! ok_build "$LMPABS" || [ -n "$E" ]; then
    # normal on the first job with the centre's lmp_icex (3Mar20): another build is chosen here and kept in lmp.auto for the next jobs
    if [ -n "$GPUMODE" ]; then has_gpu "$LMPABS" || echo "note: $LMPABS has no GPU or KOKKOS package (GPU=$GPUMODE in config.env)"
    else has_pl "$LMPABS" || echo "note: $LMPABS has no 'fix press/langevin' (LAMMPS older than 21Nov2023)"; fi
    [ -n "$E" ] && { echo "note: $LMPABS cannot start in this environment:"; echo "$E" | sed 's/^/    /'; }
    echo "looking for another LAMMPS build (with $([ -n "$GPUMODE" ] && echo "the GPU or KOKKOS package" || echo press/langevin), starting here) ..."
    # the centre's LAMMPS folders first (/home/app/lammps/<version>/..., e.g. 29Aug2024/lammps-29Aug2024_cpu/install_cpu/bin/lmp): seconds
    G=$(ls -d "$APP"/lammps*/bin/lmp* "$APP"/lammps*/*/bin/lmp* "$APP"/lammps*/*/*/bin/lmp* "$APP"/lammps*/*/*/*/bin/lmp* "$APP"/lammps*/*/*/*/*/bin/lmp* 2>/dev/null | grep -v '\.so' | sort -r)
    if ! pick $G; then
      # then a search of the disks (slow: up to 2 min per folder), not again when an earlier job searched in vain
      if [ -z "$GPUMODE" ] && [ "${PL_SEARCHED:-}" = "$LMPABS" ]; then echo "  (an earlier job searched the disks in vain: not again; delete lmp.auto to search again)"
      else for root in $SEARCH_ROOTS; do [ -d "$root" ] || continue; pick $(tmo 120 find "$root" -maxdepth 7 -type f \( -name 'lmp' -o -name 'lmp_*' \) -perm -u+x 2>/dev/null | grep -v '\.so' | sort -r | head -30) && break; done; fi
    fi
    if [ -n "$GPUMODE" ]; then
      if ! has_gpu "$LMPABS"; then
        echo "ERROR: GPU=$GPUMODE but no LAMMPS with the GPU or KOKKOS package found (ls /home/app/lammps*; module avail lammps).  Give it to go.sh:"
        echo "  env GPU_LMP=/path/to/lmp GPU=$GPUMODE ... ./go.sh ...   (or GPU=0: back to the CPU)"
        echo "stage - job ${PBS_JOBID:-} $(date): GPU=$GPUMODE: no LAMMPS with the GPU or KOKKOS package (GPU_LMP=...)" > chain.failed; exit 1
      fi
    else
      has_pl "$LMPABS" || echo "note: no LAMMPS with press/langevin found: all stages run (the shear uses fix nph on z, flow 1); only FLOW=0 (SLLOD) needs press/langevin"
    fi
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
  # keep what was found for the next jobs (no search, no check): the full path of LAMMPS, also of a build without press/langevin
  if [ -z "$(load_errors "$LMPABS")" ]; then
    {
      echo "# written by job.sh $(date): the LAMMPS build and MPI launcher found and checked in this folder.  Delete this file to search again."
      echo "LMP=\"$LMPABS\""
      has_pl "$LMPABS" || echo "PL_SEARCHED=\"$LMPABS\""
      echo "LMP_LIBDIRS=\"${ADDED_LIBDIRS}${LMP_LIBDIRS:+${ADDED_LIBDIRS:+:}$LMP_LIBDIRS}\""
      [ -n "${MPIRUN_OK}" ] && echo "MPIRUN_FROM=\"${MPIRUN_FROM_NOW}\"; MPIRUN_OK=\"${MPIRUN_OK}\""
      [ -n "$GPUMODE" ] && [ -n "${GPUTRI:-}" ] && echo "GPUTRI=$GPUTRI"
    } > $AUTO
  fi
fi
# ---- GPU: the package of this build, the MPI ranks, and whether its pppm runs a triclinic box (the shear) on the GPU: a tiny LAMMPS run
#      tells it once (GPUTRI in lmp_gpu.auto; LAMMPS 3Mar2020: no, newer: yes); no: in.all runs pppm on the CPU for a triclinic box
GPUOPT=""; GVAR=""
if [ -n "$GPUMODE" ]; then
  GPKG=${GPU_PKG:-}
  if [ -z "$GPKG" ]; then if [ -n "${LMPABS:-}" ] && in_lmp "$LMPABS" 'pppm/gpu'; then GPKG=gpu; else GPKG=kokkos; fi; fi
  NN=${NODES:-1}
  if [ "$GPKG" = kokkos ]; then
    GPUOPT="-k on g $GPUMODE -sf kk"; G1="-k on g 1 -sf kk"; NP=$((NN * GPUMODE))                     # KOKKOS: one rank per GPU
  else
    GPUOPT="-sf gpu -pk gpu $GPUMODE"; G1="-sf gpu -pk gpu 1"                                        # GPU package: GPU_RPG ranks share a GPU
    NPG=$((NN * GPUMODE * ${GPU_RPG:-8})); [ "$NPG" -lt "$NP" ] && NP=$NPG
  fi
  gprobe() {   # $1 = 1: a triclinic box
    local R="region b block 0 24 0 24 0 24"; [ "$1" = 1 ] && R="region b prism 0 24 0 24 0 24 3 0 0"
    printf 'units real\natom_style charge\n%s\ncreate_box 2 b\ncreate_atoms 1 random 20 11 NULL\ncreate_atoms 2 random 20 12 NULL\nmass * 10.0\nset type 1 charge 1.0\nset type 2 charge -1.0\npair_style lj/cut/coul/long 8.0\npair_coeff * * 0.1 3.0\nkspace_style pppm 1.0e-4\nrun 0\n' "$R" > gprobe.in
    rm -f gprobe.out; tmo 300 ${MPIRUN} -np 1 ${LMPCMD} ${G1} -in gprobe.in -log none -screen gprobe.out >/dev/null 2>&1 </dev/null
  }
  if [ -z "${GPUTRI:-}" ]; then
    if ! gprobe 0; then
      echo "ERROR: LAMMPS does not run on the GPU here (${LMPCMD} ${G1}):"; grep -E "ERROR|rror" gprobe.out 2>/dev/null | head -5; tail -5 gprobe.out 2>/dev/null
      echo "  GPU_MODULES='module load cuda' (the modules of the GPU build), GPU_LMP=..., GPU_PKG=gpu|kokkos to go.sh; GPU=0: back to the CPU"
      echo "stage - job ${PBS_JOBID:-} $(date): LAMMPS does not run on the GPU: $(grep -m1 ERROR gprobe.out 2>/dev/null)" > chain.failed; exit 1
    fi
    if gprobe 1; then GPUTRI=1; else GPUTRI=0; fi
    echo "GPUTRI=$GPUTRI   # pppm of a triclinic box on the GPU (1) or the CPU (0)" >> $AUTO
  fi
  rm -f gprobe.in gprobe.out
  GVAR="-var gpu 1 -var gputri $GPUTRI"
  echo "GPU: package $GPKG, $GPUMODE GPU(s) per node x $NN node(s), MPI $NP ranks; pppm of the triclinic box (shear) on the $([ "$GPUTRI" = 1 ] && echo GPU || echo CPU)"
fi
# ---- OpenMP (OMP): LAMMPS with its OPENMP package (-sf omp: the /omp styles, OMP threads per rank); HPE MPT places the threads of each
#      rank on cores of their own with omplace (without it, all threads of a rank share one core); in.all: plain pppm for a triclinic box
OMPOPT=""; OMPW=""; OVAR=""
if [ -n "$OMPN" ]; then
  if [ -z "${SIF}" ] && ! in_lmp "$LMPABS" 'lj/cut/coul/long/omp'; then
    echo "ERROR: OMP=$OMPN, but $LMPABS has no OPENMP package (no lj/cut/coul/long/omp).  env OMP=1 ./go.sh (the same command): MPI only"
    echo "stage - job ${PBS_JOBID:-} $(date): OMP=$OMPN, but LAMMPS has no OPENMP package; env OMP=1 ./go.sh (the same command) runs MPI only" > chain.failed; exit 1
  fi
  OMPOPT="-sf omp -pk omp $OMPN"; OVAR="-var omp 1"
  command -v omplace >/dev/null 2>&1 && OMPW="omplace -nt $OMPN"
  export I_MPI_PIN_DOMAIN=omp
  echo "OpenMP: MPI ${NP} ranks (${RPN:-?} per node) x ${OMPN} threads$([ -n "$OMPW" ] && echo ", placed with omplace")"
fi
# the normal stress of the shear (04): fix press/langevin if this LAMMPS has it, else none (-var pl 0: Lz of stage 03 kept; EXTRA="-var pl 1" overrides)
PLVAR=""
if [ -z "${SIF}" ] && ! has_pl "$LMPABS"; then PLVAR="-var pl 0"; echo "note: $LMPABS has no press/langevin: stage 04 shears at the Lz of stage 03 (constant volume; the normal stress is written in shear_<s>.dat)"; fi
# the screen output of rank 0 goes to lmp.screen (each LAMMPS start overwrites it): an error that only one rank sees ("ERROR on proc N",
# e.g. Out of range atoms, Bond atoms missing) is not in the log file, only on the screen (rank 0: lmp.screen; other ranks: this job's output)
RUN="${MPIRUN} -np ${NP} ${MPI_OPTS} ${HOSTOPT} ${OMPW} ${LMPCMD} ${GPUOPT} ${OMPOPT} -screen lmp.screen"
SIGMA=${SIGMA:-5.6}; EXTRA="${EXTRA:-} ${GVAR} ${OVAR}"
[ -f ./chain.env ] && source ./chain.env          # SIGMAS and STAGES of the chain (written by submit_chain.sh)
[ -f ./lammps_seed ] && EXTRA="-var seed $(cat ./lammps_seed) ${EXTRA}"   # an independent run (go.sh, 6th argument)
SIGMAS=${SIGMAS:-$SIGMA}; STAGES=${STAGES:-all}
echo "MPI ${NP} (${RPN:-${C:-?}} per node)${OMPN:+ x ${OMPN} OpenMP threads}  $(date)"
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

donefile() { case $1 in one) echo all.done;; water) echo water.done;; 01) echo 01_equil.done;; 02) echo 02_dehydrate.done;; 03) echo 03_consolidate_$2.done;; 04) echo shear_$2.done;; 05) echo heal_$2.done;; esac; }
want() { [ "$STAGES" = all ] || [[ " $STAGES " == *" $1 "* ]]; }
next_work() {   # "STAGE SIGMA" of the first unfinished stage: 01, 02, then 03 04 05 for each sigma in turn; nothing when all is done
  local st S
  # STAGES=one: everything in one input file (in.all), every run nrun steps
  if [ "$STAGES" = one ]; then [ -f all.done ] || echo "one -"; return; fi
  if [ "$STAGES" = water ]; then [ -f water.done ] || echo "water -"; return; fi
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
  local V="-var tlimit $TL ${EXTRA} ${PLVAR}"
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
  one|water)   # in.all (one) / in.water (water: a series of water contents): every stage in one input file; it goes on from <name>.state itself.
      # With a time limit it does one run per start: started again while time is left, until all is done or the time limit stopped a run
      local IN=in.all STF=all.state; [ $st = water ] && { IN=in.water; STF=water.state; }
      set -- $SIGMAS
      [ -n "${3:-}" ] && echo "note: $IN takes two normal stresses; $SIGMAS: only $1 and $2"
      rc=0
      while :; do
        RES=0; [ -f $STF ] && RES=1
        ${RUN} -in $IN -var data gouge.data -var settings gouge.settings -var sig1 $1 -var sig2 ${2:-0} -var resume $RES -var tlimit $TL ${EXTRA} ${PLVAR}; rc=$?
        { [ $rc -ne 0 ] || [ -f "$(donefile $st)" ] || grep -q "mid equal 1" $STF 2>/dev/null; } && break
        TL=$(tleft)
        if [ "$TL" != off ] && [ "$TL" -lt "${MIN_LEFT:-180}" ]; then echo "only ${TL} s left: the next job goes on ($(head -1 $STF 2>/dev/null))"; break; fi
        echo "-- $IN: next run $(head -1 $STF 2>/dev/null | awk '{print $4}')  $(date)  LAMMPS time limit ${TL} s"
      done
      (exit $rc) ;;
  *)  echo "unknown STAGE $st"; exit 1 ;;
  esac
  rc=$?
  [ -f procs.txt ] && echo "MPI ranks per node: $(awk -F': ' '$1 ~ /^[0-9]/{print $3}' procs.txt | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"
  # how far this job got: the last thermo line of the stage (step ... last column = time steps per second)
  local L; L=$(ls -t log.${st}* 2>/dev/null | head -1)
  [ -n "$L" ] && echo "LAMMPS ended $(date) with status $rc.  last thermo line of $L (step temp pe ... steps/s):" && grep -E "^ +[0-9]+ +[0-9.-]+ " "$L" | tail -1 | awk '{print "   step", $1, "  T", $2, "  steps/s", $NF}'
  [ $rc -eq 0 ] && rm -f .killed
  if [ $rc -eq 137 ] || [ $rc -eq 143 ] || [ $rc -eq 134 ]; then
    # killed from outside (137 = SIGKILL, 143 = SIGTERM, 134 = abort): mostly the memory of a node (the kernel kills the biggest process), else the
    # node or MPI.  Not an error of the input: the next job goes on from the last restart; three kills in a row stop the chain.
    local n; n=$(( $(cat .killed 2>/dev/null || echo 0) + 1 )); echo $n > .killed
    echo "LAMMPS was killed from outside (status $rc = signal $((rc - 128))): out of memory on a node, or the node / MPI.  Kill $n in a row."
    echo "  memory: $(qstat -f ${PBS_JOBID:-none} 2>/dev/null | grep -E 'resources_used.mem|resources_used.vmem|Resource_List.mem|Resource_List.select' | tr -s ' ' | tr '\n' ' ')"
    echo "  node: $(free -g 2>/dev/null | awk '/^Mem/{print "total " $2 " GB, used " $3 " GB, free " $4 " GB"}')"
    for f in $(ls -t log.${st}* 2>/dev/null | head -1); do echo "--- last lines of $f"; grep -E "^(in\.(all|water)|ERROR|WARNING)" "$f" | tail -3; tail -2 "$f"; done
    if [ $n -lt 3 ]; then echo "the next job goes on from the last restart"; return 1; fi
  fi
  if [ $rc -ne 0 ]; then
    # LAMMPS stopped with an error (a stop by the time limit ends with status 0): do not let the following jobs repeat it
    local EM; EM=$( { grep -h -E "^ERROR|Last command" lmp.screen 2>/dev/null; tail -400 jobs.log 2>/dev/null | grep -E "ERROR on proc|MPT ERROR|MPI_ABORT|Segmentation|signal [0-9]|Killed|Bus error|out of memory|oom" | grep -v "^ERROR: LAMMPS ended"; } | tail -4 | cut -c1-300)
    { echo "stage $st sigma $S job ${PBS_JOBID:-} $(date): LAMMPS ended with status $rc"; [ -n "$EM" ] && echo "$EM"; } > chain.failed
    echo "ERROR: LAMMPS ended with status $rc.  The chain stops here (chain.failed); the following jobs end at once."
    echo "--- the error (lmp.screen = the screen of rank 0; the other ranks print to this output):"; echo "${EM:-  (no error message found: see lmp.screen and the lines above)}"
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
  if [ -f chain.left ] && [ "$(cat chain.left 2>/dev/null || echo 0)" -le 0 ]; then
    echo "MAX_JOBS / MAX_HOURS reached: no more jobs.  ./go.sh again (the same command) continues from the last restart"; return 0
  fi
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
# MAX_JOBS / MAX_HOURS (submit_chain.sh: chain.left = the jobs left): this job counts; none left: it ends at once
if [ -f chain.left ]; then
  L=$(cat chain.left 2>/dev/null || echo 0)
  if [ "${L:-0}" -le 0 ]; then echo "MAX_JOBS / MAX_HOURS used up: this job ends.  ./go.sh again (the same command) continues from the last restart"; exit 0; fi
  echo $((L - 1)) > chain.left; echo "MAX_JOBS / MAX_HOURS: $((L - 1)) job(s) left after this one (chain.left)"
fi

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
# the figures of the shear (hot_spots.py, copied here by submit_chain.sh) in figs/ after every job that got to the shear (stage 04), in the
# time left before the wall time (a trajectory frame for fig9 is analysed once and kept in figs/.cache, so the next job goes on)
if [ -f hot_spots.py ] && ls shear_*.dat >/dev/null 2>&1; then
  PYF=${PYTHON:-python3}
  FT=$(( WALLSEC > 0 ? WALLSEC - SECONDS - 60 : 900 ))
  if [ "$FT" -lt 60 ]; then echo "figures: no time left in this job (the next one makes them)"
  elif $PYF -c "import numpy, matplotlib" >/dev/null 2>&1; then
    mkdir -p figs
    if tmo "$FT" $PYF hot_spots.py "$PWD" --outdir figs --budget "$FT" > figs/hot_spots.log 2>&1; then echo "figures: $PWD/figs  ($(tail -1 figs/hot_spots.log))"
    else echo "note: hot_spots.py did not finish (figs/hot_spots.log): $(tail -2 figs/hot_spots.log | tr '\n' ' ')"; fi
  else echo "note: $PYF on this node has no numpy/matplotlib: no figures here.  On the login node: cd $(dirname $PWD); ./figs.sh $(basename $PWD)"; fi
fi
# a series of water contents (in.water): figures of the kind of Lin & Hatano (2024) against the porosity (lin_figs.py, figs/lin_*), then the hot-spot
# figures of each water content (hot_spots.py, w<m_w/m_c>/figs) in the time left
if [ "$STAGES" = water ] && [ -f lin_figs.py ] && [ -f hot_spots.py ] && ls w*/shear_*.dat >/dev/null 2>&1; then
  PYF=${PYTHON:-python3}
  FT=$(( WALLSEC > 0 ? WALLSEC - SECONDS - 60 : 900 ))
  if [ "$FT" -lt 60 ]; then echo "figures: no time left in this job (the next one makes them)"
  elif $PYF -c "import numpy, matplotlib" >/dev/null 2>&1; then
    mkdir -p figs
    if tmo "$FT" $PYF lin_figs.py "$PWD" --outdir figs > figs/lin_figs.log 2>&1; then echo "figures: $PWD/figs/lin_*.png"
    else echo "note: lin_figs.py did not finish (figs/lin_figs.log): $(tail -2 figs/lin_figs.log | tr '\n' ' ')"; fi
    for w in w*/; do
      w=${w%/}; ls $w/shear_*.dat >/dev/null 2>&1 || continue
      FT=$(( WALLSEC > 0 ? WALLSEC - SECONDS - 60 : 900 )); [ "$FT" -lt 60 ] && break
      for f in gouge.data config.env; do [ -e $w/$f ] || ln -s ../$f $w/$f; done
      mkdir -p $w/figs
      tmo "$FT" $PYF hot_spots.py "$PWD/$w" --outdir $w/figs --budget "$FT" > $w/figs/hot_spots.log 2>&1 && echo "figures: $PWD/$w/figs" \
        || echo "note: hot_spots.py $w did not finish ($w/figs/hot_spots.log)"
    done
  else echo "note: $PYF on this node has no numpy/matplotlib: no figures here.  On the login node: cd $(dirname $PWD); ./figs.sh $(basename $PWD)"; fi
fi
echo "end $(date)"
