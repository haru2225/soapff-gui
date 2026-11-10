#!/bin/bash
#PBS -j oe
# One job of the chain.  Everything is passed by submit_chain.sh:  STAGE (01..05)  SIGMA (normal stress, MPa)  [EXTRA]
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
QSUB_OPTS=""; MODULES=""; SIF=""; SING_OPTS=""; LMP="lmp"; LMP_IN_IMAGE="lmp"; MPIRUN="mpirun"; MPI_OPTS=""; HOSTFILE_FMT="none"; CORES_PER_NODE="${CORES_PER_NODE:-}"
[ -f ./config.auto ] && source ./config.auto
[ -f ./config.env ] && source ./config.env
# safety: MPI_OPTS must hold launcher options only; cut at the first word that is the program name or a redirect (old config.auto files)
if [ -n "${MPI_OPTS}" ]; then
  MPI_OPTS=$(echo "${MPI_OPTS}" | awk '{o=""; for(i=1;i<=NF;i++){ if ($i ~ /(^|\/)(lmp[_a-z0-9]*|lammps[_a-z0-9]*)$/ || $i ~ /^[<>]/ || $i ~ /^\$/ ) break; o=o" "$i } sub(/^ /,"",o); print o}')
fi
if [ -n "${MODULES}" ]; then set +u; eval "${MODULES}"; set -u; fi
export OMP_NUM_THREADS=1
HOSTOPT=""
NPROC=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN)
if [ -n "${PBS_NODEFILE:-}" ]; then
  # cores per node = the largest of: slots listed per host in $PBS_NODEFILE, the cores this job may use (nproc), the cores requested (NCPUS)
  SLOTS=$(sort "$PBS_NODEFILE" | uniq -c | awk '{print $1}' | sort -rn | head -1)
  NREQ=${NCPUS:-${PBS_NCPUS:-0}}
  C=${CORES_PER_NODE:-$(echo "$SLOTS $NPROC $NREQ" | tr ' ' '\n' | sort -rn | head -1)}
  NODES=$(sort -u "$PBS_NODEFILE" | wc -l)
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
    echo "WARNING: '$LB' is not in PATH here. PATH=${PATH}"; echo "searching for it ..."
    LMPABS=$(find /home/app /opt /usr/local /work "$HOME" -maxdepth 6 -type f -name "$(basename "$LB")" 2>/dev/null | head -1)
  fi
  if [ -z "$LMPABS" ]; then
    echo "ERROR: LAMMPS executable '$LB' not found. Run on the login node:  which $LB;  module avail 2>&1 | grep -i lammps;  cat /home/app/lammps/run.sh"
    echo "then put the module/PATH lines in config.env (MODULES=...) or the full path (LMP=/full/path/$LB)"; exit 1
  fi
  echo "LAMMPS: $LMPABS"; LMPCMD="$LMPABS"
  # shared libraries a binary needs but cannot find (e.g. libfftw3.so.3): look for them in the software area and add their directories
  fixlibs() {
    local BIN=$1 MISS LMPDIR l f root
    MISS=$(ldd "$BIN" 2>/dev/null | awk '/not found/{print $1}'); [ -z "$MISS" ] && return 0
    echo "missing shared libraries of $BIN: $(echo $MISS)"
    LMPDIR=$(dirname "$(dirname "$BIN")")
    for l in $MISS; do
      f=""
      for root in "$LMPDIR" "$(dirname "$LMPDIR")" /home/app /opt /usr/lib64 /usr/local /work "$HOME"; do
        [ -d "$root" ] && f=$(find "$root" -maxdepth 8 -name "${l}*" 2>/dev/null | head -1) && [ -n "$f" ] && break
      done
      if [ -n "$f" ]; then export LD_LIBRARY_PATH="$(dirname "$f")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"; echo "  $l -> $(dirname "$f")"
      else echo "  $l: not found anywhere; try 'module avail 2>&1 | grep -i ${l#lib}' and add the module to config.env (MODULES=...)"; fi
    done
    ldd "$BIN" 2>/dev/null | grep "not found" && echo "WARNING: still missing libraries (see above)"
    return 0
  }
  fixlibs "$LMPABS"
  # the shear stage needs "fix press/langevin" (LAMMPS 21Nov2023 or newer): if this build lacks it, use another LAMMPS found in the software area
  has_pl() { "$1" -h 2>/dev/null | grep -q "press/langevin"; }
  if ! has_pl "$LMPABS"; then
    echo "WARNING: $LMPABS has no 'fix press/langevin' (LAMMPS older than 21Nov2023); looking for a newer build ..."
    for c in $(find /home/app /opt /work "$HOME" -maxdepth 7 -type f \( -name 'lmp' -o -name 'lmp_*' \) -perm -u+x 2>/dev/null | grep -v '\.so' | sort -r | head -20); do
      fixlibs "$c" >/dev/null 2>&1
      if has_pl "$c"; then LMPABS=$c; LMPCMD=$c; echo "using $c (has press/langevin)"; break; fi
    done
    has_pl "$LMPABS" || echo "WARNING: no LAMMPS with press/langevin found; stages 01-03 run, stage 04 (shear) will fail. Put LMP=/path/to/newer/lmp in config.env"
  fi
fi
RUN="${MPIRUN} -np ${NP} ${MPI_OPTS} ${HOSTOPT} ${LMPCMD} -screen none"
SIGMA=${SIGMA:-5.6}; EXTRA=${EXTRA:-}
echo "stage ${STAGE}  sigma ${SIGMA}  MPI ${NP} (${C:-?} per node)  $(date)"
echo "${RUN}"
case "${STAGE}" in
01) [ -f 01_equil.done ] && exit 0
    ${RUN} -in in.01_equil -var data gouge.data -var settings gouge.settings ${EXTRA} ;;
02) [ -f 02_dehydrate.done ] && exit 0
    if [ -f dehyd_current.data ] && [ -f dehyd_k.txt ]; then
      ${RUN} -in in.02_dehydrate -var data dehyd_current.data -var kstart "$(cat dehyd_k.txt)" -var settings gouge.settings ${EXTRA}
    else
      ${RUN} -in in.02_dehydrate -var data 01_equil.data -var settings gouge.settings ${EXTRA}
    fi ;;
03) [ -f 03_consolidate_${SIGMA}.done ] && exit 0
    ${RUN} -in in.03_consolidate -var data 02_dehydrate.data -var settings gouge.settings -var sigma ${SIGMA} ${EXTRA} ;;
04) [ -f shear_${SIGMA}.done ] && exit 0
    R=$(ls -t shear_${SIGMA}.restart.1 shear_${SIGMA}.restart.2 2>/dev/null | head -1)
    if [ -n "${R}" ]; then
      ${RUN} -in in.04_shear -var resume 1 -var rfile "${R}" -var sigma ${SIGMA} ${EXTRA}
    else
      ${RUN} -in in.04_shear -var data 03_consolidate_${SIGMA}.data -var settings gouge.settings -var sigma ${SIGMA} ${EXTRA}
    fi ;;
05) [ -f heal_${SIGMA}.done ] && exit 0
    ${RUN} -in in.05_heal -var data shear_${SIGMA}_final.data -var settings gouge.settings -var sigma ${SIGMA} ${EXTRA} ;;
*)  echo "unknown STAGE ${STAGE}"; exit 1 ;;
esac
echo "end $(date)"
