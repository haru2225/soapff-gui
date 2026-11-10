#!/bin/bash
#PBS -j oe
# One job of the chain.  Everything is passed by submit_chain.sh:  STAGE (01..05)  SIGMA (normal stress, MPa)  [EXTRA]
# It runs in the submit directory (the run folder holding the in.* files, gouge.data and gouge.settings).
set -u
cd "${PBS_O_WORKDIR:-$PWD}"
QSUB_OPTS=""; MODULES=""; SIF=""; SING_OPTS=""; LMP="lmp"; LMP_IN_IMAGE="lmp"; MPIRUN="mpirun"; MPI_OPTS=""; HOSTFILE_FMT="none"; CORES_PER_NODE="${CORES_PER_NODE:-}"
[ -f ./config.auto ] && source ./config.auto
[ -f ./config.env ] && source ./config.env
[ -n "${MODULES}" ] && eval "${MODULES}"
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
  if ! command -v "${LMP%% *}" >/dev/null 2>&1 && [ ! -x "${LMP%% *}" ]; then
    echo "WARNING: ${LMP} not found in PATH (${PATH}); searching /home/app ..."
    F=$(find /home/app /opt /usr/local -maxdepth 5 -type f -name "$(basename ${LMP%% *})" 2>/dev/null | head -1)
    if [ -n "$F" ]; then echo "found $F"; LMPCMD="$F"; else echo "ERROR: LAMMPS executable '${LMP}' not found. Put the module/PATH lines of the sample script in config.env (MODULES=...)"; exit 1; fi
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
