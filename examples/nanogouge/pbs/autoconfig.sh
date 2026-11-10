#!/bin/bash
# Usage: autoconfig.sh RUNDIR      writes RUNDIR/config.auto from the centre's sample script (nothing to edit by hand).
#   source of the settings: $CONFIG_SRC, else /home/app/lammps/run.sh, else ~/run.sh / ./run.sh of the cwd if it has "lmp" in it.
#   taken from it:  every "#PBS" line (queue, resources, walltime; not name/output/join/mail) -> QSUB_OPTS
#                   "module ..." lines -> MODULES;   the "singularity ... <image>.sif <lmp>" line -> SIF, SING_OPTS, LMP_IN_IMAGE
#                   the mpirun / mpiexec word -> MPIRUN (mpiexec.hydra / Intel MPI -> HOSTFILE_FMT=mpich)
#   NODES=<n> ./go.sh   rescales "select=<n>:" (PBS Pro) or "nodes=<n>" in the resources.
RUNDIR=${1:-.}; mkdir -p "$RUNDIR"
SRC=${CONFIG_SRC:-}
for c in "$SRC" /home/app/lammps/run.sh "$HOME/run.sh" ./run.sh; do [ -n "$c" ] && [ -f "$c" ] && grep -q "lmp\|lammps" "$c" && { SRC=$c; break; }; done
OUT="$RUNDIR/config.auto"; : > "$OUT"
if [ -z "$SRC" ] || [ ! -f "$SRC" ]; then
  echo "autoconfig: no sample script found (set CONFIG_SRC=/path/run.sh); defaults are used" >&2
  {
    echo 'QSUB_OPTS=""'
    # default resources: one node with all cores, if the batch system tells us
    NC=$(pbsnodes -a 2>/dev/null | awk -F= '/resources_available.ncpus/{gsub(/ /,"",$2); print $2}' | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
    [ -n "$NC" ] && echo "QSUB_OPTS=\"-l select=${NODES:-1}:ncpus=$NC:mpiprocs=$NC -l walltime=24:00:00\""
    if command -v singularity >/dev/null 2>&1; then
      SIF=$(ls /home/app/*/*.sif /home/app/*/*/*.sif "$HOME"/*.sif "$HOME"/*/*.sif 2>/dev/null | grep -i lammps | head -1)
      [ -n "$SIF" ] && echo "SIF=\"$SIF\""
    fi
  } >> "$OUT"
else
  echo "autoconfig: reading $SRC" >&2
  PBS=$(grep -E '^[[:space:]]*#PBS' "$SRC" | sed -E 's/^[[:space:]]*#PBS[[:space:]]*//' | grep -vE '^-(N|o|e|j|m|M|S|V|v|J|k|r)( |$)' | sed -E 's/[[:space:]]+#.*$//')
  if [ -n "${NODES:-}" ]; then PBS=$(echo "$PBS" | sed -E "s/select=[0-9]+/select=${NODES}/; s/nodes=[0-9]+/nodes=${NODES}/"); fi
  Q=""; while IFS= read -r l; do [ -n "$l" ] && Q="$Q $l"; done <<< "$PBS"
  echo "QSUB_OPTS=\"$Q\"" >> "$OUT"
  # environment of the sample script: module / source / export / ulimit / PATH lines (not PBS directives, not the job's own cd/INPUT lines)
  MOD=$(grep -vE '^[[:space:]]*(#|$)' "$SRC" | grep -E '^[[:space:]]*(module|source|\.|export|ulimit|unset|set|PATH=|LD_LIBRARY_PATH=|MPI[A-Z_]*=|OMP[A-Z_]*=)[[:space:]=]' | grep -vE 'PBS_O_WORKDIR|INPUT|OMP_NUM_THREADS' | sed -E 's/^[[:space:]]+//' | paste -sd';' -)
  [ -n "$MOD" ] && echo "MODULES=\"$MOD\"" >> "$OUT"
  SL=$(grep -E 'singularity[[:space:]]+(exec|run)' "$SRC" | grep -vE '^[[:space:]]*#' | head -1)
  if [ -n "$SL" ]; then
    SIF=$(echo "$SL" | tr ' ' '\n' | grep -E '\.(sif|simg|img)$' | head -1)
    ex=$(echo "$SL" | sed -E 's/.*singularity[[:space:]]+(exec|run)[[:space:]]+//; s/'"$(echo "$SIF" | sed 's/[][\.*^$/]/\\&/g')"'.*//')
    BIN=$(echo "$SL" | sed -E 's#.*\.(sif|simg|img)[[:space:]]+##' | awk '{print $1}')
    [ -n "$SIF" ] && echo "SIF=\"$SIF\"" >> "$OUT"
    [ -n "$ex" ] && echo "SING_OPTS=\"$(echo $ex | sed -E 's/--bind[ =][^ ]+//g')\"" >> "$OUT"
    [ -n "$BIN" ] && echo "LMP_IN_IMAGE=\"$BIN\"" >> "$OUT"
  else
    LL=$(grep -E '(^|[[:space:]/])(lmp[_a-z0-9]*|lammps)([[:space:]]|$)' "$SRC" | grep -vE '^[[:space:]]*#' | head -1)
    B=$(echo "$LL" | tr ' ' '\n' | grep -E '(^|/)(lmp[_a-z0-9]*|lammps)$' | head -1); [ -n "$B" ] && echo "LMP=\"$B\"" >> "$OUT"
  fi
  LL=$(grep -vE '^[[:space:]]*#' "$SRC" | grep -E '(lmp[_a-z0-9]*|lammps[_a-z0-9]*)([[:space:]]|$)' | grep -E '(-in|-i)[[:space:]]|INPUT' | head -1)
  MP=$(echo "$LL" | tr ' ' '\n' | grep -E '^(mpirun|mpiexec[_a-z.]*|mpijob|srun|aprun|mpirun_mpt|mpiexec_mpt)$' | head -1)
  [ -n "$MP" ] && echo "MPIRUN=\"$MP\"" >> "$OUT"
  # options between the launcher and the program, without -np/-n N (the number of ranks is set from the allocation)
  if [ -n "$MP" ]; then
    OPTS=$(echo "$LL" | sed -E "s/.*(^|[[:space:]])$MP[[:space:]]+//" | sed -E 's/(singularity[[:space:]]+(exec|run)).*//; s/[[:space:]](lmp[_a-z0-9]*|lammps[_a-z0-9]*|[^ ]*\/lmp[_a-z0-9]*)([[:space:]].*)?$//' | sed -E 's/(-np|-n|-c)[[:space:]]+[^ ]+//g' | sed -E 's/[[:space:]]+/ /g; s/^ //')
    [ -n "$OPTS" ] && echo "MPI_OPTS=\"$OPTS\"" >> "$OUT"
  fi
  # a host file option is passed only if the sample script itself uses one (otherwise the PBS-aware mpirun reads the allocation)
  grep -qE '(-machinefile|-f)[[:space:]]+[^ ]*(NODEFILE|hostfile)' "$SRC" && echo 'HOSTFILE_FMT="mpich"' >> "$OUT"
  grep -qE '(--hostfile|-hostfile)[[:space:]]+' "$SRC" && echo 'HOSTFILE_FMT="openmpi"' >> "$OUT"
fi
grep -q "ncpus=1\b\|ncpus=1:" "$OUT" && echo "autoconfig WARNING: the sample script asks for ncpus=1; use NODES / edit QSUB_OPTS in config.auto to ask for all cores" >&2
echo "---- detected settings ($OUT) ----" >&2; cat "$OUT" >&2; echo "-----------------------------------" >&2
