#!/bin/bash
#PBS -S /bin/bash
# The figures from the results so far (hot_spots.py: Zhou et al. figs 6-10; a water series also lin_figs.py: Lin & Hatano figs against the
# porosity), by hand on the login node:
#   ./figs.sh                                  every run folder RUN_*/ with a shear (shear_*.dat, or w*/shear_*.dat of go_water.sh)
#   ./figs.sh RUN_N22_one RUN_N22_water        these folders
#   env JOB=1 ./figs.sh RUN_N22_one            the same as a short PBS job (when the login node has no singularity or limits the time)
#   env BUDGET=1800 ./figs.sh ...              s for the water analysis of the trajectory (fig9; default 600, 0 = no limit); a frame done is
#                                              kept in figs/.cache, so the next ./figs.sh goes on
# Also stacks.py: are the particles stacked (01, 02, 03, 04)?  RUN_*/figs/stacks.png, stacks.csv
# Out: RUN_N22_one/figs/*.png + hot_spots_summary.txt;  RUN_N22_water/figs/lin_*.png and RUN_N22_water/w<m_w/m_c>/figs/*.png
# Python with numpy and matplotlib (scipy for fig9): ./sifpy (Singularity, make_sif.sh: made now when it is missing), else any python of this
# machine that has them (PATH, modules, /home/app, other .sif images), else the python3 of this node + matplotlib from pip in pyfigs_lib/.
# The one found is written into the run folders' config.env (PYTHON=...): their next jobs draw the figures themselves.
# What was tried is in figs_python.log.
set -u
if [ -n "${PBS_JOBID:-}" ] && [ -n "${FIGS_DIR:-}" ]; then cd "$FIGS_DIR"; else cd "$(dirname "$0")"; fi
HERE=$(pwd -P)
tmo() { local t=$1; shift; if command -v timeout >/dev/null 2>&1; then timeout "$t" "$@"; else "$@"; fi; }
BUDGET=${BUDGET:-600}

# ---- the run folders
if [ -n "${PBS_JOBID:-}" ] && [ -f .figs_args ]; then set -- $(cat .figs_args); fi
RUNS=""
if [ $# -eq 0 ]; then
  for r in RUN_*/; do r=${r%/}; { ls $r/shear_*.dat || ls $r/w*/shear_*.dat || ls $r/0*.dump; } >/dev/null 2>&1 && RUNS="$RUNS $r"; done
  [ -n "$RUNS" ] || { echo "figs.sh: no run folder RUN_*/ here has results yet (*.dump, shear_*.dat)"; exit 1; }
else
  for r in "$@"; do
    r=${r%/}; case "$r" in "$HERE"/*) r=${r#$HERE/};; esac
    [ -d "$r" ] || { echo "figs.sh: no folder $r in $HERE"; exit 1; }
    { ls $r/shear_*.dat || ls $r/w*/shear_*.dat || ls $r/0*.dump; } >/dev/null 2>&1 || { echo "figs.sh: $r has no results yet (*.dump, shear_*.dat)"; exit 1; }
    RUNS="$RUNS $r"
  done
fi
echo "figs.sh: run folders:$RUNS"

# ---- as a PBS job: the qsub options of the first run folder (queue, nodes), one node, FIGS_WALLTIME (default 1 h)
if [ "${JOB:-0}" = 1 ] && [ -z "${PBS_JOBID:-}" ]; then
  R=$(echo $RUNS | awk '{print $1}')
  QSUB_OPTS=""
  QSUB_OPTS=$( cd "$R" && { [ -f config.auto ] && . ./config.auto; [ -f config.env ] && . ./config.env; } >/dev/null 2>&1; echo "${QSUB_OPTS:-}" )
  QO=$(echo "$QSUB_OPTS" | sed -E 's/walltime=[0-9:]+,//g; s/,walltime=[0-9:]+//g; s/-l[[:space:]]*walltime=[0-9:]+//g; s/select=[0-9]+/select=1/')
  echo $RUNS > .figs_args
  J=$(qsub -N ng_figs $QO -l walltime=${FIGS_WALLTIME:-1:00:00} -j oe -o "$HERE/figs_job.log" -v "FIGS_DIR=$HERE,BUDGET=${BUDGET}" figs.sh 2>&1) \
    || { echo "qsub refused it: $J   (qsub options: $QO)"; exit 1; }
  echo "submitted $J (qstat -u \$USER); its output: $HERE/figs_job.log"
  exit 0
fi

# ---- a python with numpy and matplotlib
LOG=$HERE/figs_python.log; : > "$LOG"
PYTEST="import numpy, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot, matplotlib.colors, matplotlib.ticker"
py_ok() { tmo 120 "$@" -c "$PYTEST" >>"$LOG" 2>&1 </dev/null; }
init_module() {
  type module >/dev/null 2>&1 && return 0
  local f
  for f in "${MODULESHOME:-/nonexistent}/init/bash" "${LMOD_PKG:-/nonexistent}/init/bash" /etc/profile.d/modules.sh /etc/profile.d/lmod.sh \
           /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash /opt/modules/default/init/bash; do
    [ -f "$f" ] && { set +u; . "$f" >/dev/null 2>&1; set -u; type module >/dev/null 2>&1 && return 0; }
  done
  return 1
}
# a wrapper ./pyfigs for a python of this machine (its absolute path and libraries; a module is not needed in the jobs then)
wrap() {  # python [PYTHONPATH]
  { echo "#!/bin/bash"; echo "# python with numpy, matplotlib for the figures (written by figs.sh)"
    echo "export LD_LIBRARY_PATH=\"${LD_LIBRARY_PATH:-}\"${2:+ PYTHONPATH=\"$2\"} MPLCONFIGDIR=\"$HERE/.mplconfig\""
    [ "${NUS:-0}" = 1 ] && echo "export PYTHONNOUSERSITE=1   # not the packages in ~/.local (a matplotlib there that does not fit the numpy)"
    echo "exec \"$1\" \"\$@\""; } > pyfigs
  chmod +x pyfigs
}
PY=""
find_py() {
  local p m s
  echo "== PYTHON / ./sifpy / ./pyfigs" >>"$LOG"
  [ -n "${PYTHON:-}" ] && py_ok "$PYTHON" && { PY=$PYTHON; return 0; }
  [ -x sifpy ] && py_ok ./sifpy && { PY=$HERE/sifpy; return 0; }
  [ -x pyfigs ] && py_ok ./pyfigs && { PY=$HERE/pyfigs; return 0; }
  # 1. Singularity: python 3.11 + numpy, scipy, matplotlib (make_sif.sh)
  echo "-- 1. Singularity image (./make_sif.sh) ..."
  echo "== make_sif.sh" >>"$LOG"
  ./make_sif.sh 2>&1 | tee -a "$LOG" | sed 's/^/   /'
  [ -x sifpy ] && py_ok ./sifpy && { PY=$HERE/sifpy; return 0; }
  # 2. a python of this machine that has numpy and matplotlib
  echo "-- 2. a python of this machine with numpy and matplotlib ..."
  echo "== pythons of this machine" >>"$LOG"
  for p in $(command -v python3 python 2>/dev/null) $(ls -d /home/app/*/bin/python3 /home/app/*/*/bin/python3 /home/app/*/*/*/bin/python3 \
             /opt/*/bin/python3 /opt/*/*/bin/python3 /usr/local/bin/python3 /usr/bin/python3.* "$HOME"/*conda*/bin/python3 "$HOME"/.local/bin/python3 2>/dev/null | grep -v -- '-config'); do
    echo "   $p" >>"$LOG"
    PYTHONNOUSERSITE=1 py_ok "$p" && { NUS=1 wrap "$p"; PY=$HERE/pyfigs; echo "   found $p"; return 0; }
    py_ok "$p" && { wrap "$p"; PY=$HERE/pyfigs; echo "   found $p"; return 0; }
  done
  if init_module; then
    for m in $(module -t avail 2>&1 | grep -iE '^(python|anaconda|miniconda|miniforge|conda|intelpython|intel-python|py)[^ ]*' | sort -r | head -12); do
      echo "   module $m" >>"$LOG"
      p=$( set +u; module purge >/dev/null 2>&1; module load "$m" >/dev/null 2>&1 || exit 1
           for q in python3 python; do c=$(command -v $q 2>/dev/null) || continue
             PYTHONNOUSERSITE=1 tmo 120 "$c" -c "$PYTEST" >/dev/null 2>&1 </dev/null && { echo "$c|$LD_LIBRARY_PATH"; exit 0; }; done; exit 1 ) || continue
      NUS=1 LD_LIBRARY_PATH=${p#*|} wrap "${p%%|*}"; PY=$HERE/pyfigs; echo "   found module $m: ${p%%|*}"; return 0
    done
  fi
  # 3. another Singularity image of this machine with python, numpy, matplotlib
  SING=$( { init_module && { set +u; module load singularity || module load apptainer; set -u; } ; } >/dev/null 2>&1
          command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null )
  if [ -n "$SING" ]; then
    echo "-- 3. another image (.sif) of this machine with python, numpy, matplotlib ..."
    echo "== .sif images" >>"$LOG"
    for s in $(ls /home/app/*.sif /home/app/*/*.sif /home/app/*/*/*.sif /home/app/*/*/*/*.sif "$HOME"/*.sif "$HOME"/*/*.sif 2>/dev/null | head -40); do
      [ "$s" = "$HERE/sifpy.sif" ] && continue
      echo "   $s" >>"$LOG"
      tmo 120 $SING exec -e "$s" python3 -c "$PYTEST" >>"$LOG" 2>&1 </dev/null || continue
      { echo "#!/bin/bash"; echo "# python with numpy, matplotlib in $s (written by figs.sh)"
        echo "B=\"-B $HERE\"; case \"\$(pwd -P)/\" in \"$HERE\"/*|\"\${HOME:-/nonexistent}\"/*) ;; *) B=\"\$B -B \$(pwd -P)\";; esac"
        echo "export SINGULARITYENV_MPLCONFIGDIR=\"$HERE/.mplconfig\" APPTAINERENV_MPLCONFIGDIR=\"$HERE/.mplconfig\""
        echo "exec $SING exec -e \$B \"$s\" python3 \"\$@\""; } > pyfigs
      chmod +x pyfigs; PY=$HERE/pyfigs; echo "   found $s"; return 0
    done
  fi
  # 4. the python3 of this node (it has numpy: go.sh built the gouge with it) + matplotlib from pip (binary wheels only) in pyfigs_lib/
  echo "-- 4. python3 of this node + matplotlib from pip into pyfigs_lib/ (binary wheels, at most 10 min) ..."
  echo "== pip into pyfigs_lib" >>"$LOG"
  p=$(command -v python3 2>/dev/null) || return 1
  "$p" --version >>"$LOG" 2>&1
  rm -rf pyfigs_lib; mkdir -p pyfigs_lib
  PYTHONNOUSERSITE=1 tmo 600 "$p" -m pip install --no-cache-dir --disable-pip-version-check -q --only-binary=:all: --timeout 30 --retries 2 \
      --target "$HERE/pyfigs_lib" numpy matplotlib scipy >>"$LOG" 2>&1 \
    || PYTHONNOUSERSITE=1 tmo 600 "$p" -m pip install --no-cache-dir --disable-pip-version-check -q --only-binary=:all: --timeout 30 --retries 2 \
      --target "$HERE/pyfigs_lib" numpy matplotlib >>"$LOG" 2>&1
  PYTHONNOUSERSITE=1 PYTHONPATH=$HERE/pyfigs_lib py_ok "$p" && { NUS=1 wrap "$p" "$HERE/pyfigs_lib"; PY=$HERE/pyfigs; echo "   found $p + pyfigs_lib"; return 0; }
  return 1
}
if ! find_py; then
  echo
  echo "figs.sh: no python with numpy and matplotlib here.  What was tried: $LOG  (tail -30 $LOG)"
  echo "  - singularity only on the compute nodes?   env JOB=1 ./figs.sh$RUNS"
  echo "  - no internet on this node?   a proxy (tcsh): setenv https_proxy http://<proxy>:<port>; ./figs.sh$RUNS"
  echo "  - else copy the results to a computer with python (numpy, matplotlib) and run there:  python3 analysis/hot_spots.py RUN_N22_one --outdir RUN_N22_one/figs"
  exit 1
fi
echo "python: $PY  ($($PY -c "import sys, numpy, matplotlib; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'matplotlib', matplotlib.__version__)" 2>/dev/null))"
$PY -c "import scipy" >/dev/null 2>&1 || echo "  (no scipy: everything except fig9)"
# the jobs of the run folders draw the figures themselves from now on
for r in RUN_*/; do
  [ -f "$r/config.env" ] || continue
  { grep -v "# go.sh: python" "$r/config.env" || true; echo "PYTHON=\"$PY\"   # go.sh: python"; } > "$r/config.env.tmp" && mv -f "$r/config.env.tmp" "$r/config.env"
done

# ---- the figures (analysis/ of this folder: the newest scripts after git pull)
A=analysis; [ -f $A/hot_spots.py ] || A=.
NICE=""; command -v nice >/dev/null 2>&1 && NICE="nice -n 10"
hs() {  # folder
  local r=$1; mkdir -p "$r/figs"
  echo "== hot_spots.py $r  (log: $r/figs/hot_spots.log)"
  $NICE $PY -u $A/hot_spots.py "$HERE/$r" --outdir "$HERE/$r/figs" --budget "$BUDGET" 2>&1 | tee "$r/figs/hot_spots.log" | sed 's/^/   /'
  [ "${PIPESTATUS[0]}" = 0 ] || echo "   hot_spots.py $r failed (see above; $r/figs/hot_spots.log)"
}
for r in $RUNS; do
  # are the particles stacked? (from the trajectories of 01, 02, 03, 04; quick)
  if ls $r/0*.dump $r/shear_*.dump >/dev/null 2>&1; then
    mkdir -p "$r/figs"
    echo "== stacks.py $r  (are the particles stacked?)"
    $NICE $PY -u $A/stacks.py "$HERE/$r" --outdir "$HERE/$r/figs" 2>&1 | tee "$r/figs/stacks.log" | sed 's/^/   /'
  fi
  if ls $r/w*/shear_*.dat >/dev/null 2>&1; then
    mkdir -p "$r/figs"
    echo "== lin_figs.py $r  (log: $r/figs/lin_figs.log)"
    $NICE $PY -u $A/lin_figs.py "$HERE/$r" --outdir "$HERE/$r/figs" 2>&1 | tee "$r/figs/lin_figs.log" | sed 's/^/   /'
    [ "${PIPESTATUS[0]}" = 0 ] || echo "   lin_figs.py $r failed (see above; $r/figs/lin_figs.log)"
    for w in $r/w*/; do
      w=${w%/}; ls $w/shear_*.dat >/dev/null 2>&1 || continue
      for f in gouge.data config.env; do [ -e $w/$f ] || ln -s ../$f $w/$f; done
      hs $w
    done
  fi
  ls $r/shear_*.dat >/dev/null 2>&1 && hs $r
  # a group of folders to compare (go_nro.sh: lin_group.txt): Lin & Hatano's figures of all of them, in this folder's figs/
  if [ -f $r/lin_group.txt ] && ls $r/shear_*.dat >/dev/null 2>&1; then
    LG=$(for d in $(cat $r/lin_group.txt); do ls $d/shear_*.dat >/dev/null 2>&1 && echo "$HERE/$d"; done | tr '\n' ' ')
    mkdir -p "$r/figs"
    echo "== lin_figs.py $LG (log: $r/figs/lin_figs.log)"
    $NICE $PY -u $A/lin_figs.py $LG --outdir "$HERE/$r/figs" 2>&1 | tee "$r/figs/lin_figs.log" | sed 's/^/   /'
    [ "${PIPESTATUS[0]}" = 0 ] || echo "   lin_figs.py failed (see above; $r/figs/lin_figs.log)"
  fi
done
echo
echo "figures:"
for r in $RUNS; do for d in $r/figs $r/w*/figs; do [ -d $d ] && echo "  $HERE/$d: $(ls $d/*.png 2>/dev/null | wc -l | tr -d ' ') png"; done; done
echo "(./figs.sh again later: the new results and trajectory frames are added)"
