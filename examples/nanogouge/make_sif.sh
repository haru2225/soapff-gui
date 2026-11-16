#!/bin/bash
# A python for the figures (hot_spots.py: numpy, scipy, matplotlib) in a Singularity / Apptainer image; the runs do not need it.
#   ./make_sif.sh            (on the login node; go.sh calls it when no python with matplotlib is found)
# It pulls python:3.11-slim (no root needed) to figs.sif, installs numpy, scipy, matplotlib with the image's pip into figs_pylib/ (binary wheels:
# a minute or two; the old pip of a system python 3.6 compiles matplotlib from source, which looks like a hang) and writes the command figs_python.
# Every run folder RUN_*/ gets PYTHON=<this folder>/figs_python in its config.env: the next job of a running chain draws the figures (RUN_*/figs).
# Without internet on this node: build figs.sif and figs_pylib/ on a computer with internet and Singularity/Apptainer with this script, copy both
# here (scp -r figs.sif figs_pylib <user>@<host>:<this folder>/), then run ./make_sif.sh again (it uses them).
set -u
cd "$(dirname "$0")"
HERE=$(pwd -P)
tmo() { local t=$1; shift; if command -v timeout >/dev/null 2>&1; then timeout "$t" "$@"; else "$@"; fi; }
SING=$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true)
if [ -z "$SING" ]; then
  for m in singularity apptainer; do command -v module >/dev/null 2>&1 && module load $m >/dev/null 2>&1; done
  SING=$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true)
fi
[ -n "$SING" ] || { echo "no singularity / apptainer on this node (also not after 'module load singularity'): no image"; exit 1; }
echo "using $SING ($($SING --version 2>/dev/null | head -1))"
IMG=${FIGS_IMAGE:-docker://python:3.11-slim}
if [ ! -s figs.sif ]; then
  echo "pulling $IMG to figs.sif (at most 15 min) ..."
  if ! tmo 900 $SING pull figs.sif "$IMG"; then
    rm -f figs.sif
    echo "could not pull $IMG (no internet on this node? a proxy: setenv https_proxy http://<proxy>:<port> in tcsh, then again)."
    echo "Else build figs.sif and figs_pylib/ with this script on a computer with internet and Singularity/Apptainer, copy them here, run ./make_sif.sh again."
    exit 1
  fi
fi
# the packages, with the image's pip, into figs_pylib/ (a folder here: no root, no writable image needed)
ok() { SINGULARITYENV_PYTHONPATH="$HERE/figs_pylib" APPTAINERENV_PYTHONPATH="$HERE/figs_pylib" SINGULARITYENV_PYTHONNOUSERSITE=1 APPTAINERENV_PYTHONNOUSERSITE=1 \
       tmo 120 $SING exec -e -B "$HERE" figs.sif python3 -c "import numpy, scipy, matplotlib; print('numpy', numpy.__version__, 'scipy', scipy.__version__, 'matplotlib', matplotlib.__version__)"; }
if ! ok >/dev/null 2>&1; then
  echo "installing numpy, scipy, matplotlib into figs_pylib/ with the pip of the image (at most 20 min) ..."
  mkdir -p figs_pylib
  PYTHONNOUSERSITE=1 tmo 1200 $SING exec -B "$HERE" figs.sif python3 -m pip install --no-cache-dir --disable-pip-version-check -q --target "$HERE/figs_pylib" numpy scipy matplotlib \
    || { echo "pip in the image failed (no internet? see above).  Without internet: see the head of this script."; exit 1; }
fi
ok || { echo "the image does not import numpy, scipy, matplotlib (see above)"; exit 1; }
# the command the jobs call: figs_python hot_spots.py ... (the run folder is bound when it is not below the home folder)
cat > figs_python <<EOF
#!/bin/bash
# python 3.11 with numpy, scipy, matplotlib in $HERE/figs.sif (make_sif.sh)
B="-B $HERE"; case "\$(pwd -P)/" in "$HERE"/*|"\${HOME:-/nonexistent}"/*) ;; *) B="\$B -B \$(pwd -P)";; esac
export SINGULARITYENV_PYTHONPATH="$HERE/figs_pylib" APPTAINERENV_PYTHONPATH="$HERE/figs_pylib" SINGULARITYENV_PYTHONNOUSERSITE=1 APPTAINERENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_MPLCONFIGDIR="$HERE/figs_pylib/.mpl" APPTAINERENV_MPLCONFIGDIR="$HERE/figs_pylib/.mpl"
exec $SING exec -e \$B "$HERE/figs.sif" python3 "\$@"
EOF
chmod +x figs_python
./figs_python -c "import matplotlib; print('figs_python works: matplotlib', matplotlib.__version__)" || { echo "figs_python does not start (see above)"; exit 1; }
# every run folder here: its jobs use it from their next job on
for r in RUN_*/; do
  [ -f "$r/config.env" ] || continue
  { grep -v "# go.sh: python for the figures" "$r/config.env" || true; echo "PYTHON=\"$HERE/figs_python\"   # go.sh: python for the figures"; } > "$r/config.env.tmp" && mv -f "$r/config.env.tmp" "$r/config.env"
  rm -f "$r/.figpy_none"; echo "  ${r%/}: its next jobs draw the figures in ${r%/}/figs"
done
HS=analysis/hot_spots.py; [ -f hot_spots.py ] && HS=hot_spots.py
echo "done: $HERE/figs_python (figures by hand: ./figs_python $HS RUN_N67 --outdir RUN_N67/figs)"
