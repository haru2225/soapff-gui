#!/bin/bash
# The python of this folder: python 3.11 + numpy, scipy, matplotlib in a Singularity / Apptainer image, for everything that needs them
# (build_gouge.py, the figures hot_spots.py, the analysis scripts); LAMMPS does not.  No pip on the supercomputer's own python.
#   ./make_sif.sh            (on the login node; go.sh calls it once when ./sifpy is missing)
# It pulls python:3.11-slim (no root needed) to sifpy.sif, installs numpy, scipy, matplotlib with the image's pip into sifpy_lib/ (binary wheels,
# a minute or two) and writes the command ./sifpy (= python3 of the image):  ./sifpy build_gouge.py ...,  ./sifpy analysis/hot_spots.py RUN_N22
# Every run folder RUN_*/ gets PYTHON=<this folder>/sifpy in its config.env: the next job of a running chain draws the figures (RUN_*/figs).
# Without internet on this node: run this script on a computer with internet and Singularity/Apptainer, copy sifpy.sif and sifpy_lib/ here
# (scp -r sifpy.sif sifpy_lib <user>@<host>:<this folder>/), then run ./make_sif.sh again (it uses them).
set -u
cd "$(dirname "$0")"
HERE=$(pwd -P)
tmo() { local t=$1; shift; if command -v timeout >/dev/null 2>&1; then timeout "$t" "$@"; else "$@"; fi; }
# made by an earlier version under other names
[ -s figs.sif ] && [ ! -s sifpy.sif ] && mv figs.sif sifpy.sif
[ -d figs_pylib ] && [ ! -d sifpy_lib ] && mv figs_pylib sifpy_lib
rm -f figs_python
SING=$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true)
if [ -z "$SING" ]; then
  for m in singularity apptainer; do command -v module >/dev/null 2>&1 && module load $m >/dev/null 2>&1; done
  SING=$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true)
fi
[ -n "$SING" ] || { echo "make_sif.sh: no singularity / apptainer on this node (also not after 'module load singularity')"; exit 1; }
echo "make_sif.sh: using $SING ($($SING --version 2>/dev/null | head -1))"
IMG=${SIF_IMAGE:-docker://python:3.11-slim}
if [ ! -s sifpy.sif ]; then
  echo "pulling $IMG to sifpy.sif (at most 10 min) ..."
  if ! tmo 600 $SING pull sifpy.sif "$IMG"; then
    rm -f sifpy.sif
    echo "could not pull $IMG (no internet on this node? a proxy: setenv https_proxy http://<proxy>:<port> in tcsh, then ./make_sif.sh again)."
    echo "Else run make_sif.sh on a computer with internet and Singularity/Apptainer, copy sifpy.sif and sifpy_lib/ here, run ./make_sif.sh again."
    exit 1
  fi
fi
# the packages, with the image's pip, into sifpy_lib/ (a folder here: no root, no writable image needed)
ok() { SINGULARITYENV_PYTHONPATH="$HERE/sifpy_lib" APPTAINERENV_PYTHONPATH="$HERE/sifpy_lib" SINGULARITYENV_PYTHONNOUSERSITE=1 APPTAINERENV_PYTHONNOUSERSITE=1 \
       tmo 120 $SING exec -e -B "$HERE" sifpy.sif python3 -c "import numpy, scipy, matplotlib; print('numpy', numpy.__version__, 'scipy', scipy.__version__, 'matplotlib', matplotlib.__version__)"; }
if ! ok >/dev/null 2>&1; then
  echo "installing numpy, scipy, matplotlib into sifpy_lib/ with the pip of the image (at most 10 min) ..."
  mkdir -p sifpy_lib
  PYTHONNOUSERSITE=1 tmo 600 $SING exec -B "$HERE" sifpy.sif python3 -m pip install --no-cache-dir --disable-pip-version-check -q --target "$HERE/sifpy_lib" numpy scipy matplotlib \
    || { echo "pip in the image failed (no internet? see above).  Without internet: see the head of make_sif.sh."; exit 1; }
fi
ok || { echo "the image does not import numpy, scipy, matplotlib (see above)"; exit 1; }
# the command: ./sifpy script.py ... (this folder is bound, and the current folder when it is not below this one or the home folder)
cat > sifpy <<EOF
#!/bin/bash
# python 3.11 with numpy, scipy, matplotlib in $HERE/sifpy.sif (made by make_sif.sh)
B="-B $HERE"; case "\$(pwd -P)/" in "$HERE"/*|"\${HOME:-/nonexistent}"/*) ;; *) B="\$B -B \$(pwd -P)";; esac
export SINGULARITYENV_PYTHONPATH="$HERE/sifpy_lib" APPTAINERENV_PYTHONPATH="$HERE/sifpy_lib" SINGULARITYENV_PYTHONNOUSERSITE=1 APPTAINERENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_MPLCONFIGDIR="$HERE/sifpy_lib/.mpl" APPTAINERENV_MPLCONFIGDIR="$HERE/sifpy_lib/.mpl"
exec $SING exec -e \$B "$HERE/sifpy.sif" python3 "\$@"
EOF
chmod +x sifpy
./sifpy -c "import numpy, scipy, matplotlib; print('./sifpy works: numpy', numpy.__version__, 'scipy', scipy.__version__, 'matplotlib', matplotlib.__version__)" \
  || { echo "./sifpy does not start (see above)"; exit 1; }
rm -f .sifpy_failed
# every run folder here: its jobs use it from their next job on
for r in RUN_*/; do
  [ -f "$r/config.env" ] || continue
  { grep -v "# go.sh: python" "$r/config.env" || true; echo "PYTHON=\"$HERE/sifpy\"   # go.sh: python"; } > "$r/config.env.tmp" && mv -f "$r/config.env.tmp" "$r/config.env"
  rm -f "$r/.figpy_none"; echo "  ${r%/}: its next jobs draw the figures in ${r%/}/figs"
done
HS=analysis/hot_spots.py; [ -f hot_spots.py ] && HS=hot_spots.py
echo "done: $HERE/sifpy  (by hand: ./sifpy $HS RUN_N22 --outdir RUN_N22/figs)"
