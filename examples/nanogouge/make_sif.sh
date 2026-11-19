#!/bin/bash
# The python of this folder: python 3.11 + numpy, scipy, matplotlib in a Singularity / Apptainer image, for everything that needs them
# (build_gouge.py, the figures hot_spots.py, lin_figs.py); LAMMPS does not.  No pip on the supercomputer's own python.
#   ./make_sif.sh            (on the login node; go.sh and figs.sh call it when ./sifpy is missing)
# It pulls python:3.11-slim (no root needed) to sifpy.sif (or, when that fails, unpacks it into the folder sifpy_sandbox/), installs numpy,
# scipy, matplotlib with the image's pip into sifpy_lib/ (binary wheels only, a minute or two) and writes the command ./sifpy (= python3 of
# the image):  ./sifpy build_gouge.py ...,  ./sifpy analysis/hot_spots.py RUN_N22
# Every run folder RUN_*/ gets PYTHON=<this folder>/sifpy in its config.env: the next job of a running chain draws the figures (RUN_*/figs).
# singularity / apptainer: on the PATH, else from the module system (module load singularity / apptainer / singularityce, also from a bash
# started by tcsh), else in /usr/local, /opt, /home/app.  Its cache and temporary files go into .sifcache/ here (not the home quota, not /tmp).
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
# the module command in this bash (tcsh's 'module' alias is not passed to a bash script)
init_module() {
  type module >/dev/null 2>&1 && return 0
  local f
  for f in "${MODULESHOME:-/nonexistent}/init/bash" "${LMOD_PKG:-/nonexistent}/init/bash" /etc/profile.d/modules.sh /etc/profile.d/lmod.sh \
           /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash /opt/modules/default/init/bash; do
    [ -f "$f" ] && { set +u; . "$f" >/dev/null 2>&1; set -u; type module >/dev/null 2>&1 && return 0; }
  done
  return 1
}
find_sing() {
  SING=$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true); [ -n "$SING" ] && return 0
  if init_module; then
    local m
    for m in singularity apptainer singularityce singularity-ce $(module -t avail 2>&1 | grep -iE '^(singularity|apptainer)' | sort -r | head -3); do
      set +u; module load "$m" >/dev/null 2>&1; set -u
      SING=$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true); [ -n "$SING" ] && { echo "make_sif.sh: module load $m"; return 0; }
    done
  fi
  SING=$(ls /usr/local/bin/singularity /usr/local/bin/apptainer /opt/singularity*/bin/singularity /opt/apptainer*/bin/apptainer \
            /usr/local/singularity*/bin/singularity /home/app/singularity*/bin/singularity /home/app/singularity*/*/bin/singularity \
            /home/app/apptainer*/bin/apptainer /home/app/apptainer*/*/bin/apptainer 2>/dev/null | head -1)
  [ -n "$SING" ]
}
find_sing || { echo "make_sif.sh: no singularity / apptainer on this node (PATH, module avail, /usr/local, /opt, /home/app).  'module avail' shows the modules;"
               echo "  on some machines it is only on the compute nodes: env JOB=1 ./figs.sh RUN_N22_one (a short job)"; exit 1; }
SV=$($SING --version 2>/dev/null | head -1)
echo "make_sif.sh: using $SING ($SV)"
# its cache and temporary files here
mkdir -p .sifcache/tmp
export SINGULARITY_CACHEDIR="$HERE/.sifcache" APPTAINER_CACHEDIR="$HERE/.sifcache" SINGULARITY_TMPDIR="$HERE/.sifcache/tmp" APPTAINER_TMPDIR="$HERE/.sifcache/tmp"
export SINGULARITY_PULLFOLDER="$HERE" TMPDIR="${TMPDIR:-$HERE/.sifcache/tmp}"
IMG=${SIF_IMAGE:-docker://python:3.11-slim}
IMGF=""
[ -s sifpy.sif ] && IMGF="$HERE/sifpy.sif"
[ -z "$IMGF" ] && [ -x sifpy_sandbox/usr/local/bin/python3 ] && IMGF="$HERE/sifpy_sandbox"
if [ -z "$IMGF" ]; then
  echo "pulling $IMG to sifpy.sif (at most 10 min) ..."
  case "$SV" in
    *" 2."*|2.*) tmo 600 $SING pull --name sifpy.sif "$IMG" ;;     # Singularity 2.x
    *)           tmo 600 $SING pull sifpy.sif "$IMG" ;;
  esac
  if [ -s sifpy.sif ]; then IMGF="$HERE/sifpy.sif"
  else
    rm -f sifpy.sif
    echo "the pull did not make sifpy.sif; trying an unpacked image (a folder, no squashfs needed): sifpy_sandbox/ ..."
    rm -rf sifpy_sandbox
    tmo 600 $SING build --sandbox sifpy_sandbox "$IMG" && [ -x sifpy_sandbox/usr/local/bin/python3 ] && IMGF="$HERE/sifpy_sandbox"
  fi
  if [ -z "$IMGF" ]; then
    rm -rf sifpy_sandbox
    echo "could not get $IMG (no internet on this node? a proxy: setenv https_proxy http://<proxy>:<port> in tcsh, then ./make_sif.sh again)."
    echo "Else run make_sif.sh on a computer with internet and Singularity/Apptainer, copy sifpy.sif and sifpy_lib/ here, run ./make_sif.sh again."
    exit 1
  fi
fi
echo "image: $IMGF"
# the packages, with the image's pip, into sifpy_lib/ (a folder here: no root, no writable image needed)
ok() { SINGULARITYENV_PYTHONPATH="$HERE/sifpy_lib" APPTAINERENV_PYTHONPATH="$HERE/sifpy_lib" SINGULARITYENV_PYTHONNOUSERSITE=1 APPTAINERENV_PYTHONNOUSERSITE=1 \
       tmo 120 $SING exec -e -B "$HERE" "$IMGF" python3 -c "import numpy, scipy, matplotlib; print('numpy', numpy.__version__, 'scipy', scipy.__version__, 'matplotlib', matplotlib.__version__)"; }
if ! ok >/dev/null 2>&1; then
  echo "installing numpy, scipy, matplotlib into sifpy_lib/ with the pip of the image (binary wheels only, at most 10 min) ..."
  mkdir -p sifpy_lib
  PYTHONNOUSERSITE=1 tmo 600 $SING exec -B "$HERE" "$IMGF" python3 -m pip install --no-cache-dir --disable-pip-version-check -q --only-binary=:all: \
      --timeout 30 --retries 2 --target "$HERE/sifpy_lib" numpy scipy matplotlib \
    || { echo "pip in the image failed (no internet to pypi.org? a proxy: setenv https_proxy http://<proxy>:<port>).  Without internet: see the head of make_sif.sh."; exit 1; }
fi
ok || { echo "the image does not import numpy, scipy, matplotlib (see above)"; exit 1; }
# the command: ./sifpy script.py ... (this folder is bound, and the current folder when it is not below this one or the home folder)
cat > sifpy <<EOF
#!/bin/bash
# python 3.11 with numpy, scipy, matplotlib in $IMGF (made by make_sif.sh)
B="-B $HERE"; case "\$(pwd -P)/" in "$HERE"/*|"\${HOME:-/nonexistent}"/*) ;; *) B="\$B -B \$(pwd -P)";; esac
export SINGULARITYENV_PYTHONPATH="$HERE/sifpy_lib" APPTAINERENV_PYTHONPATH="$HERE/sifpy_lib" SINGULARITYENV_PYTHONNOUSERSITE=1 APPTAINERENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_MPLCONFIGDIR="$HERE/sifpy_lib/.mpl" APPTAINERENV_MPLCONFIGDIR="$HERE/sifpy_lib/.mpl"
exec $SING exec -e \$B "$IMGF" python3 "\$@"
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
echo "done: $HERE/sifpy  (figures now: ./figs.sh RUN_N22_one)"
