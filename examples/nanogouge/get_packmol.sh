#!/bin/bash
# packmol for rewater.py (the water put into a fabric: go.sh FROMFAB=...), made once in ./bin/packmol from its source
# (github.com/m3g/packmol) with the Fortran compiler of this computer: gfortran, ifx or ifort (also tried from the module system).
#   ./get_packmol.sh                   on the login node, about a minute; needs internet once (git clone or curl)
#   env FC=/path/to/ifort ./get_packmol.sh     another compiler
# Without it rewater.py puts the water on a grid instead (numpy only).  A packmol of your own: PACKMOL=/path/packmol, or packmol on the PATH.
set -u
cd "$(dirname "$0")"
[ -x bin/packmol ] && { echo "bin/packmol is there already (rm bin/packmol to make it again)"; exit 0; }
init_module() {
  type module >/dev/null 2>&1 && return 0
  local f
  for f in "${MODULESHOME:-/nonexistent}/init/bash" "${LMOD_PKG:-/nonexistent}/init/bash" /etc/profile.d/modules.sh /etc/profile.d/lmod.sh \
           /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash /opt/modules/default/init/bash; do
    [ -f "$f" ] && { set +u; . "$f" >/dev/null 2>&1; set -u; type module >/dev/null 2>&1 && return 0; }
  done
  return 1
}
find_fc() { local c; for c in gfortran ifx ifort; do command -v $c >/dev/null 2>&1 && { FC=$(command -v $c); return 0; }; done; return 1; }
FC=${FC:-}
if [ -z "$FC" ] && ! find_fc; then
  if init_module; then
    for m in gcc gnu gfortran intel intel-oneapi oneapi intel-compiler compiler; do
      module load $m >/dev/null 2>&1 && find_fc && { echo "(module load $m)"; break; }
    done
  fi
fi
[ -n "$FC" ] || { echo "no Fortran compiler found (gfortran, ifx, ifort; also not from 'module load'): env FC=/path/to/compiler ./get_packmol.sh"; exit 1; }
echo "Fortran compiler: $FC"
mkdir -p bin .packmol_src
cd .packmol_src
if [ ! -f packmol/Makefile ]; then
  rm -rf packmol
  git clone -q --depth 1 https://github.com/m3g/packmol.git packmol 2>/dev/null \
  || { curl -fsSL -o packmol.tar.gz https://github.com/m3g/packmol/archive/refs/heads/master.tar.gz && mkdir -p packmol \
       && tar xzf packmol.tar.gz -C packmol --strip-components 1; } \
  || { echo "could not get the packmol source (no internet on this node?).  On a computer with internet: git clone https://github.com/m3g/packmol.git,"
       echo "copy the folder to $(pwd)/packmol and run ./get_packmol.sh again"; exit 1; }
fi
cd packmol
make clean >/dev/null 2>&1
if make FORTRAN="$FC" FLAGS="-O2" > ../make.log 2>&1 && [ -x packmol ]; then
  cp packmol ../../bin/packmol && echo "made bin/packmol: rewater.py uses it ($(./packmol < /dev/null 2>&1 | grep -m1 Version | tr -s ' '))"
else
  echo "make failed: $(pwd)/../make.log"; tail -5 ../make.log; exit 1
fi
