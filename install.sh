#!/bin/sh
# Install soapff-GUI from this checkout and (macOS) make ~/Applications/soapff-GUI.app.
#   sh install.sh                    a new environment in ~/.soapff-gui/venv (Python: $PYTHON, default python3)
#   PYTHON=/path/to/python sh install.sh --here   into that Python's own environment (e.g. one that already has LAMMPS)
# Needs Python >= 3.9 with tkinter.  LAMMPS (lmp, mpirun) is not installed here.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
"$PY" -c "import sys, tkinter; assert sys.version_info >= (3, 9)" 2>/dev/null || {
  echo "$PY: Python 3.9 or newer with tkinter is needed (set PYTHON=/path/to/python)"; exit 1; }
if [ "$1" = "--here" ]; then
  "$PY" -m pip install -q -e "$HERE"
  BIN=$(dirname "$("$PY" -c 'import sys; print(sys.executable)')")
else
  VENV=${SOAPFF_VENV:-$HOME/.soapff-gui/venv}
  "$PY" -m venv "$VENV"
  "$VENV/bin/python" -m pip install -q --upgrade pip
  "$VENV/bin/python" -m pip install -q -e "$HERE"
  BIN=$VENV/bin
fi
echo "installed: $BIN/soapff-gui"
if [ "$(uname)" = "Darwin" ]; then
  "$BIN/soapff-gui" --install-app ${SOAPFF_APP_DIR:+"$SOAPFF_APP_DIR"}
  echo "start it from Finder / Launchpad (soapff-GUI), or with: $BIN/soapff-gui"
else
  echo "start it with: $BIN/soapff-gui"
fi
