#!/bin/sh
# Install soapff-GUI from this checkout and (macOS) make ~/Applications/soapff-GUI.app (with its icon).
#   git clone git@github.com:haru2225/soapff-gui.git && sh soapff-gui/install.sh
#   sh install.sh                    a new environment in ~/.soapff-gui/venv
#   PYTHON=/path/to/python sh install.sh --here   into that Python's own environment (e.g. one that already has LAMMPS)
# Python: $PYTHON, else the first of python3, python3.13 ... python3.9 (PATH, Homebrew, python.org) that has tkinter with
# Tk 8.6 or newer; when there is none, uv (if installed) makes the environment with a Python of its own (with tkinter).
# LAMMPS (lmp, mpirun) is not installed here.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ok() { "$1" -c "import sys, tkinter; assert sys.version_info >= (3, 9) and tkinter.TkVersion >= 8.6" 2>/dev/null; }
PY=""
if [ -n "$PYTHON" ]; then
  ok "$PYTHON" || { echo "$PYTHON: Python 3.9 or newer with tkinter (Tk 8.6 or newer) is needed"; exit 1; }
  PY=$PYTHON
else
  for c in python3 python3.13 python3.12 python3.11 python3.10 python3.9 \
           /opt/homebrew/bin/python3 /usr/local/bin/python3 /Library/Frameworks/Python.framework/Versions/Current/bin/python3; do
    p=$(command -v "$c" 2>/dev/null) || continue
    if ok "$p"; then PY=$p; break; fi
  done
fi
NOPY="No Python with tkinter (Tk 8.6 or newer) was found. Install one, e.g. 'brew install python-tk@3.12' or uv
('brew install uv'), or give it with PYTHON=/path/to/python, and run this again."
if [ "$1" = "--here" ]; then
  [ -n "$PY" ] || { echo "$NOPY"; exit 1; }
  "$PY" -m pip install -q -e "$HERE"
  BIN=$(dirname "$("$PY" -c 'import sys; print(sys.executable)')")
else
  VENV=${SOAPFF_VENV:-$HOME/.soapff-gui/venv}
  if [ -n "$PY" ]; then
    echo "Python: $PY"
    "$PY" -m venv "$VENV"
  elif command -v uv >/dev/null 2>&1; then
    echo "No Python with tkinter here: uv makes the environment with a Python 3.12 of its own"
    uv venv -q --seed --clear --python-preference only-managed --python 3.12 "$VENV"
  else
    echo "$NOPY"; exit 1
  fi
  ok "$VENV/bin/python" || { echo "$VENV/bin/python has no tkinter. $NOPY"; exit 1; }
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
