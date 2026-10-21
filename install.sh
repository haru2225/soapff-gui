#!/bin/sh
# Install soapff-GUI from this checkout and (macOS) make ~/Applications/soapff-GUI.app (with its icon).
#   git clone git@github.com:haru2225/soapff-gui.git && sh soapff-gui/install.sh
#   sh install.sh                    a new environment in ~/.soapff-gui/venv
#   PYTHON=/path/to/python sh install.sh --here   into that Python's own environment (e.g. one that already has LAMMPS)
# Python: $PYTHON, else the first of python3, python3.13 ... python3.9 (PATH, Homebrew, python.org) that has tkinter with
# Tk 8.6 or newer; when there is none, uv (if installed) makes the environment with a Python of its own (with tkinter).
# LAMMPS + Open MPI and OVITO (conda-forge / conda.ovito.org, no sudo) and SOVA (pip) are installed too: see the steps below; SOAPFF_NO_LAMMPS=1 / SOAPFF_NO_OVITO=1 skip LAMMPS / OVITO.
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
# --- SOVA (PDF, RDF, S(Q)) into the same environment ---
if "$BIN/python" -c "import sys; assert sys.version_info >= (3, 11)" 2>/dev/null; then
  "$BIN/python" -m pip install -q sovapy && echo "SOVA (sovapy) installed" \
    || echo "WARNING: pip install sovapy failed: PDF / RDF / S(Q) are not available (try 'pip install sovapy' by hand)"
else
  echo "WARNING: SOVA (sovapy) needs Python 3.11 or newer: PDF / RDF / S(Q) are not available (re-run with PYTHON=python3.12 sh install.sh)"
fi

# --- micromamba (for LAMMPS and OVITO below): an existing one, else fetched into ~/.soapff-gui/bin (no sudo) ---
SH=${SOAPFF_HOME:-$HOME/.soapff-gui}
get_mm() {
  MM=$(command -v micromamba || command -v mamba || true)
  [ -n "$MM" ] && return 0
  [ -x "$SH/bin/micromamba" ] && MM="$SH/bin/micromamba" && return 0
  case "$(uname)-$(uname -m)" in
    Darwin-arm64) P=osx-arm64;; Darwin-x86_64) P=osx-64;; Linux-x86_64) P=linux-64;; Linux-aarch64) P=linux-aarch64;; *) return 1;;
  esac
  echo "fetching micromamba ($P)"
  mkdir -p "$SH/bin" && curl -fsSL "https://micro.mamba.pm/api/micromamba/$P/latest" | tar -xj -C "$SH" bin/micromamba 2>/dev/null \
    && MM="$SH/bin/micromamba" && return 0
  MM=""; return 1
}

# --- LAMMPS and Open MPI (conda-forge, MPI build with KSPACE and MOLECULE) into ~/.soapff-gui/lammps ---
if command -v lmp >/dev/null 2>&1 || command -v lmp_mpi >/dev/null 2>&1 || [ -x "$SH/lammps/bin/lmp" ]; then
  echo "LAMMPS found: $(command -v lmp || command -v lmp_mpi || echo "$SH/lammps/bin/lmp")"
elif [ "${SOAPFF_NO_LAMMPS:-}" = "1" ]; then
  echo "LAMMPS: skipped (SOAPFF_NO_LAMMPS=1)"
elif get_mm; then
  echo "LAMMPS: installing lammps + openmpi from conda-forge into $SH/lammps (about 1 GB)"
  "$MM" create -y -q -p "$SH/lammps" -c conda-forge "lammps=*=*openmpi*" openmpi >/dev/null 2>&1 \
    && echo "LAMMPS installed: $SH/lammps/bin/lmp" \
    || echo "WARNING: the LAMMPS install failed; install it by hand (conda install -c conda-forge lammps openmpi, brew install lammps open-mpi, apt install lammps)"
else
  echo "WARNING: LAMMPS was not installed (micromamba could not be fetched); install it by hand:"
  echo "         conda install -c conda-forge lammps openmpi   |   brew install lammps open-mpi   |   sudo apt install lammps openmpi-bin"
fi

# --- OVITO (the application, for 'OVITO で開く'): macOS Homebrew cask; elsewhere conda.ovito.org into ~/.soapff-gui/ovito ---
if [ -d /Applications/Ovito.app ] || [ -d "$HOME/Applications/Ovito.app" ] || command -v ovito >/dev/null 2>&1 || [ -x "$SH/ovito/bin/ovito" ]; then
  echo "OVITO found"
elif [ "${SOAPFF_NO_OVITO:-}" = "1" ]; then
  echo "OVITO: skipped (SOAPFF_NO_OVITO=1)"
elif [ "$(uname)" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
  echo "OVITO: brew install --cask ovito"
  brew install --cask ovito >/dev/null 2>&1 && echo "OVITO installed" || echo "WARNING: brew install --cask ovito failed; get OVITO from https://www.ovito.org"
elif get_mm; then
  echo "OVITO: installing from conda.ovito.org into $SH/ovito (about 1.5 GB)"
  "$MM" create -y -q -p "$SH/ovito" -c https://conda.ovito.org -c conda-forge ovito >/dev/null 2>&1 \
    && echo "OVITO installed: $SH/ovito/bin/ovito" || echo "WARNING: the OVITO install failed; get OVITO from https://www.ovito.org"
else
  echo "WARNING: OVITO was not installed; get it from https://www.ovito.org (the Python module for the supercell is installed already)"
fi

echo "installed: $BIN/soapff-gui"
if [ "$(uname)" = "Darwin" ]; then
  "$BIN/soapff-gui" --install-app ${SOAPFF_APP_DIR:+"$SOAPFF_APP_DIR"}
  echo "start it from Finder / Launchpad (soapff-GUI), or with: $BIN/soapff-gui"
else
  echo "start it with: $BIN/soapff-gui"
fi
