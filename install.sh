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
# --- SOVA (PDF, RDF, S(Q)) into the same environment ---
if "$BIN/python" -c "import sys; assert sys.version_info >= (3, 11)" 2>/dev/null; then
  "$BIN/python" -m pip install -q sovapy && echo "SOVA (sovapy) installed" \
    || echo "WARNING: pip install sovapy failed: PDF / RDF / S(Q) are not available (try 'pip install sovapy' by hand)"
else
  echo "WARNING: SOVA (sovapy) needs Python 3.11 or newer: PDF / RDF / S(Q) are not available (re-run with PYTHON=python3.12 sh install.sh)"
fi

# --- Debyer (XRD) built into ~/.soapff-gui/bin/debyer (SOAPFF_HOME to change) ---
SH=${SOAPFF_HOME:-$HOME/.soapff-gui}
if command -v debyer >/dev/null 2>&1 || [ -x "$SH/bin/debyer" ]; then
  echo "Debyer found: $(command -v debyer || echo "$SH/bin/debyer")"
else
  need=""
  for t in autoreconf automake gengetopt make cc git; do command -v $t >/dev/null 2>&1 || need="$need $t"; done
  if [ -n "$need" ] && [ "$(uname)" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
    echo "Debyer: installing the build tools with Homebrew (autoconf automake gengetopt libtool)"
    brew install autoconf automake gengetopt libtool >/dev/null 2>&1 || true
    need=""
    for t in autoreconf automake gengetopt make cc git; do command -v $t >/dev/null 2>&1 || need="$need $t"; done
  fi
  if [ -n "$need" ]; then
    echo "WARNING: Debyer (XRD) was not built; missing:$need"
    echo "         Ubuntu/Debian: sudo apt install git build-essential autoconf automake libtool gengetopt   (then run this again)"
    echo "         macOS: brew install autoconf automake gengetopt libtool   (then run this again)"
  else
    echo "Debyer: building (git clone, autoreconf, configure, make)"
    mkdir -p "$SH/bin" && rm -rf "$SH/debyer-src" \
      && git clone -q --depth 1 https://github.com/wojdyr/debyer "$SH/debyer-src" \
      && (cd "$SH/debyer-src" && autoreconf -i >/dev/null 2>&1 && ./configure -q --without-zlib --without-bzlib CFLAGS="-O3 -ffast-math" >/dev/null \
          && make -j2 >/dev/null 2>&1) \
      && cp "$SH/debyer-src/debyer/debyer" "$SH/bin/debyer" \
      && echo "Debyer installed: $SH/bin/debyer" \
      || echo "WARNING: building Debyer failed (XRD is not available); see https://github.com/wojdyr/debyer"
  fi
fi

echo "installed: $BIN/soapff-gui"
if [ "$(uname)" = "Darwin" ]; then
  "$BIN/soapff-gui" --install-app ${SOAPFF_APP_DIR:+"$SOAPFF_APP_DIR"}
  echo "start it from Finder / Launchpad (soapff-GUI), or with: $BIN/soapff-gui"
else
  echo "start it with: $BIN/soapff-gui"
fi
