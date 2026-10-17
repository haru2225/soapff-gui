#!/bin/sh
# double-click in Finder (macOS) to start soapff-GUI.  PY: a Python with soapff-cui's dependencies (numpy, ase, dscribe) and
# tkinter (soapff-cui itself is found next to this folder or in $SOAPFF_HOME); LAMMPS (lmp) and mpirun are looked for next to it first, then on PATH.
PY=${PY:-$HOME/ScoreMD/.pixi/envs/default/bin/python}
HERE=$(cd "$(dirname "$0")" && pwd)
PYTHONPATH="$HERE:$PYTHONPATH" exec "$PY" -m soapff_gui
