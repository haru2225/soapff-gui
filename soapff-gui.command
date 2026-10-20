#!/bin/sh
# double-click in Finder (macOS) to start soapff-GUI.  PY: a Python with numpy, scipy, ase, dscribe, matplotlib and tkinter
# (soapff-cui is bundled in this folder as soapff_cui); LAMMPS (lmp) and mpirun are looked for next to it first, then on PATH.
PY=${PY:-$HOME/ScoreMD/.pixi/envs/default/bin/python}
HERE=$(cd "$(dirname "$0")" && pwd)
PYTHONPATH="$HERE:$PYTHONPATH" exec "$PY" -m soapff_gui
