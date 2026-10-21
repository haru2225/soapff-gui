"""Swelling curve from gui_run/swelling.dat (written by the soapff swelling block of soapff-gui): for each water content,
d001 = lz averaged over the last 20 ps of its NPT (and the standard deviation), against the water content in g H2O / g clay
(clay = every atom that is not water, Na included).  -> gui_run/swelling.csv, swelling.png
    python analyze_swelling.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../src"))
from soapff_cui.analysis import swelling  # noqa: E402

out = os.path.join(HERE, "gui_run")
r = swelling(out, out)
for s in r["states"]:
    print(f"{s['waters']:3d} waters  {s['water_g_per_g']:.3f} g/g  d001 {s['d001']:.2f} +- {s['d001_sd']:.2f} A")
