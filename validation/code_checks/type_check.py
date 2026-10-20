"""Check of the SOAP-type / bond-graph-type cross-check of soapff_cui.core.make (report.json "type_check"):
no false alarm on structures whose types are right, and detection of the obts -> ob error of muscovite typed with the
5-structure library without micas.  Paths are those of the author's machine.
    python type_check.py"""
import contextlib, glob, io, os, sys
sys.path.insert(0, "/Users/harutokono/clay/soapff-gui/src")
from soapff_cui.core import make

OUT = "/tmp/soapff_type_check"
EX = "/Users/harutokono/clay/soapff-cui/examples"
NEW = "/Users/harutokono/clay/soapff-gui/src/soapff_gui/data/library_clayff_bulk.npz"          # 13 structures
OLD = f"{EX}/minerals/library_clayff_bulk.npz"                                                  # 5 structures
MUSC = f"{EX}/minerals/muscovite_mech/muscovite_2M1_221.extxyz"
cases = [("montmorillonite (JAEA)", "/Users/harutokono/clay/soapff-gui/examples/montmorillonite_2W/initial.xyz", NEW),
         ("muscovite, 13-structure library", MUSC, NEW), ("muscovite, 5-structure library", MUSC, OLD)]
cases += [(os.path.basename(f)[:-7], f, NEW) for f in sorted(glob.glob(f"{EX}/minerals/eval_sheet/tests/*.extxyz"))]
cases += [(os.path.basename(f)[:-7], f, NEW) for f in sorted(glob.glob(f"{EX}/minerals/eval/tests/*_0.extxyz"))]
cases += [(os.path.basename(f)[:-7], f, NEW) for f in sorted(glob.glob(f"{EX}/minerals/eval/tests/*_001_[23].extxyz"))]
for name, st, lib in cases:
    with contextlib.redirect_stdout(io.StringIO()):
        rep = make(st, f"{OUT}/{name.replace(' ', '_')}", lib, md=False, check=False, log=lambda m: None)
    print(f"{name:34s} fallback {rep['fallback']:4d}  type_check {rep['type_check']}")
# result (2026-10-20): 0 differences in the 18 right structures; muscovite with the 5-structure library: obts->ob x 64;
# phlogopite MD frame with the library without phlogopite (eval_sheet/loo): obts->ob x 16
