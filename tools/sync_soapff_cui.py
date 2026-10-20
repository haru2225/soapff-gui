#!/usr/bin/env python3
"""Copy the soapff-cui package of a soapff-cui checkout into this repository (soapff_cui/, bundled so that installing
soapff-gui alone is enough).  The writer/typing helpers come from soapff_cui/_clayff (a copy of mincharge), so imports of
mincharge are rewritten.  Also copies the library.
    python tools/sync_soapff_cui.py ../soapff-cui"""
import os, shutil, sys
SRC = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..", "soapff-cui"))
DST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "soapff_cui")
for f in ("core.py", "clayff.py", "analysis.py", "sova_structure.py", "__main__.py"):
    t = open(os.path.join(SRC, "soapff_cui", f)).read()
    t = t.replace("from mincharge.chem import", "from ._clayff.chem import").replace("from mincharge.ff import", "from ._clayff.ff import")
    if "import mincharge" in t or "from mincharge" in t:
        sys.exit(f"{f}: unexpected mincharge import, update this script")
    open(os.path.join(DST, f), "w").write(t)
mc = os.path.join(os.path.dirname(SRC), "mincharge", "mincharge")
if os.path.isdir(mc):
    for f in ("bonds.py", "chem.py", "ff.py"):
        shutil.copy(os.path.join(mc, f), os.path.join(DST, "_clayff", f))
lib = os.path.join(SRC, "examples", "minerals", "library_clayff_bulk.npz")
if os.path.isfile(lib):
    shutil.copy(lib, os.path.join(os.path.dirname(DST), "soapff_gui", "data", "library_clayff_bulk.npz"))
print("copied soapff-cui from", SRC)
