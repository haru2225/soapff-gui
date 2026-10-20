"""Leave-one-out: for each sheet silicate, a library without its reference (refs_clayff/<name>_bulk.data); its last MD
frame (md_sheet.py) typed with that library.  -> results_loo.json"""
import sys, os, io, glob, json, contextlib, collections, csv
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui"); sys.path.insert(0, "/Users/harutokono/clay/mincharge")
sys.path.insert(0, "/Users/harutokono/clay/soapff-gui/src")
from ase import Atoms
from ase.io import read
from mincharge.ff import ff_types
from soapff_cui.core import make, build_library
from soapff_cui.clayff import Q
from soapff_cui.sova_structure import frames
W = os.path.dirname(os.path.abspath(__file__)); R = "/Users/harutokono/clay/soapff-cui/examples/minerals/refs_clayff"
names = sys.argv[1:] or ["muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "kaolinite", "lizerdite", "talc", "pyrophyllite_2"]
out = json.load(open(f"{W}/results_loo.json")) if os.path.isfile(f"{W}/results_loo.json") else {}
for name in names:
    if not os.path.isfile(f"{W}/md/{name}/final.data"):
        continue
    refs = [f for f in sorted(glob.glob(f"{R}/*.data")) if os.path.basename(f) != f"{name}_bulk.data"]
    lib = f"{W}/loo/lib_without_{name}.npz"; os.makedirs(f"{W}/loo", exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        build_library(refs, (), lib)
    test = read(f"{W}/tests/{name}.extxyz"); truth = ff_types(test)
    st, cell, pos = list(frames(f"{W}/md/{name}/traj.lammpstrj"))[-1]
    last = Atoms(test.get_chemical_symbols(), positions=pos, cell=cell, pbc=True)
    d = f"{W}/loo/{name}"
    err = None
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            make(last, d, lib, nvt_ps=1, check=True, log=lambda m: None)
    except ValueError as e:
        err = str(e).splitlines()[0]
        with contextlib.redirect_stdout(io.StringIO()):
            make(last, d, lib, nvt_ps=1, check=False, log=lambda m: None)
    rows = list(csv.DictReader(open(f"{d}/charges.csv")))
    bad = collections.Counter(f"{truth[i]}->{r['clayff_type']}" for i, r in enumerate(rows)
                              if r["clayff_type"] != truth[i] or abs(float(r["charge"]) - Q[truth[i]]) > 0.01)
    out[name] = dict(n=len(rows), wrong=sum(bad.values()), confusions=dict(bad), stopped=err,
                     fallback=sum(r["source"] == "fallback" for r in rows), minsim=min(float(r["similarity"]) for r in rows))
    print(name, json.dumps(out[name]), flush=True)
json.dump(out, open(f"{W}/results_loo.json", "w"), indent=1)
