"""(eval_sheet copy: the earlier tests of ../eval with library_clayff_sheet.npz; part E not run)
Evaluate library_clayff_bulk.npz (montmorillonite 2W initial structure + bulk quartz, albite, brucite, pyrophyllite):
 A one-to-one check (done by check-library), B JAEA initial structure -> labels, C bulk supercells, D interfaces with
 random water (20 per mineral), E MD (minimise, 10 ps relax NVT, 10 ps NVT, flexible SPC, 0.5 fs)."""
import sys, os, glob, io, contextlib, csv, json, collections, subprocess, numpy as np
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui"); sys.path.insert(0, "/Users/harutokono/clay/mincharge")
from ase.io import read
from soapff_cui.core import make
T = "/Users/harutokono/clay/soapff-cui/examples/minerals/eval"; W = os.path.dirname(os.path.abspath(__file__)); EX = "/Users/harutokono/clay/soapff-cui/examples"
LIB = f"{EX}/minerals/library_clayff_sheet.npz"; BIN = "/Users/harutokono/ScoreMD/.pixi/envs/default/bin"
out = {}
def typed(at, d, qt, tt):
    with contextlib.redirect_stdout(io.StringIO()):
        rep = make(at.copy(), d, LIB, nvt_ps=1, check=True, log=lambda m: None) if False else None
    return None
def run(at, d, qt, tt):
    err = None
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rep = make(at.copy(), d, LIB, nvt_ps=1, check=True, log=lambda m: None)
    except ValueError as e:
        err = str(e)
        with contextlib.redirect_stdout(io.StringIO()):
            rep = make(at.copy(), d, LIB, nvt_ps=1, check=False, log=lambda m: None)
    rows = list(csv.DictReader(open(d + "/charges.csv")))
    bad = [(tt[i], round(float(qt[i]), 4), x["clayff_type"], round(float(x["charge"]), 4)) for i, x in enumerate(rows)
           if x["clayff_type"] != tt[i] or abs(float(x["charge"]) - qt[i]) > 0.01]
    return dict(n=len(rows), bad=bad, stopped=err is not None, fallback=sum(x["source"] == "fallback" for x in rows),
                minsim=min(float(x["similarity"]) for x in rows), net=rep["charge_sum_before"])
# B JAEA
at = read(f"{EX}/jaea/initial.xyz"); f0 = f"{EX}/minerals/refs_clayff/montmorillonite_2W.data"
from soapff_cui.core import read_lammps_data
a0, q0, names = read_lammps_data(f0)
from soapff_cui.clayff import Q
qt = np.array([Q[n] for n in names]); r = run(at, f"{W}/out_old/jaea", qt, names)
out["jaea"] = dict(n=r["n"], wrong=len(r["bad"]), fallback=r["fallback"], minsim=r["minsim"], net=r["net"]); print("jaea", out["jaea"], flush=True)
# C bulk, D interface
for kind, pat in (("bulk", "*_0.extxyz"), ("iface", "*_[0-9]*.extxyz")):
    res = collections.defaultdict(list)
    for f in sorted(glob.glob(f"{T}/tests/{pat}")):
        name = os.path.basename(f).rsplit("_", 1)[0]
        if (kind == "bulk") != f.endswith("_0.extxyz"): continue
        at = read(f); r = run(at, f"{W}/out_old/{os.path.basename(f)[:-7]}", at.arrays["q_true"], [str(x) for x in at.arrays["t_true"]])
        res[name].append(r)
    out[kind] = {k: dict(structures=len(v), atoms=[min(x["n"] for x in v), max(x["n"] for x in v)], wrong=sum(len(x["bad"]) for x in v),
                         structs_wrong=sum(1 for x in v if x["bad"]), stopped=sum(x["stopped"] for x in v),
                         fallback=sum(x["fallback"] for x in v), minsim=min(x["minsim"] for x in v),
                         kinds=dict(collections.Counter(f"{a}({b})->{c}({e})" for x in v for a, b, c, e in x["bad"])))
                 for k, v in res.items()}
    print(kind, json.dumps(out[kind]), flush=True)
json.dump(out, open(f"{W}/results_old_tests.json", "w"), indent=1)
raise SystemExit
# E MD
for lab, src in (("pyrophyllite_001", "pyrophyllite_001_2"), ("brucite_001", "brucite_001_2"), ("quartz_bulk", "quartz_0"), ("albite_bulk", "albite_0"), ("montmorillonite", None)):
    d = f"{W}/md/{lab}"
    with contextlib.redirect_stdout(io.StringIO()):
        make(read(f"{W}/tests/{src}.extxyz") if src else f"{EX}/jaea/initial.xyz", d, LIB, nvt_ps=10, relax_ps=10, log=lambda m: None)
    subprocess.run(f"OMP_NUM_THREADS=1 FI_PROVIDER=tcp {BIN}/mpirun -np 4 {BIN}/lmp -in in.lammps -log log.lammps -screen none", shell=True, cwd=d)
    print("md", lab, os.path.isfile(f"{d}/final.data"), flush=True)
