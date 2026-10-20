"""Results of md_sheet.py: (A) SOAP typing of the last MD frame (300 K, after NPT; truth = the bond-graph ClayFF type of
the test structure, same atom order); (D) lattice parameters averaged over the last 30 ps of the NPT (tri, 1 atm,
300 K) per unit cell, against the input (experimental) structure.  -> results_sheet.json"""
import sys, os, io, contextlib, json, glob
import numpy as np
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui"); sys.path.insert(0, "/Users/harutokono/clay/mincharge")
sys.path.insert(0, "/Users/harutokono/clay/soapff-gui/src")
from ase import Atoms
from ase.io import read
from mincharge.ff import ff_types
from soapff_cui.core import make
from soapff_cui.sova_structure import frames
W = os.path.dirname(os.path.abspath(__file__)); LIB = "/Users/harutokono/clay/soapff-cui/examples/minerals/library_clayff_sheet.npz"
import md_sheet as M

out = {}
for d in sorted(glob.glob(f"{W}/md/*/")):
    name = os.path.basename(d.rstrip("/"))
    if not os.path.isfile(f"{d}/eval.json") or not os.path.isfile(f"{d}/traj.lammpstrj"):
        continue
    ev = json.load(open(f"{d}/eval.json"))
    test = read(f"{W}/tests/{name}.extxyz")
    truth = ff_types(test)
    fr = list(frames(f"{d}/traj.lammpstrj"))
    st, cell, pos = fr[-1]
    last = Atoms(test.get_chemical_symbols(), positions=pos, cell=cell, pbc=True)
    dd = f"{W}/frame/{name}"
    with contextlib.redirect_stdout(io.StringIO()):
        make(last, dd, LIB, nvt_ps=1, log=lambda m: None)
    ev["md_frame"] = M.check(dd, truth)
    # lattice: frames every 1 ps; production = 10 ps NVT + 40 ps NPT; average t >= 20 ps
    rep = np.array(ev["repeat"], float)
    cp = np.array([Atoms(cell=c, pbc=True).cell.cellpar() for s, c, p in fr if s * 0.5 / 1000.0 >= 20.0 - 1e-9])
    mean = cp.mean(axis=0); mean[:3] /= rep
    ref = np.array(ev["cell0"]); ref[:3] /= rep
    ev["lattice_md"] = mean.round(3).tolist(); ev["lattice_ref"] = ref.round(3).tolist()
    ev["lattice_dev_percent"] = ((mean[:3] - ref[:3]) / ref[:3] * 100).round(2).tolist()
    ev["angle_dev_deg"] = (mean[3:] - ref[3:]).round(2).tolist()
    ev["volume_dev_percent"] = round(float((np.prod(mean[:3]) / np.prod(ref[:3]) * np.sqrt(1 - sum(np.cos(np.radians(mean[3:]))**2) + 2*np.prod(np.cos(np.radians(mean[3:])))) / np.sqrt(1 - sum(np.cos(np.radians(ref[3:]))**2) + 2*np.prod(np.cos(np.radians(ref[3:])))) - 1) * 100), 2)
    out[name] = ev
    print(f"{name:15s} test {ev['structure']['wrong']}/{ev['structure']['n']} wrong, frame {ev['md_frame']['wrong']}/{ev['md_frame']['n']} "
          f"(min sim {ev['md_frame']['minsim']:.3f}, fallback {ev['md_frame']['fallback']}); a b c {ev['lattice_dev_percent']} %, "
          f"angles {ev['angle_dev_deg']} deg, V {ev['volume_dev_percent']} %", flush=True)
json.dump(out, open(f"{W}/results_sheet.json", "w"), indent=1)
