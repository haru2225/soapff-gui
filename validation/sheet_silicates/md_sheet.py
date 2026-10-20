"""Held-out tests and MD of the bulk sheet silicates with library_clayff_sheet.npz (13 references):
 test structures = supercells different from the library's (one more cell along a), truth = bond-graph ClayFF type;
 MD = soapff-cui make (SOAP charges) -> minimise, 10 ps relaxation NVT, 10 ps NVT, 40 ps NPT (tri, 1 atm), 300 K, 0.5 fs.
    python md_sheet.py [names]"""
import sys, os, io, contextlib, subprocess, json, csv
import numpy as np
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui"); sys.path.insert(0, "/Users/harutokono/clay/mincharge")
from ase.io import read, write
from mincharge.ff import ff_types
from soapff_cui.clayff import Q
from soapff_cui.core import make
W = os.path.dirname(os.path.abspath(__file__)); EX = "/Users/harutokono/clay/soapff-cui/examples"
LIB = f"{EX}/minerals/library_clayff_sheet.npz"; BIN = "/Users/harutokono/ScoreMD/.pixi/envs/default/bin"
sys.path.insert(0, f"{EX}/bulk")
NAMES = ("muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "kaolinite", "lizerdite", "talc", "pyrophyllite_2")


def unit(name):
    import make_bulk_sheet as M
    return M.read(name if name != "pyrophyllite_2" else "pyrophyllite")


def check(d, truth):
    rows = list(csv.DictReader(open(f"{d}/charges.csv")))
    bad = [(truth[i], r["clayff_type"]) for i, r in enumerate(rows)
           if r["clayff_type"] != truth[i] or abs(float(r["charge"]) - Q[truth[i]]) > 0.01]
    return dict(n=len(rows), wrong=len(bad), bad=sorted(set(map(str, bad)))[:10],
                fallback=sum(r["source"] == "fallback" for r in rows), minsim=min(float(r["similarity"]) for r in rows))


if __name__ == "__main__":
    names = sys.argv[1:] or NAMES
    for name in names:
        at = unit(name)
        V = at.get_volume()
        w = [V / np.linalg.norm(np.cross(at.cell[(k + 1) % 3], at.cell[(k + 2) % 3])) for k in range(3)]
        rep = [max(1, int(np.ceil(10.0 / x - 1e-6))) for x in w]
        rep[0] += 1                                             # different from the library's supercell
        test = at.repeat(rep); test.wrap()
        os.makedirs(f"{W}/tests", exist_ok=True)
        write(f"{W}/tests/{name}.extxyz", test, format="extxyz")
        truth = ff_types(test)
        d = f"{W}/md/{name}"
        with contextlib.redirect_stdout(io.StringIO()):
            make(test, d, LIB, nvt_ps=10, npt_ps=40, baro="tri", relax_ps=10, dump_ps=1.0, log=lambda m: None)
        res = dict(name=name, atoms=len(test), repeat=rep, cell0=test.cell.cellpar().tolist(), structure=check(d, truth))
        subprocess.run(f"OMP_NUM_THREADS=1 FI_PROVIDER=tcp {BIN}/mpirun -np 2 {BIN}/lmp -in in.lammps -log log.lammps "
                       "-screen none", shell=True, cwd=d)
        res["md_done"] = os.path.isfile(f"{d}/final.data")
        json.dump(res, open(f"{d}/eval.json", "w"), indent=1)
        print(json.dumps({k: v for k, v in res.items() if k != "cell0"}), flush=True)
