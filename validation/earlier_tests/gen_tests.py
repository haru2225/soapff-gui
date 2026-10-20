"""test structures: (1) pyrophyllite (001) and brucite (001) interfaces, slab + water placed at random (>= 3 A from the
mineral, centres >= 2 A above the surface), seeds 2-21; no termination needed (no bond is cut);
(2) bulk crystals in supercells different from the library's.  Truth = bond-graph ClayFF type, ClayFF charge."""
import sys, os, math, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui"); sys.path.insert(0, "/Users/harutokono/clay/soapff-cui/examples/minerals")
sys.path.insert(0, "/Users/harutokono/clay/mincharge")
W = os.path.dirname(os.path.abspath(__file__)); os.makedirs(f"{W}/tests", exist_ok=True)
def one(job):
    kind, name, seed = job
    f = f"{W}/tests/{name}_{seed}.extxyz"
    if os.path.isfile(f): return f
    import build_ase as B
    from ase.io import write
    from mincharge.ff import ff_types
    from soapff_cui.clayff import Q
    if kind == "iface":
        mineral, hkl = B.SYSTEMS[name]; B.TET_AL = False
        slab, _ = B.cut_slab(B.bulk(mineral), hkl)
        rep = [max(1, math.ceil(10.0 / np.linalg.norm(slab.cell[k]))) for k in (0, 1)]; slab = slab.repeat((rep[0], rep[1], 1))
        at, info = B.add_water(slab, seed=seed, tol_min=3.0, z_off=2.0)
    else:
        rep = {"quartz": (4, 3, 2), "albite": (3, 2, 1), "brucite": (5, 5, 4), "pyrophyllite": (3, 2, 3)}[name]
        at = B.bulk(name).repeat(rep)
    t = ff_types(at)
    at.set_array("q_true", np.array([Q[x] for x in t], float)); at.set_array("t_true", np.array(t, dtype="U8"))
    write(f + ".tmp.extxyz", at, format="extxyz"); os.rename(f + ".tmp.extxyz", f)
    print(name, seed, len(at), flush=True); return f
if __name__ == "__main__":
    jobs = [("iface", n, s) for n in ("pyrophyllite_001", "brucite_001") for s in range(2, 22)]
    jobs += [("bulk", n, 0) for n in ("quartz", "albite", "brucite", "pyrophyllite")]
    with Pool(8) as p: p.map(one, jobs, chunksize=1)
