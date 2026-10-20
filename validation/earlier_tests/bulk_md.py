"""MD of the bulk-crystal test structures (supercells different from the library's): minimise, 10 ps relaxation NVT,
10 ps NVT, 40 ps NPT (all cell lengths and tilts, 1 atm), 300 K, 0.5 fs."""
import sys, os, io, contextlib, subprocess, json
sys.path.insert(0, "/Users/harutokono/clay/soapff-cui"); sys.path.insert(0, "/Users/harutokono/clay/mincharge")
from ase.io import read
from soapff_cui.core import make
W = os.path.dirname(os.path.abspath(__file__)); LIB = "/Users/harutokono/clay/soapff-cui/examples/minerals/library_clayff_bulk.npz"
BIN = "/Users/harutokono/ScoreMD/.pixi/envs/default/bin"
for name in ("quartz", "albite", "brucite", "pyrophyllite"):
    d = f"{W}/mdbulk/{name}"
    with contextlib.redirect_stdout(io.StringIO()):
        make(read(f"{W}/tests/{name}_0.extxyz"), d, LIB, nvt_ps=10, npt_ps=40, baro="tri", relax_ps=10, dump_ps=1.0, log=lambda m: None)
    subprocess.run(f"OMP_NUM_THREADS=1 FI_PROVIDER=tcp {BIN}/mpirun -np 2 {BIN}/lmp -in in.lammps -log log.lammps -screen none", shell=True, cwd=d)
    print("md", name, os.path.isfile(f"{d}/final.data"), flush=True)
