"""muscovite 2M1 (2x2x1, 672 atoms): 0 K elastic constants, tension x y z and shear xy xz yz at 300 K, with the commands
of soapff_cui.properties (the GUI's ＋ 物性), ClayFF + SOAP charges (library_clayff_sheet.npz)."""
import sys
sys.path.insert(0, "/Users/harutokono/clay/soapff-gui/src")
from soapff_cui import properties as P, core
d = "run"
tid = {}
L = open(f"{d}/system.data").read().splitlines()
i = L.index("Masses") + 2
while L[i].strip():
    tid[L[i].split("#")[1].strip()] = int(L[i].split()[0]); i += 1
props = ["elastic0", "tensile", "shear"]
p = dict(temp=300, press=1.0, relax_ps=10, npt_ps=20, baro="tri", e0_strain=1e-3, tensile_dirs="x y z", tensile_mode="tension",
         tensile_rate=1e9, tensile_max=0.10, shear_planes="xy xz yz", shear_mode="solid", shear_rate=1e9, shear_max=0.10,
         mech_dump_ps=2.0)
src = open(f"{d}/in.lammps").read().splitlines()
k = max(j for j, l in enumerate(src) if l.startswith("minimize"))
src = src[:k + 1] + P.triclinic_lines(d) + P.lines(props, p, tid, core.analysis_lines)
open(f"{d}/in.lammps", "w").write("\n".join(src) + "\n")
print(sum(int(l.split()[1]) for l in src if l.startswith("run ")), "MD steps")
