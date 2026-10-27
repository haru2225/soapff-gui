"""Why the mineral is typed first: the mineral atoms of the last frames of the 11 water-(0 0 1) interface MDs (refs/if_*/), typed with the
water-free mineral library from their SOAP computed (a) on the whole structure with the water, (b) on the structure without the water
molecules (the ions stay).  Counts the atoms below the similarity 0.95 and the atoms with a wrong type (the true types are those of the MD).

    python mineral_in_water.py     (with the Python of soapff-gui)  ->  mineral_in_water.json"""
import json
import os

import numpy as np

from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "refs")
lib = core.load_library(os.path.join(HERE, "..", "..", "src", "soapff_gui", "data", "library_clayff_bulk.npz"))
X=lib["X"].astype(np.float32); el=lib["element"]; ty=lib["type"]
tot=dict(n=0,low_with=0,wrong_with=0,low_wo=0,wrong_wo=0); per={}
for m in ["montmorillonite_dry_Na","muscovite_2M1","muscovite_1M","phlogopite_1M","nanpingite","pyrophyllite","pyrophyllite_2","kaolinite","lizerdite","talc","brucite"]:
    d=f"{R}/if_{m}"
    at0,q0,t0=core.read_lammps_data(f"{d}/system.data"); t0=np.array([core.TYPE_ALIAS.get(str(x),str(x)) for x in t0])
    at=core.read_dump_frames(f"{d}/traj.lammpstrj", at0.get_chemical_symbols())[-1]
    sym=np.array(at.get_chemical_symbols())
    mineral=~np.isin(t0,["ow","hw"])
    res={}
    for tag,sel in (("with",np.ones(len(at),bool)),("wo",mineral)):
        sub=at[sel]
        F=core.soap_vectors(sub,lib["species"],**lib["soap"]).astype(np.float32)
        idx=np.nonzero(mineral[sel])[0]
        e=np.array(sub.get_chemical_symbols())
        best=np.zeros(len(idx)); got=np.empty(len(idx),dtype=object)
        for k,i in enumerate(idx):
            mm=el==e[i]
            s=X[mm]@F[i]; j=s.argmax(); best[k]=s[j]; got[k]=ty[mm][j]
        truth=t0[sel][idx]
        res[tag]=(int((best<0.95).sum()),int((got!=truth).sum()),float(best.min()))
    n=int(mineral.sum()); per[m]=dict(n=n,**{f"{k}":v for k,v in res.items()})
    tot["n"]+=n; tot["low_with"]+=res["with"][0]; tot["wrong_with"]+=res["with"][1]; tot["low_wo"]+=res["wo"][0]; tot["wrong_wo"]+=res["wo"][1]
    print(m,per[m],flush=True)
print(tot)
json.dump(dict(per=per, tot=tot), open(os.path.join(HERE, "mineral_in_water.json"), "w"), indent=1)
