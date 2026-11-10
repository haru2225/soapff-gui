#!/usr/bin/env python3
"""Orientational order parameter S of the nanoparticles (Underwood & Bourg 2020: S = <P2>, 0 = isotropic, 1 = aligned),
water / clay mass ratio and dry bulk density of a LAMMPS data file written by the stages (needs image flags: write_data keeps them).
    python order_parameter.py 02_dehydrate.data [gouge.json]            (several files: one line each)
"""
import sys, json, re, numpy as np
meta = json.load(open(sys.argv[-1])) if sys.argv[-1].endswith('.json') else {}
files = [f for f in sys.argv[1:] if f.endswith('.data')]
nclay = max(meta.get('clay_types', list(range(1, 19))))
def read(fn):
    L = np.zeros((3, 2)); sec = None; masses = {}; at = []; tilt = np.zeros(3)
    for l in open(fn):
        s = l.split('#')[0].strip()
        if not s: continue
        m = re.match(r'^(\S+) (\S+) ([xyz])lo [xyz]hi', s)
        if m: L['xyz'.index(m.group(3))] = (float(m.group(1)), float(m.group(2))); continue
        if re.match(r'^\S+ \S+ \S+ xy xz yz', s): tilt = np.array(s.split()[:3], float); continue
        if s.startswith('Masses'): sec = 'm'; continue
        if s.startswith('Atoms'): sec = 'a'; continue
        if re.match(r'^[A-Z][A-Za-z ]+$', s): sec = None; continue
        p = s.split()
        if sec == 'm': masses[int(p[0])] = float(p[1])
        elif sec == 'a' and len(p) >= 7: at.append([float(v) for v in p[:7]] + ([float(v) for v in p[7:10]] if len(p) >= 10 else [0, 0, 0]))
    return L, tilt, masses, np.array(at)
for fn in files:
    L, tilt, masses, at = read(fn); Ld = L[:, 1] - L[:, 0]
    typ = at[:, 2].astype(int); mol = at[:, 1].astype(int)
    ix, iy, iz = at[:, 7], at[:, 8], at[:, 9]
    xyz = at[:, 4:7].copy()                                     # unwrapped with the image flags (and the tilt of a sheared box)
    xyz[:, 0] += ix * Ld[0] + iy * tilt[0] + iz * tilt[1]; xyz[:, 1] += iy * Ld[1] + iz * tilt[2]; xyz[:, 2] += iz * Ld[2]
    clay = typ <= nclay; mc = sum(masses[t] for t in typ[clay]); mw = sum(masses[t] for t in typ[~clay])
    V = Ld.prod(); rho = mc * 1.66054 / V
    Q = np.zeros((3, 3)); n = 0
    for m in np.unique(mol[clay]):
        x = xyz[mol == m]; x = x[typ[mol == m] <= nclay]; x = x - x.mean(0)
        w, v = np.linalg.eigh(np.cov(x.T)); nv = v[:, 0]
        Q += 1.5 * np.outer(nv, nv) - 0.5 * np.eye(3); n += 1
    ev, evec = np.linalg.eigh(Q / n)
    print(f'{fn}: particles {n}  S = {ev[-1]:.3f}  director {np.round(evec[:, -1], 2)}  mw/mc = {mw / mc:.3f}  rho_b = {rho:.3f} g/cm3')
