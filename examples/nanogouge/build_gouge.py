#!/usr/bin/env python3
"""Build a disordered assemblage of N edge-bearing Na-montmorillonite nanoparticles in water (LAMMPS data file).

After Underwood & Bourg, J. Phys. Chem. C 124, 3702 (2020): particles are added one by one with a random position and
orientation (no overlap) in a periodic box, NCONF independent configurations are generated, and the most randomly oriented
one (lowest order parameter S) is kept; it is then solvated at a water-to-clay mass ratio m_w/m_c (paper: about 3).

The particle is the one soapff-gui writes for the optional library "Montmorillonite edge"
(../montmorillonite_edge/out_clay_montmorillonite/system.data); its ClayFF/Lammers edge force field is copied as it is.
Water is the flexible SPC water of ClayFF (Cygan et al. 2004). Needs numpy and scipy only.

    python build_gouge.py -n 12 --mwmc 3.0 --nconf 300 --seed 1 -o gouge_N12
"""
import argparse, json, re, sys
import numpy as np
from scipy.spatial import cKDTree

p = argparse.ArgumentParser()
p.add_argument('--particle', default='../montmorillonite_edge/out_clay_montmorillonite', help='folder with system.data and system.settings of one particle')
p.add_argument('-n', type=int, default=12, help='number of particles')
p.add_argument('--mwmc', type=float, default=3.0, help='initial water / clay mass ratio (paper: ~3)')
p.add_argument('--aspect', type=float, default=40 / 24, help='Lz / Lx = Lz / Ly of the box (paper: 40 nm / 24 nm)')
p.add_argument('--nconf', type=int, default=300, help='number of random configurations; the lowest-S one is kept (paper: 300)')
p.add_argument('--mind', type=float, default=3.5, help='minimum atom-atom distance between two particles (A)')
p.add_argument('--water-mind', type=float, default=2.6, help='minimum distance water O to a clay atom (A)')
p.add_argument('--clay-density', type=float, default=2.7, help='g/cm3, to size the box')
p.add_argument('--seed', type=int, default=1)
p.add_argument('-o', '--out', default='gouge')
a = p.parse_args()
rng = np.random.default_rng(a.seed)

# ---------------------------------------------------------------- read the particle (LAMMPS data, atom_style full)
txt = open(f'{a.particle}/system.data').read().split('\n')
sec = None; masses = {}; atoms = []; bonds = []; angles = []; box = {}
for l in txt:
    s = l.strip()
    if not s: continue
    m = re.match(r'^(\S+)\s+(\S+)\s+([xyz])lo\s+[xyz]hi', s)
    if m: box[m.group(3)] = (float(m.group(1)), float(m.group(2))); continue
    if s.startswith('Masses'): sec = 'm'; continue
    if s.startswith('Atoms'): sec = 'a'; continue
    if s.startswith('Bonds'): sec = 'b'; continue
    if s.startswith('Angles'): sec = 'g'; continue
    if s.startswith('Velocities'): sec = None; continue
    pp = s.split()
    if sec == 'm': masses[int(pp[0])] = (float(pp[1]), s.split('#')[1].strip() if '#' in s else '')
    elif sec == 'a' and len(pp) >= 7: atoms.append((int(pp[0]), int(pp[2]), float(pp[3]), float(pp[4]), float(pp[5]), float(pp[6])))
    elif sec == 'b' and len(pp) >= 4: bonds.append(tuple(int(v) for v in pp[:4]))
    elif sec == 'g' and len(pp) >= 5: angles.append(tuple(int(v) for v in pp[:5]))
atoms.sort()
ptype = np.array([x[1] for x in atoms]); pq = np.array([x[2] for x in atoms]); pos0 = np.array([[x[3], x[4], x[5]] for x in atoms])
natp = len(atoms); ntp = max(masses)
n_bt = max(b[1] for b in bonds); n_at = max(g[1] for g in angles)
mass_p = sum(masses[t][0] for t in ptype)
settings = open(f'{a.particle}/system.settings').read().split('\n')
coef = [l for l in settings if re.match(r'^(pair_coeff|bond_coeff|angle_coeff)\s', l)]
# the particle sits whole in its box: centre it, find its disk radius and thickness
cen = pos0.mean(0); pos0 = pos0 - cen
w, v = np.linalg.eigh(np.cov(pos0.T)); normal0 = v[:, 0]
pos0 = pos0 @ v[:, ::-1]                      # x,y in the plane of the layer, z = normal
R = np.sqrt((pos0[:, :2] ** 2).sum(1)).max(); thick = np.ptp(pos0[:, 2])
nNa = int((np.array([masses[t][1] for t in ptype]) == 'Na').sum())
print(f'particle: {natp} atoms, mass {mass_p:.0f} amu, radius {R:.1f} A, thickness {thick:.1f} A, {nNa} Na, {len(bonds)} bonds, {len(angles)} angles')

# ---------------------------------------------------------------- water (flexible SPC, ClayFF)
T_HW, T_OW = ntp + 1, ntp + 2; B_W = n_bt + 1; A_W = n_at + 1
M_H, M_O = 1.00794, 15.9994; Q_H, Q_O = 0.41, -0.82; M_W = 2 * M_H + M_O
r_oh, th = 1.0, np.deg2rad(109.47)
wat0 = np.array([[0, 0, 0], [r_oh * np.sin(th / 2), 0, r_oh * np.cos(th / 2)], [-r_oh * np.sin(th / 2), 0, r_oh * np.cos(th / 2)]])
wat0 -= (wat0 * np.array([M_O, M_H, M_H])[:, None]).sum(0) / M_W          # about the centre of mass

# ---------------------------------------------------------------- box
n_w = int(round(a.mwmc * a.n * mass_p / M_W))
V = a.n * mass_p * 1.66054 / a.clay_density + n_w * 29.9                    # A^3  (1 amu/(g/cm3) = 1.66054 A^3; water 29.9 A^3)
Lx = (V / a.aspect) ** (1 / 3); L = np.array([Lx, Lx, a.aspect * Lx])
print(f'N={a.n}: {n_w} water, box {L[0]:.1f} x {L[1]:.1f} x {L[2]:.1f} A, {(a.n * natp + 3 * n_w)} atoms')

def rotmat(q):
    q = q / np.linalg.norm(q); w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])

def mic(d): return d - L * np.round(d / L)

def order_S(Rms):
    Q = np.zeros((3, 3))
    for Rm in Rms:
        n = Rm[:, 2]; Q += 1.5 * np.outer(n, n) - 0.5 * np.eye(3)
    return np.linalg.eigvalsh(Q / len(Rms))[-1]

best = None
for k in range(a.nconf):
    out = []; Rms = []
    ok = True
    for i in range(a.n):
        for t in range(5000):
            c = rng.random(3) * L; Rm = rotmat(rng.normal(size=4))
            xyz = (pos0 @ Rm.T + c) % L
            tr = cKDTree(xyz, boxsize=L); good = True
            for (c2, x2, tr2, _) in out:
                if np.linalg.norm(mic(c - c2)) > 2 * R + 2: continue
                if tr.count_neighbors(tr2, a.mind) > 0: good = False; break
            if good: out.append((c, xyz, tr, Rm)); Rms.append(Rm); break
        else: ok = False; break
    if not ok: continue
    S = order_S(Rms)
    if best is None or S < best[0]: best = (S, out, Rms, k)
    if (k + 1) % 50 == 0: print(f'  config {k + 1}: best S = {best[0]:.3f}' if best else f'  config {k + 1}: none placed')
if best is None: sys.exit('could not place the particles: use fewer particles, a larger --mwmc (bigger box) or --mind smaller')
S, out, Rms, kb = best
print(f'kept configuration {kb} with S = {S:.3f}')
clay_xyz = np.vstack([o[1] for o in out])

# ---------------------------------------------------------------- water on a jittered grid, random orientation
tree_c = cKDTree(clay_xyz, boxsize=L)
samp = rng.random((200000, 3)) * L; dd, _ = tree_c.query(samp); ffree = (dd > a.water_mind).mean()
spacing = (ffree * L.prod() / n_w) ** (1 / 3)
print(f'free volume fraction {ffree:.3f}; water grid spacing {spacing:.2f} A (liquid water: 3.10 A)')
if spacing < 2.75: sys.exit('too much water for the free volume: lower --mwmc or --water-mind')
ng = np.ceil(L / (spacing * 0.93)).astype(int); g = (np.indices(ng).reshape(3, -1).T + 0.5) * (L / ng)
g = (g + (rng.random(g.shape) - 0.5) * 0.25 * (L / ng)) % L
dd, _ = tree_c.query(g); g = g[dd > a.water_mind]
if len(g) < n_w: sys.exit(f'only {len(g)} water sites for {n_w} molecules')
g = g[rng.permutation(len(g))[:n_w]]
wxyz = np.empty((n_w, 3, 3))
for i in range(n_w): wxyz[i] = (wat0 @ rotmat(rng.normal(size=4)).T) + g[i]

# ---------------------------------------------------------------- write data + settings
ntypes = ntp + 2; nat = a.n * natp + 3 * n_w
with open(f'{a.out}.data', 'w') as f:
    f.write(f'{a.n} Na-montmorillonite edge particles + {n_w} flexible SPC water (build_gouge.py), S = {S:.3f}\n\n')
    f.write(f'{nat} atoms\n{ntypes} atom types\n{a.n * len(bonds) + 2 * n_w} bonds\n{n_bt + 1} bond types\n{a.n * len(angles) + n_w} angles\n{n_at + 1} angle types\n\n')
    for d, lab in zip(range(3), 'xyz'): f.write(f'0.0 {L[d]:.6f} {lab}lo {lab}hi\n')
    f.write('\nMasses\n\n')
    for t in range(1, ntp + 1): f.write(f'{t} {masses[t][0]:.5f}  # {masses[t][1]}\n')
    f.write(f'{T_HW} {M_H}  # hw\n{T_OW} {M_O}  # ow\n\nAtoms # full\n\n')
    aid = 0
    def wr(i, mol, t, q, x):
        f.write(f'{i} {mol} {t} {q:.8f} {x[0]:.5f} {x[1]:.5f} {x[2]:.5f} 0 0 0\n')
    for pi, (c, xyz, tr, Rm) in enumerate(out):
        # keep each particle whole: unwrapped positions about its centre, placed inside the box by shifting the whole particle
        full = pos0 @ Rm.T + c
        for j in range(natp): wr(pi * natp + j + 1, pi + 1, int(ptype[j]), pq[j], full[j])
    base = a.n * natp
    for i in range(n_w):
        for j, (t, q) in enumerate(((T_OW, Q_O), (T_HW, Q_H), (T_HW, Q_H))): wr(base + 3 * i + j + 1, a.n + 1 + i, t, q, wxyz[i, j])
    f.write('\nBonds\n\n'); bi = 0
    for pi in range(a.n):
        for b in bonds: bi += 1; f.write(f'{bi} {b[1]} {b[2] + pi * natp} {b[3] + pi * natp}\n')
    for i in range(n_w):
        o = base + 3 * i + 1
        bi += 1; f.write(f'{bi} {B_W} {o} {o + 1}\n'); bi += 1; f.write(f'{bi} {B_W} {o} {o + 2}\n')
    f.write('\nAngles\n\n'); gi = 0
    for pi in range(a.n):
        for g_ in angles: gi += 1; f.write(f'{gi} {g_[1]} {g_[2] + pi * natp} {g_[3] + pi * natp} {g_[4] + pi * natp}\n')
    for i in range(n_w):
        o = base + 3 * i + 1; gi += 1; f.write(f'{gi} {A_W} {o + 1} {o} {o + 2}\n')
# the data file keeps the particles unwrapped about their centres; LAMMPS wraps them with `change_box all ... remap` / run start
names = {t: masses[t][1] for t in masses}
with open(f'{a.out}.settings', 'w') as f:
    f.write('# force field: Na-montmorillonite edge particle (Underwood & Bourg 2020, after Lammers et al.; ClayFF) + flexible SPC water (ClayFF, Cygan et al. 2004)\n')
    f.write('# generated by build_gouge.py; the pair_coeff / bond_coeff / angle_coeff lines of the particle are those of soapff-gui\n')
    for l in coef: f.write(l + '\n')
    f.write(f'pair_coeff {T_HW} {T_HW} 0.0 0.0  # hw\npair_coeff {T_OW} {T_OW} 0.1554 3.1655413252962536  # ow\n')
    f.write(f'bond_coeff {B_W} 554.1349 1.0  # O-H (water)\nangle_coeff {A_W} 45.7696 109.47  # HOH (water)\n')
json.dump(dict(n_particles=a.n, n_water=n_w, atoms_per_particle=natp, n_atom_types=ntypes, type_hw=T_HW, type_ow=T_OW, clay_types=list(range(1, ntp + 1)),
               mass_clay_amu=a.n * mass_p, S_initial=float(S), box=L.tolist(), seed=a.seed, nconf=a.nconf, mwmc=a.mwmc), open(f'{a.out}.json', 'w'), indent=1)
print(f'wrote {a.out}.data, {a.out}.settings, {a.out}.json')
