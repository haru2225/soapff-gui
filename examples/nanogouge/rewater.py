#!/usr/bin/env python3
"""A fabric (the clay of a LAMMPS data file, e.g. RUN_N54_fab_ori/02_dehydrate.data) with a new amount of water: all its water is taken
out, then m_w/m_c = --mwmc of water is put into the space between the particles, by packmol (default when found) or on a grid.
  python3 rewater.py RUN_N54_fab_ori/02_dehydrate.data --mwmc 0.1 -o RUN_N54_ori_w0.1/gouge.data
  --method packmol   packmol (PACKMOL=/path/packmol, else ./bin/packmol (get_packmol.sh), else packmol on the PATH): the molecules anywhere in
                     the free space, at least --tol A from every atom (also across the periodic boundaries)
  --method grid      no packmol: liquid-density grid sites at least --mind A from the clay, the sites nearest to the clay first
                     (the water on the surfaces and in the interlayers first); numpy only
  --method auto      packmol when found, else grid (default)
The water molecule (types, charges, geometry, bond and angle types) is taken from the water of the data file itself.  The clay atoms (with Na)
keep their positions and image flags; the atoms are numbered again (clay first); no velocities (in.all gives them).  The box must be
orthogonal (a fabric made by 02 is).  Writes the data file and <data>.json (what was done).
"""
import argparse, json, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('data')
ap.add_argument('--mwmc', type=float, required=True, help='water / clay mass ratio of the new structure (0 = dry)')
ap.add_argument('-o', '--out', required=True)
ap.add_argument('--method', choices=('auto', 'packmol', 'grid'), default='auto')
ap.add_argument('--tol', type=float, default=2.0, help='packmol tolerance (A)')
ap.add_argument('--mind', type=float, default=2.6, help='grid: smallest distance water O - clay atom (A)')
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--tow', type=int, default=0, help='type of water O (default: from the Masses comments, else the highest type)')
ap.add_argument('--thw', type=int, default=0, help='type of water H (default: from the Masses comments, else the second highest type)')
a = ap.parse_args()
class _OldRng:                                    # numpy < 1.17 (the python3 of a login node) has no default_rng
    def __init__(self, seed): self.r = np.random.RandomState(seed)
    def random(self, size=None): return self.r.random_sample(size)
    def normal(self, size=None): return self.r.normal(size=size)
    def permutation(self, n): return self.r.permutation(n)
rng = np.random.default_rng(a.seed) if hasattr(np.random, 'default_rng') else _OldRng(a.seed)
M_WATER = 18.0153                                      # as in.all: n_water = m_w/m_c * mass(clay) / 18.0153

SECTIONS = ('Masses', 'Pair Coeffs', 'PairIJ Coeffs', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs', 'Improper Coeffs', 'Atoms', 'Velocities',
            'Bonds', 'Angles', 'Dihedrals', 'Impropers')
# ---------------------------------------------------------------- read
counts = {}; box = {}; tilt = None; sec = {}; order = []; cur = None; title = ''
with open(a.data) as f:
    title = f.readline().strip()
    for line in f:
        s = line.split('#')[0].strip()
        head = line.strip()
        name = next((n for n in SECTIONS if head == n or head.startswith(n + ' ') or head.startswith(n + '#')), None)
        if name:
            cur = name; sec[name] = {'head': head, 'lines': []}; order.append(name); continue
        if not s: continue
        if cur is None:
            m = re.match(r'^(\d+)\s+(atoms|bonds|angles|dihedrals|impropers|atom types|bond types|angle types|dihedral types|improper types)$', s)
            if m: counts[m.group(2)] = int(m.group(1)); continue
            m = re.match(r'^([-\d.eE+]+)\s+([-\d.eE+]+)\s+([xyz])lo\s+[xyz]hi$', s)
            if m: box[m.group(3)] = (float(m.group(1)), float(m.group(2))); continue
            m = re.match(r'^([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+xy\s+xz\s+yz$', s)
            if m: tilt = tuple(float(m.group(i)) for i in (1, 2, 3)); continue
            continue
        sec[cur]['lines'].append(line.rstrip('\n'))
if 'Atoms' not in sec: sys.exit(f'{a.data}: no Atoms section')
if tilt and any(abs(t) > 1e-6 for t in tilt): sys.exit(f'{a.data}: a triclinic box (xy xz yz {tilt}); rewater.py needs an orthogonal one')
lo = np.array([box[k][0] for k in 'xyz']); L = np.array([box[k][1] - box[k][0] for k in 'xyz'])
masses = {}; wcomment = {}
for line in sec['Masses']['lines']:
    w = line.split('#')[0].split()
    masses[int(w[0])] = float(w[1])
    if '#' in line: wcomment[line.split('#')[1].strip().lower()] = int(w[0])
tow = a.tow or wcomment.get('ow') or max(masses)
thw = a.thw or wcomment.get('hw') or sorted(masses)[-2]
wtypes = {tow, thw}

# atoms: id mol type q x y z [ix iy iz]
A = []
for line in sec['Atoms']['lines']:
    w = line.split('#')[0].split()
    img = tuple(int(x) for x in w[7:10]) if len(w) >= 10 else (0, 0, 0)
    A.append((int(w[0]), int(w[1]), int(w[2]), float(w[3]), float(w[4]), float(w[5]), float(w[6])) + img)
A.sort()
ids = np.array([r[0] for r in A]); mol = np.array([r[1] for r in A]); typ = np.array([r[2] for r in A]); q = np.array([r[3] for r in A])
xyz = np.array([r[4:7] for r in A]); img = np.array([r[7:10] for r in A])
iswat = np.isin(typ, list(wtypes))
def readcon(name):
    return [tuple(int(x) for x in l.split('#')[0].split()) for l in sec[name]['lines']] if name in sec else []
bonds, angles = readcon('Bonds'), readcon('Angles')
dihed, impro = readcon('Dihedrals'), readcon('Impropers')
pos_of = {i: k for k, i in enumerate(ids)}

# the water molecule of this data file: O, H, H with charges and relative positions; its bond and angle types
tmpl = None
if iswat.any():
    io = int(np.nonzero(typ == tow)[0][0]); oid = ids[io]
    hb = [b for b in bonds if oid in b[2:4]]
    hs = [b[2] if b[3] == oid else b[3] for b in hb]
    ang = [an for an in angles if oid in an[2:5]]
    if len(hs) == 2 and ang:
        rel = np.array([xyz[pos_of[h]] - xyz[io] for h in hs]); rel -= L * np.round(rel / L)
        tmpl = dict(q=(q[io], q[pos_of[hs[0]]], q[pos_of[hs[1]]]), bt=hb[0][1], at=ang[0][1],
                    geo=np.vstack([np.zeros(3), rel]))
mclay = sum(masses[t] for t in typ[~iswat])
n_old = int(iswat.sum() // 3)
n_w = int(round(a.mwmc * mclay / M_WATER))
print(f'{a.data}: {len(ids)} atoms, {n_old} water molecules (m_w/m_c {n_old * M_WATER / mclay:.3f}), box {L[0]:.1f} x {L[1]:.1f} x {L[2]:.1f} A')
print(f'new: {n_w} water molecules (m_w/m_c {a.mwmc})')
if n_w > 0 and tmpl is None: sys.exit('no water molecule in the data file to take as the template (types, charges, bond and angle types)')

# ---------------------------------------------------------------- the clay alone (wrapped into the box: the free space)
cid = ids[~iswat]
cxyz = xyz[~iswat] - lo; cxyz -= L * np.floor(cxyz / L)                 # 0 <= x < L
newid = np.zeros(ids.max() + 1, dtype=np.int64); newid[cid] = np.arange(1, len(cid) + 1)
keep = lambda con, k0: [c for c in con if all(newid[x] > 0 for x in c[k0:])]
bonds_c, angles_c, dihed_c, impro_c = keep(bonds, 2), keep(angles, 2), keep(dihed, 2), keep(impro, 2)

def rand_rot(n):
    qq = rng.normal(size=(n, 4)); qq /= np.linalg.norm(qq, axis=1)[:, None]; w, x, y, z = qq.T
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
                     np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
                     np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], 1)

def find_packmol():
    for c in (os.environ.get('PACKMOL'), os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin', 'packmol'), shutil.which('packmol')):
        if c and os.path.isfile(c) and os.access(c, os.X_OK): return c
    return None

def place_packmol(exe):
    """packmol: the clay fixed (with its periodic images near the faces), the water inside the box less tol/2 at each face"""
    m = a.tol + 0.5
    imgs = [cxyz]
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                if sx == sy == sz == 0: continue
                p = cxyz + np.array([sx, sy, sz]) * L
                sel = np.all((p > -m) & (p < L + m), axis=1)
                if sel.any(): imgs.append(p[sel])
    fixed = np.vstack(imgs)
    d = tempfile.mkdtemp(prefix='rewater_', dir=os.path.dirname(os.path.abspath(a.out)) or '.')
    with open(os.path.join(d, 'clay.xyz'), 'w') as f:
        f.write(f'{len(fixed)}\nclay (fixed)\n'); f.writelines(f'C {x:.4f} {y:.4f} {z:.4f}\n' for x, y, z in fixed)
    g = tmpl['geo'] - tmpl['geo'].mean(0)
    with open(os.path.join(d, 'water.xyz'), 'w') as f:
        f.write('3\nwater\n'); f.writelines(f'{e} {x:.4f} {y:.4f} {z:.4f}\n' for e, (x, y, z) in zip('OHH', g))
    h = a.tol / 2
    with open(os.path.join(d, 'pack.inp'), 'w') as f:
        f.write(f'tolerance {a.tol}\nfiletype xyz\noutput out.xyz\nseed {a.seed}\n'
                f'structure clay.xyz\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n'
                f'structure water.xyz\n  number {n_w}\n  inside box {h} {h} {h} {L[0] - h:.4f} {L[1] - h:.4f} {L[2] - h:.4f}\nend structure\n')
    print(f'packmol ({exe}) in {d}: {n_w} water molecules, {len(fixed)} fixed atoms (clay and its images); log: {d}/packmol.log')
    with open(os.path.join(d, 'pack.inp')) as fi, open(os.path.join(d, 'packmol.log'), 'w') as fo:
        rc = subprocess.call([exe], stdin=fi, stdout=fo, stderr=subprocess.STDOUT, cwd=d)
    out = os.path.join(d, 'out.xyz')
    if not os.path.isfile(out) and os.path.isfile(out + '_FORCED'):
        out += '_FORCED'; print('packmol: no full solution; its best one is taken (the minimisation of in.all removes the small overlaps)')
    if not os.path.isfile(out): sys.exit(f'packmol failed (status {rc}): see {d}/packmol.log')
    with open(out) as f:
        lines = f.readlines()[2:]
    w = np.array([[float(v) for v in l.split()[1:4]] for l in lines[len(fixed):len(fixed) + 3 * n_w]]).reshape(n_w, 3, 3)
    # the fixed clay must be where it was (packmol keeps a "fixed" structure without "center" in place)
    c0 = np.array([float(v) for v in lines[0].split()[1:4]])
    if np.abs(c0 - fixed[0]).max() > 1e-3: sys.exit('packmol moved the fixed clay: not expected; see ' + d)
    with open(os.path.join(d, 'packmol.log')) as f: log = f.read()
    ok = 'Success!' in log
    print('packmol: ' + ('converged (every distance at least the tolerance)' if ok else 'not fully converged; the best solution is used'))
    shutil.rmtree(d) if ok else None
    return w, ('packmol' if ok else 'packmol (not converged)')

def place_grid():
    """liquid-density grid sites (3.10 A, jittered) at least mind from the clay, the sites nearest to the clay first; numpy only (voxels)"""
    h = 0.5; nv = np.maximum(1, np.round(L / h).astype(int)); hv = L / nv
    lev = np.full(nv, 255, dtype=np.uint8)
    iv = np.floor(cxyz / hv).astype(np.int64) % nv
    shells = (a.mind, a.mind + 0.6, a.mind + 1.2, a.mind + 2.0, a.mind + 3.5)
    for k, R in reversed(list(enumerate(shells))):
        r = int(math.ceil((R + 1.0) / hv.min()))
        o = np.indices((2 * r + 1,) * 3).reshape(3, -1).T - r
        Re = R + (np.sqrt((hv ** 2).sum()) if k == 0 else 0.0)     # the free space: no site closer than mind to a clay atom (voxel error)
        o = o[((o * hv) ** 2).sum(1) <= Re * Re]
        for dv in o:
            j = (iv + dv) % nv
            lev[j[:, 0], j[:, 1], j[:, 2]] = k
    ng = np.ceil(L / 3.10).astype(int); g = (np.indices(ng).reshape(3, -1).T + 0.5) * (L / ng)
    g = (g + (rng.random(g.shape) - 0.5) * 0.15 * (L / ng)) % L
    gl = lev[tuple((np.floor(g / hv).astype(np.int64) % nv).T)].astype(float)
    ok = gl > 0; g, gl = g[ok], gl[ok]
    if len(g) < n_w: sys.exit(f'grid: only {len(g)} water sites in the free space for {n_w} molecules: less water (--mwmc) or packmol')
    o = np.argsort(gl + rng.random(len(gl)) * 0.5)[:n_w]
    print(f'grid: {len(g)} free sites, {n_w} taken, the nearest to the clay first (up to shell {int(gl[o].max())} of {len(shells)}: '
          f'{shells[min(int(gl[o].max()), len(shells) - 1)]:.1f} A)')
    geo = tmpl['geo'] - tmpl['geo'].mean(0)
    return np.einsum('nij,aj->nai', rand_rot(n_w), geo) + g[o][:, None, :], 'grid'

os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
wpos = np.zeros((0, 3, 3)); method = 'dry'
if n_w > 0:
    exe = find_packmol() if a.method in ('auto', 'packmol') else None
    if a.method == 'packmol' and not exe: sys.exit('no packmol: PACKMOL=/path/packmol, ./get_packmol.sh, or --method grid')
    if exe: wpos, method = place_packmol(exe)
    else:
        if a.method == 'auto': print('(no packmol: the grid method; ./get_packmol.sh makes packmol in ./bin)')
        wpos, method = place_grid()

# ---------------------------------------------------------------- write
nc = len(cid); molmax = int(mol[~iswat].max()) if nc else 0
out_atoms = []
for k in np.nonzero(~iswat)[0]:
    x = xyz[k]; im = img[k]
    out_atoms.append(f'{newid[ids[k]]} {mol[k]} {typ[k]} {q[k]:.6f} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f} {im[0]} {im[1]} {im[2]}')
wb, wa = [], []
nb0, na0 = len(bonds_c), len(angles_c)
for i in range(len(wpos)):
    p = wpos[i] + lo                                           # the molecule whole; each atom wrapped into the box with its image flag
    im = np.floor((p - lo) / L).astype(int); pw = p - im * L
    a0 = nc + 3 * i + 1; m = molmax + i + 1
    for j, (t, qq) in enumerate(zip((tow, thw, thw), tmpl['q'] if tmpl else (0, 0, 0))):
        out_atoms.append(f'{a0 + j} {m} {t} {qq:.6f} {pw[j][0]:.6f} {pw[j][1]:.6f} {pw[j][2]:.6f} {im[j][0]} {im[j][1]} {im[j][2]}')
    wb += [(tmpl['bt'], a0, a0 + 1), (tmpl['bt'], a0, a0 + 2)]
    wa.append((tmpl['at'], a0 + 1, a0, a0 + 2))
counts.update({'atoms': len(out_atoms), 'bonds': nb0 + len(wb), 'angles': na0 + len(wa)})
if dihed: counts['dihedrals'] = len(dihed_c)
if impro: counts['impropers'] = len(impro_c)
with open(a.out, 'w') as f:
    f.write(f'rewater.py: {a.data} with m_w/m_c {a.mwmc} ({n_w} water, {method})\n\n')
    for k in ('atoms', 'atom types', 'bonds', 'bond types', 'angles', 'angle types', 'dihedrals', 'dihedral types', 'impropers', 'improper types'):
        if k in counts: f.write(f'{counts[k]} {k}\n')
    f.write('\n' + ''.join(f'{box[k][0]} {box[k][1]} {k}lo {k}hi\n' for k in 'xyz'))
    if tilt: f.write(f'{tilt[0]} {tilt[1]} {tilt[2]} xy xz yz\n')
    for name in order:
        if name == 'Velocities': continue
        f.write(f'\n{sec[name]["head"]}\n\n')
        if name == 'Atoms': f.write('\n'.join(out_atoms) + '\n')
        elif name == 'Bonds':
            f.writelines(f'{i + 1} {b[1]} {newid[b[2]]} {newid[b[3]]}\n' for i, b in enumerate(bonds_c))
            f.writelines(f'{nb0 + i + 1} {t} {i1} {i2}\n' for i, (t, i1, i2) in enumerate(wb))
        elif name == 'Angles':
            f.writelines(f'{i + 1} {an[1]} {newid[an[2]]} {newid[an[3]]} {newid[an[4]]}\n' for i, an in enumerate(angles_c))
            f.writelines(f'{na0 + i + 1} {t} {i1} {i2} {i3}\n' for i, (t, i1, i2, i3) in enumerate(wa))
        elif name in ('Dihedrals', 'Impropers'):
            cc = dihed_c if name == 'Dihedrals' else impro_c
            f.writelines(f'{i + 1} {c[1]} ' + ' '.join(str(newid[x]) for x in c[2:]) + '\n' for i, c in enumerate(cc))
        else: f.write('\n'.join(sec[name]['lines']) + '\n')
info = dict(source=os.path.abspath(a.data), mwmc=a.mwmc, n_water=n_w, n_water_before=n_old, method=method, tol=a.tol, seed=a.seed,
            atoms=len(out_atoms), box=[float(v) for v in L])
with open(os.path.splitext(a.out)[0] + '.json', 'w') as f: json.dump(info, f, indent=1)
print(f'wrote {a.out}: {len(out_atoms)} atoms ({nc} clay and Na, {n_w} water molecules, {method})')
