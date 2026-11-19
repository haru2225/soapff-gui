#!/usr/bin/env python3
"""Are the clay particles stacked?  From the trajectories of a run folder (clay atoms: id mol type xu yu zu), at the end of 01, along the
dehydration 02, at the end of each consolidation 03, and at the start and end of each shear 04:

  stacked pairs     two particles face to face: normals within --angle, water film < --hmax A, lateral offset < --overlap (R1 + R2)
                    (the same test as hot_spots.py for the interfaces of fig6)
  in stacks         particles that belong to a stack (a chain of stacked pairs);  largest = particles in the largest stack
  film              mean water film between the stacked particles (A);  S = orientational order of the normals (0 random, 1 aligned)

  foliation (04)    along each shear (--frames04 frames): S, the tilt of the mean particle plane from the shear plane (deg, in the
                    plane of the shear xz; + towards +x), the share of particles within 20 deg of the shear plane, stacks; against the strain
  python3 stacks.py RUN_N22_one [--outdir RUN_N22_one/figs]     -> table on the screen, stacks.csv, stacks.png, foliation.png, foliation.csv
Needs numpy only (matplotlib for the figure).
"""
import argparse, glob, mmap, os, re, sys
import numpy as np


def frame_offsets(fn):
    with open(fn, 'rb') as f:
        if os.fstat(f.fileno()).st_size == 0:
            return []
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        offs = [m.start() for m in re.finditer(rb'ITEM: TIMESTEP', mm)]
        mm.close()
    return offs


def read_frame(fn, off):
    with open(fn) as f:
        f.seek(off)
        f.readline(); step = int(f.readline())
        f.readline(); n = int(f.readline())
        hb = f.readline()
        b = [list(map(float, f.readline().split())) for _ in range(3)]
        cols = f.readline().split()[2:]
        txt = ''.join(f.readline() for _ in range(n))
    a = np.array(txt.split(), float)
    if a.size != n * len(cols):
        raise ValueError('frame not complete')
    a = a.reshape(n, len(cols))
    if 'xy' in hb:
        xy, xz, yz = b[0][2], b[1][2], b[2][2]
        xlo = b[0][0] - min(0, xy, xz, xy + xz); xhi = b[0][1] - max(0, xy, xz, xy + xz)
        ylo = b[1][0] - min(0, yz); yhi = b[1][1] - max(0, yz); zlo, zhi = b[2][0], b[2][1]
    else:
        xy = xz = yz = 0.0; (xlo, xhi), (ylo, yhi), (zlo, zhi) = b[0][:2], b[1][:2], b[2][:2]
    H = np.array([[xhi - xlo, xy, xz], [0, yhi - ylo, yz], [0, 0, zhi - zlo]])
    return step, H, {k: a[:, i] for i, k in enumerate(cols)}


def na_types(run):
    """atom types named Na in the Masses section of gouge.data (the counter-ions: not part of the particles)"""
    t = set(); f = os.path.join(run, 'gouge.data')
    if os.path.exists(f):
        inm = False
        with open(f) as fh:
            for l in fh:
                s = l.strip()
                if s.startswith('Masses'): inm = True; continue
                if inm and s and not s[0].isdigit(): break
                if inm and '#' in s and s.split('#')[1].strip().lower() == 'na': t.add(int(s.split()[0]))
    return t


def analyse(fn, off, tna, args):
    step, H, c = read_frame(fn, off)
    mol = c['mol'].astype(int); typ = c['type'].astype(int)
    r = np.c_[c['xu'], c['yu'], c['zu']]
    keep = ~np.isin(typ, list(tna)) if tna else np.ones(len(typ), bool)
    mols = np.unique(mol[keep])
    # the particles: the molecules with many atoms (water, when in the trajectory, has 3)
    com, nrm, thk, rad = [], [], [], []
    for m in mols:
        x = r[keep & (mol == m)]
        if len(x) < 50:
            continue
        cc = x.mean(0); w, v = np.linalg.eigh(np.cov((x - cc).T, bias=True))
        com.append(cc); nrm.append(v[:, 0]); thk.append(np.sqrt(12 * max(w[0], 0))); rad.append(2 * np.sqrt(max(0.5 * (w[1] + w[2]), 0)))
    com, nrm, thk, rad = map(np.array, (com, nrm, thk, rad))
    N = len(com)
    Hi = np.linalg.inv(H)
    pairs, films = [], []
    for i in range(N):
        for j in range(i + 1, N):
            if abs(nrm[i] @ nrm[j]) < np.cos(np.radians(args.angle)):
                continue
            d = com[j] - com[i]; d = d - H @ np.round(Hi @ d)
            nb = nrm[i] + np.sign(nrm[i] @ nrm[j]) * nrm[j]; nb /= np.linalg.norm(nb)
            hc = abs(d @ nb); lat = np.linalg.norm(d - (d @ nb) * nb); tt = 0.5 * (thk[i] + thk[j])
            if hc < tt + args.hmax and lat < args.overlap * (rad[i] + rad[j]):
                pairs.append((i, j)); films.append(max(hc - tt, 0.0))
    # stacks = connected groups of stacked pairs
    parent = list(range(N))
    def root(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]; k = parent[k]
        return k
    for i, j in pairs:
        parent[root(i)] = root(j)
    sizes = {}
    for k in range(N):
        sizes[root(k)] = sizes.get(root(k), 0) + 1
    instack = sum(s for s in sizes.values() if s > 1)
    Q = sum(1.5 * np.outer(n, n) - 0.5 * np.eye(3) for n in nrm) / max(N, 1)
    if N:
        w, v = np.linalg.eigh(Q); S = float(w[-1]); dvec = v[:, -1] * (1 if v[2, -1] >= 0 else -1)
        tilt = float(np.degrees(np.arctan2(dvec[0], dvec[2])))          # mean normal from z in the xz plane = plane tilt from the shear plane
        flat = float(np.mean(np.abs(nrm[:, 2]) > np.cos(np.radians(20))))
    else:
        S = tilt = flat = np.nan
    return dict(step=step, N=N, pairs=len(pairs), instack=instack, largest=max(sizes.values()) if sizes else 0,
                film=float(np.mean(films)) if films else np.nan, S=S, tilt=tilt, flat=flat)


def pick(n, k):
    """k frame indices spread over n frames, always the last one"""
    if n <= k:
        return list(range(n))
    return sorted(set(np.linspace(0, n - 1, k).round().astype(int).tolist()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('--hmax', type=float, default=12.0, help='A: largest water film between two stacked particles')
    ap.add_argument('--angle', type=float, default=25.0, help='deg: largest angle between the normals of stacked particles')
    ap.add_argument('--overlap', type=float, default=0.75, help='lateral offset of stacked particles < overlap * (R1 + R2)')
    ap.add_argument('--frames02', type=int, default=11, help='frames of the dehydration (02) to look at')
    ap.add_argument('--frames04', type=int, default=16, help='frames of each shear (04) to look at')
    ap.add_argument('--outdir', default=None, help='default: <run>/figs')
    args = ap.parse_args()
    run = args.run.rstrip('/'); out = args.outdir or os.path.join(run, 'figs'); os.makedirs(out, exist_ok=True)
    tna = na_types(run)
    todo = []                                                   # (stage label, file, frame indices)
    f = os.path.join(run, '01_equil.dump')
    if os.path.exists(f): todo.append(('01 end', f, [-1]))
    f = os.path.join(run, '02_dehydrate.dump')
    if os.path.exists(f): todo.append(('02', f, None))
    for f in sorted(glob.glob(os.path.join(run, '03_consolidate_*.dump')), key=lambda x: float(re.findall(r'_([0-9.]+)\.dump$', x)[0])):
        todo.append((f'03 {re.findall(r"_([0-9.]+)[.]dump$", f)[0]} MPa end', f, [-1]))
    for f in sorted(glob.glob(os.path.join(run, 'shear_*.dump')), key=lambda x: x):
        m = re.match(r'shear_([0-9.]+)\.dump$', os.path.basename(f))
        if m: todo.append((f'04 {m.group(1)} MPa', f, 'shear'))
    if not todo:
        sys.exit(f'{run}: no trajectory (01_equil.dump, 02_dehydrate.dump, 03_consolidate_*.dump, shear_*.dump) yet')
    rows = []
    print(f'{"stage":<18} {"step":>9} {"stacked pairs":>13} {"in stacks":>10} {"largest":>8} {"film (A)":>9} {"S":>6}'
          '   (tilt: deg of the mean particle plane from the shear plane; flat: share within 20 deg of it)')
    for lab, f, idx in todo:
        offs = frame_offsets(f)
        if not offs:
            continue
        ks = pick(len(offs), args.frames02) if idx is None else pick(len(offs), args.frames04) if idx == 'shear' else sorted(set(i % len(offs) for i in idx))
        gam = None
        if idx == 'shear':                                       # strain of each frame from shear_<s>.dat
            fd = f[:-5] + '.dat'
            try:
                D = np.atleast_2d(np.array([[float(x) for x in l.split()[:2]] for l in open(fd) if l.strip() and not l.startswith('#')]))
                if D.size == 0: raise ValueError('no lines')
                gam = lambda st, D=D: float(np.interp(st, np.r_[0, D[:, 0]], np.r_[0, D[:, 1]]))
            except Exception:
                gam = None
        for k in ks:
            try:
                d = analyse(f, offs[k], tna, args)
            except (ValueError, IndexError) as e:
                print(f'  {os.path.basename(f)} frame {k}: not readable ({e})'); continue
            lb = lab if idx is not None or len(ks) == 1 else f'{lab} {k + 1}/{len(offs)}'
            if idx == 'shear':
                d['gamma'] = gam(d['step']) if gam else np.nan; d['sig'] = lab.split()[1]
                lb = f'{lab} g {d["gamma"]:.2f}'
            rows.append((lb, d))
            print(f'{lb:<18} {d["step"]:>9} {d["pairs"]:>13} {d["instack"]:>6} / {d["N"]:<3} {d["largest"]:>6} {d["film"]:>9.1f} {d["S"]:>6.3f}'
                  f'   tilt {d["tilt"]:6.1f}  flat {d["flat"]:.2f}', flush=True)
    with open(os.path.join(out, 'stacks.csv'), 'w') as fh:
        fh.write('stage,step,particles,stacked_pairs,particles_in_stacks,largest_stack,mean_film_A,S,tilt_deg,flat_share,gamma\n')
        for lb, d in rows:
            fh.write(f'{lb},{d["step"]},{d["N"]},{d["pairs"]},{d["instack"]},{d["largest"]},{d["film"]:.2f},{d["S"]:.4f},{d["tilt"]:.2f},{d["flat"]:.3f},{d.get("gamma", np.nan):.4f}\n')
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import hot_spots as hs
        plt = hs.plt
    except Exception:
        print(f'numbers in {os.path.join(out, "stacks.csv")} (no matplotlib: no figure)'); return
    x = np.arange(len(rows))
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, key, yl in zip(axs, ('pairs', 'instack', 'S'), ('stacked pairs', 'particles in stacks', 'order parameter S')):
        y = [d[key] for _, d in rows]
        ax.plot(x, y, color=hs.SLOT[0], lw=hs.LW, marker='o', ms=hs.MS)
        ax.set_xticks(x); ax.set_xticklabels([lb for lb, _ in rows], rotation=60, ha='right', fontsize=8)
        ax.set_ylabel(yl); hs.style(ax, yl)
        if key == 'instack' and rows: ax.axhline(rows[0][1]['N'], color=hs.MUTED, lw=1)
    fig.suptitle(f'Stacking of the clay particles along the run · {os.path.basename(os.path.abspath(run))}', x=0.01, ha='left', fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'stacks.png'), dpi=130); plt.close(fig)
    print(f'stacks.png and stacks.csv in {os.path.abspath(out)}')
    # the foliation forming in the shear: against the strain, one line per normal stress
    sh = [d for _, d in rows if 'gamma' in d and np.isfinite(d['gamma'])]
    if not sh:
        return
    sigs = sorted({d['sig'] for d in sh}, key=float)
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.4))
    keys = (('S', 'order parameter S'), ('tilt', 'tilt of the particles from the shear plane (deg)'),
            ('flat', 'share of particles within 20 deg of the shear plane'), ('instack', 'particles in stacks'))
    for ax, (key, yl) in zip(axs, keys):
        for k, sg in enumerate(sigs):
            dd = sorted([d for d in sh if d['sig'] == sg], key=lambda d: d['gamma'])
            ax.plot([d['gamma'] for d in dd], [d[key] for d in dd], color=hs.SLOT[k % len(hs.SLOT)], lw=hs.LW, marker='o', ms=hs.MS, label=f'{sg} MPa')
        if key == 'tilt': ax.axhline(0, color=hs.MUTED, lw=1)
        ax.set_xlabel('shear strain'); ax.set_ylabel(yl); hs.style(ax, yl)
    axs[0].legend(frameon=False)
    fig.suptitle(f'Foliation forming in the shear · {os.path.basename(os.path.abspath(run))}', x=0.01, ha='left', fontweight='bold')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'foliation.png'), dpi=130); plt.close(fig)
    print(f'foliation.png in {os.path.abspath(out)}')


if __name__ == '__main__':
    main()
