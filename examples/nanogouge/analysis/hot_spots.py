#!/usr/bin/env python3
"""Slip between clay particles during the shear (stage 04): are particle-particle interfaces "hot spots" (Zhou et al., Am. Mineral. 108, 87, 2023)?
Figures of the same kind as theirs, from the files that in.04_shear writes in each run folder:

  fig6_interface_slip.png    relative in-plane displacement s(t) of stacked particle pairs (their Fig. 6: stationary, stepwise, steady sliding)
                             + how much faster the water films shear than the gouge: (v / h) / gamma_dot of every interface
  fig7_velocity_profile.png  vx(z) of the water and of the clay across the gouge, and T(z) (their Fig. 7)
  fig8_tau_v.png             shear stress against slip velocity, log-log, ln tau = k ln v + ln tau_i (their Fig. 8; one point per run folder)
  fig9_water_structure.png   hydrogen bonds per water molecule and g_OO(r) of water between two particles / at one surface / in pores
                             (their Fig. 9; needs the water in the trajectory: env DUMP=all ./go.sh ...)
  fig10_viscosity.png        effective viscosity tau / shear rate against the shear rate, of the gouge and of the sliding films (their Fig. 10)
  hot_spots_summary.txt      the numbers

  python3 hot_spots.py RUN_N22                                   one run
  python3 hot_spots.py RUN_N22 RUN_N22_r3e8 RUN_N22_r1e9 RUN_N22_r3e9 --sigma 5.6    a series of shear rates (env FROM=RUN_N22 RATE=... ./go.sh)
Options: --gmin  strain from which the steady state is averaged (default: half of the strain reached);  --outdir  where the figures go
numpy, matplotlib; scipy for fig9.
"""
import argparse, glob, mmap, os, re, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

A2NM = 0.1
V_MS = 1e5          # 1 A/fs = 1e5 m/s


def run_params(run):
    """dt (fs) and the -var options of the folder (config.env lines written by go.sh)."""
    opts = {}
    f = os.path.join(run, 'config.env')
    if os.path.exists(f):
        for m in re.finditer(r'-var\s+(\S+)\s+(\S+?)["\s]', open(f).read()):
            opts[m.group(1)] = m.group(2)
    return float(opts.get('dt', 1.0)), opts


def read_shear(fn):
    d = np.loadtxt(fn, comments='#', ndmin=2)
    if d.shape[0] < 2:
        return None
    cols = ['step', 'gamma', 'tau', 'sn', 'mu', 'lz', 'rho', 'T', 'lx', 'ly']
    return {c: d[:, i] for i, c in enumerate(cols) if i < d.shape[1]}


def read_avetime(files):
    """fix ave/time mode vector files (one per job) -> steps (T,), data (T, rows, cols); a step written twice is kept once."""
    blocks = {}
    for fn in files:
        with open(fn) as f:
            lines = [l for l in f if not l.startswith('#')]
        i = 0
        while i < len(lines):
            h = lines[i].split()
            if len(h) < 2:
                i += 1; continue
            st, nr = int(h[0]), int(h[1])
            rows = np.array([l.split()[1:] for l in lines[i + 1:i + 1 + nr]], float)
            if rows.shape[0] == nr:
                blocks[st] = rows
            i += 1 + nr
    if not blocks:
        return None, None
    steps = np.array(sorted(blocks))
    return steps, np.stack([blocks[s] for s in steps])


def job_files(run, kind, sig):
    fs = glob.glob(os.path.join(run, f'{kind}_{sig}.*.dat'))
    return sorted(fs, key=lambda f: int(f.rsplit('.', 2)[-2]))


def tilt_of(sh, rate_dt):
    """cumulative tilt X(step) = integral of rate * lz d(step*dt) (the box flips do not matter for unwrapped positions)."""
    st, lz = sh['step'], sh['lz']
    st0 = np.concatenate([[0.0], st]); lz0 = np.concatenate([[lz[0]], lz])
    X = np.concatenate([[0.0], np.cumsum(0.5 * (lz0[1:] + lz0[:-1]) * np.diff(st0) * rate_dt)])
    return lambda s: np.interp(s, st0, X)


def particle_frame(gyr):
    """normal (eigenvector of the smallest eigenvalue), thickness and radius of each particle from its gyration tensor xx yy zz xy xz yz."""
    n = gyr.shape[0]
    G = np.zeros((n, 3, 3))
    G[:, 0, 0], G[:, 1, 1], G[:, 2, 2] = gyr[:, 0], gyr[:, 1], gyr[:, 2]
    G[:, 0, 1] = G[:, 1, 0] = gyr[:, 3]; G[:, 0, 2] = G[:, 2, 0] = gyr[:, 4]; G[:, 1, 2] = G[:, 2, 1] = gyr[:, 5]
    w, v = np.linalg.eigh(G)
    normal = v[:, :, 0]
    thick = np.sqrt(12 * np.clip(w[:, 0], 0, None))           # slab: Rg^2 = t^2 / 12
    radius = 2 * np.sqrt(np.clip(0.5 * (w[:, 1] + w[:, 2]), 0, None))   # disk: Rg^2 in plane = R^2 / 4
    return normal, thick, radius


def interfaces(run, sig, sh, dt, gmin, args):
    """stacked pairs at the first frame, their relative in-plane displacement in time and slip velocity in the steady state."""
    steps, P = read_avetime(job_files(run, 'particles', sig))
    if steps is None or len(steps) < 3:
        return None
    rate_dt = np.polyfit(sh['step'], sh['gamma'], 1)[0]            # gamma = rate * dt * step
    X = tilt_of(sh, rate_dt)
    lx = np.interp(steps, sh['step'], sh['lx']) if 'lx' in sh else None
    if lx is None:
        print(f'  {run}: shear_{sig}.dat has no Lx/Ly columns (older in.04_shear): no interface analysis'); return None
    ly = np.interp(steps, sh['step'], sh['ly']); lz = np.interp(steps, sh['step'], sh['lz'])
    com, gyr = P[:, :, 0:3], P[:, :, 6:12]
    n0, t0, R0 = particle_frame(gyr[0])
    H = lambda k: np.array([[lx[k], 0, X(steps[k])], [0, ly[k], 0], [0, 0, lz[k]]])
    H0inv = np.linalg.inv(H(0))
    N = com.shape[1]
    pairs = []
    for i in range(N):
        for j in range(i + 1, N):
            ci = abs(n0[i] @ n0[j])
            if ci < np.cos(np.radians(args.angle)):
                continue
            d = com[0, j] - com[0, i]
            img = np.round(H0inv @ d); d = d - H(0) @ img
            nb = n0[i] + np.sign(n0[i] @ n0[j]) * n0[j]; nb /= np.linalg.norm(nb)
            hc = abs(d @ nb); lat = np.linalg.norm(d - (d @ nb) * nb)
            tt = 0.5 * (t0[i] + t0[j])
            if hc < tt + args.hmax and lat < args.overlap * (R0[i] + R0[j]):
                pairs.append((i, j, img, nb, d, max(hc - tt, 2.0)))
    if not pairs:
        print(f'  {run} sigma {sig}: no stacked particle pairs (face to face, film < {args.hmax} A) at the first frame'); return None
    time_ns = steps * dt * 1e-6
    gam = np.interp(steps, sh['step'], sh['gamma'])
    win = gam >= gmin
    if win.sum() < 3:
        win = np.ones_like(win)
    out = []
    for (i, j, img, nb, d0, h) in pairs:
        d = np.array([com[k, j] - com[k, i] - H(k) @ img for k in range(len(steps))])
        dd = d - d0
        s = dd - np.outer(dd @ nb, nb)                           # in-plane relative displacement (A)
        # affine part: a homogeneous shear moves the two centres apart by (X(t)-X(0))/lz * dz in x
        aff = np.zeros_like(s); aff[:, 0] = (np.array([X(st) for st in steps]) - X(steps[0])) / lz * d0[2]
        aff = aff - np.outer(aff @ nb, nb)
        smag = np.linalg.norm(s, axis=1)
        v = np.polyfit(time_ns[win], smag[win], 1)[0] if win.sum() > 2 else np.nan   # A/ns
        out.append(dict(i=i, j=j, h=h, s=s, na=s - aff, smag=smag, v=v * 1e-6 * V_MS))     # v in m/s
    gdot = rate_dt / dt * 1e15                                   # 1/s
    return dict(t=time_ns, gam=gam, pairs=out, gdot=gdot, win=win)


def fig6(res, fn):
    runs = [(k, r) for k, r in res.items() if r.get('ifc')]
    if not runs:
        return
    fig, ax = plt.subplots(len(runs), 4, figsize=(17, 3.6 * len(runs)), squeeze=False)
    for row, (key, r) in enumerate(runs):
        I = r['ifc']; ps = sorted(I['pairs'], key=lambda p: -p['v'])
        pick = [('fastest', ps[0]), ('median', ps[len(ps) // 2]), ('slowest', ps[-1])]
        for c, (lab, p) in enumerate(pick):
            a = ax[row, c]
            a.plot(I['t'], p['s'][:, 0] * A2NM, label='x (flow)'); a.plot(I['t'], p['s'][:, 1] * A2NM, label='y')
            a.plot(I['t'], p['smag'] * A2NM, 'k', lw=0.8, label='|s|')
            a.set_title(f'{key}: {lab} interface ({p["i"]+1}-{p["j"]+1}), v = {p["v"]:.3g} m/s, film {p["h"]:.1f} A', fontsize=8)
            a.set_xlabel('t (ns)'); a.set_ylabel('in-plane slip s (nm)'); a.legend(fontsize=7)
        # hot spot factor: shear rate in the film over the shear rate of the gouge
        f = np.array([abs(p['v']) / (p['h'] * 1e-10) / I['gdot'] for p in ps])
        a = ax[row, 3]
        a.hist(np.log10(np.clip(f, 1e-3, None)), bins=20)
        a.axvline(0, color='k', ls='--', lw=0.8)
        a.set_xlabel('log10[(v/h) / gamma_dot]  (>0: the film shears faster than the gouge)'); a.set_ylabel('interfaces')
        a.set_title(f'{key}: {len(ps)} interfaces; {np.mean(f > 1) * 100:.0f} % faster than the gouge', fontsize=8)
    plt.tight_layout(); plt.savefig(fn, dpi=110); plt.close()


def profiles(run, sig, sh, gmin):
    steps, P = read_avetime(job_files(run, 'prof', sig))
    if steps is None:
        return None
    gam = np.interp(steps, sh['step'], sh['gamma']); m = gam >= gmin
    if m.sum() == 0:
        m[:] = True
    lz = np.interp(steps[m], sh['step'], sh['lz']).mean()
    D = P[m]                                                     # coord, n_water, vx_water, n_clay, vx_clay, T
    nw, nc = D[:, :, 1], D[:, :, 3]
    vw = (D[:, :, 2] * nw).sum(0) / np.maximum(nw.sum(0), 1e-9)
    vc = (D[:, :, 4] * nc).sum(0) / np.maximum(nc.sum(0), 1e-9)
    T = D[:, :, 5].mean(0)
    rate_dt = np.polyfit(sh['step'], sh['gamma'], 1)[0]
    return dict(z=D[0, :, 0] * lz, vw=vw * V_MS, vc=vc * V_MS, nw=nw.mean(0), nc=nc.mean(0), T=T, lz=lz)


def fig7(res, fn):
    runs = [(k, r) for k, r in res.items() if r.get('prof')]
    if not runs:
        return
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for key, r in runs:
        p = r['prof']
        l, = ax[0].plot(p['vw'], p['z'] * A2NM, 'o-', ms=3, label=f'{key} water')
        ax[0].plot(p['vc'], p['z'] * A2NM, 's--', ms=3, color=l.get_color(), label=f'{key} clay')
        ax[1].plot(p['nw'], p['z'] * A2NM, color=l.get_color(), label=f'{key} water'); ax[1].plot(p['nc'], p['z'] * A2NM, '--', color=l.get_color(), label=f'{key} clay')
        ax[2].plot(p['T'], p['z'] * A2NM, color=l.get_color(), label=key)
    ax[0].set_xlabel('vx (m/s)'); ax[0].set_ylabel('z (nm)'); ax[0].legend(fontsize=7); ax[0].set_title('velocity profile (steady state)', fontsize=9)
    ax[1].set_xlabel('atoms per layer'); ax[1].legend(fontsize=7); ax[1].set_title('where the water and the clay are', fontsize=9)
    ax[2].set_xlabel('T of the layer (K)'); ax[2].legend(fontsize=7); ax[2].set_title('temperature (flow of the layer removed)', fontsize=9)
    plt.tight_layout(); plt.savefig(fn, dpi=110); plt.close()


def fig8_10(res, fn8, fn10, summary):
    by_sig = {}
    for key, r in res.items():
        by_sig.setdefault(r['sig'], []).append((key, r))
    fig8, ax8 = plt.subplots(1, len(by_sig), figsize=(5.5 * len(by_sig), 4.2), squeeze=False)
    fig10, ax10 = plt.subplots(1, len(by_sig), figsize=(5.5 * len(by_sig), 4.2), squeeze=False)
    for c, (sig, lst) in enumerate(sorted(by_sig.items(), key=lambda kv: float(kv[0]))):
        V, tau, vh, gd, gdf = [], [], [], [], []
        for key, r in lst:
            V.append(r['V']); tau.append(r['tau']); gd.append(r['gdot'])
            I = r.get('ifc')
            if I:
                vs = np.sort([abs(p['v']) for p in I['pairs']])[::-1]
                hs = np.array([p['h'] for p in sorted(I['pairs'], key=lambda p: -abs(p['v']))])
                k = max(1, len(vs) // 10)                            # the fastest 10 % (at least one): the hot spots
                vh.append(vs[:k].mean()); gdf.append((vs[:k] / (hs[:k] * 1e-10)).mean())
            else:
                vh.append(np.nan); gdf.append(np.nan)
        V, tau, vh, gd, gdf = map(np.array, (V, tau, vh, gd, gdf))
        a = ax8[0, c]
        a.loglog(V, tau, 'o', label='gouge: V = gamma_dot Lz')
        ok = np.isfinite(vh) & (vh > 0)
        if ok.any():
            a.loglog(vh[ok], tau[ok], 's', label='hot interfaces (fastest 10 %)')
        for x, lab in ((V, 'gouge'), (vh, 'interfaces')):
            m = np.isfinite(x) & (x > 0) & (tau > 0)
            if m.sum() >= 2:
                k, lti = np.polyfit(np.log(x[m]), np.log(tau[m]), 1)
                xx = np.logspace(np.log10(x[m].min()), np.log10(x[m].max()), 20)
                a.loglog(xx, np.exp(lti) * xx ** k, '--', lw=0.8)
                summary.append(f'sigma {sig} {lab}: ln tau = {k:.3f} ln v + {lti:.3f}  (k = {k:.3f}, tau_i = {np.exp(lti):.3g} MPa at 1 m/s)')
        a.set_xlabel('v (m/s)'); a.set_ylabel('shear stress (MPa)'); a.set_title(f'sigma_n = {sig} MPa', fontsize=9); a.legend(fontsize=7)
        a = ax10[0, c]
        a.loglog(gd, tau * 1e6 / gd * 1e3, 'o', label='gouge: tau / gamma_dot')
        m = np.isfinite(gdf) & (gdf > 0)
        if m.any():
            a.loglog(gdf[m], tau[m] * 1e6 / gdf[m] * 1e3, 's', label='sliding films: tau / (v/h)')
        a.set_xlabel('shear rate (1/s)'); a.set_ylabel('effective viscosity (mPa s)'); a.set_title(f'sigma_n = {sig} MPa', fontsize=9); a.legend(fontsize=7)
    fig8.tight_layout(); fig8.savefig(fn8, dpi=110); fig10.tight_layout(); fig10.savefig(fn10, dpi=110)
    plt.close(fig8); plt.close(fig10)


# ---------------- fig 9: water structure from the trajectory ----------------
def type_names(run):
    names = {}
    f = os.path.join(run, 'gouge.data')
    if os.path.exists(f):
        with open(f) as fh:
            inm = False
            for l in fh:
                if l.strip() == 'Masses':
                    inm = True; continue
                if inm:
                    if l.strip() and l.strip()[0].isalpha():
                        break
                    p = l.split('#')
                    if len(p) == 2 and p[0].split():
                        names[int(p[0].split()[0])] = p[1].strip()
    return names


def dump_frames(fn):
    with open(fn, 'rb') as f:
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
        a = np.loadtxt([f.readline() for _ in range(n)], ndmin=2)
    if 'xy' in hb:
        xy, xz, yz = b[0][2], b[1][2], b[2][2]
        xlo = b[0][0] - min(0, xy, xz, xy + xz); xhi = b[0][1] - max(0, xy, xz, xy + xz)
        ylo = b[1][0] - min(0, yz); yhi = b[1][1] - max(0, yz); zlo, zhi = b[2][0], b[2][1]
    else:
        xy = xz = yz = 0.0; (xlo, xhi), (ylo, yhi), (zlo, zhi) = b[0][:2], b[1][:2], b[2][:2]
    H = np.array([[xhi - xlo, xy, xz], [0, yhi - ylo, yz], [0, 0, zhi - zlo]])
    c = {k: a[:, i] for i, k in enumerate(cols)}
    return step, np.array([xlo, ylo, zlo]), H, c


def images(pos, lo, H, rc):
    """positions in the box and their periodic images within rc of the faces: (points, index of the original atom)."""
    Hi = np.linalg.inv(H)
    s = (pos - lo) @ Hi.T; s -= np.floor(s)
    V = abs(np.linalg.det(H)); a, b, c = H[:, 0], H[:, 1], H[:, 2]
    w = np.array([V / np.linalg.norm(np.cross(b, c)), V / np.linalg.norm(np.cross(c, a)), V / np.linalg.norm(np.cross(a, b))])
    mrg = rc / w
    P, I = [], []
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                t = s + np.array([sx, sy, sz])
                m = np.all((t > -mrg) & (t < 1 + mrg), axis=1)
                P.append(lo + t[m] @ H.T); I.append(np.nonzero(m)[0])
    return lo + s @ H.T, np.concatenate(P), np.concatenate(I)


def water_structure(run, sig, args):
    from scipy.spatial import cKDTree
    fn = os.path.join(run, f'shear_{sig}.dump')
    if not os.path.exists(fn):
        return None
    names = type_names(run)
    tow = [t for t, n in names.items() if n.lower() == 'ow'] or [20]
    thw = [t for t, n in names.items() if n.lower() == 'hw'] or [19]
    tna = [t for t, n in names.items() if n.lower() == 'na']
    offs = dump_frames(fn)
    if not offs:
        return None
    sel = sorted(set([0, len(offs) - 1] + [k for k in args.frames if 0 <= k < len(offs)]))
    out = []
    for k in sel:
        step, lo, H, c = read_frame(fn, offs[k])
        typ = c['type'].astype(int); mol = c['mol'].astype(int)
        r = np.c_[c['xu'], c['yu'], c['zu']]
        isO = np.isin(typ, tow); isH = np.isin(typ, thw)
        if isO.sum() == 0:
            print(f'  {fn}: no water in the trajectory (env DUMP=all ./go.sh ... writes it): no fig9'); return None
        clay = ~(isO | isH | np.isin(typ, tna))
        # the two H of each water O (same molecule ID), as vectors from the O (unwrapped: the molecule is whole)
        iO = np.nonzero(isO)[0]; iH = np.nonzero(isH)[0]
        order = np.argsort(mol[iH], kind='stable'); iH = iH[order]
        mH = mol[iH]; first = np.searchsorted(mH, mol[iO]); okH = (first + 1 < len(mH))
        okH[okH] &= (mH[first[okH]] == mol[iO[okH]]) & (mH[first[okH] + 1] == mol[iO[okH]])
        iO = iO[okH]; first = first[okH]
        vOH = np.stack([r[iH[first]] - r[iO], r[iH[first + 1]] - r[iO]], 1)    # (nO, 2, 3)
        O, Oimg, Oidx = images(r[iO], lo, H, max(args.rmax, 3.5))
        tO = cKDTree(Oimg)
        # class of each water: clay atoms of how many particles within rc
        Cw, Cimg, Cidx = images(r[clay], lo, H, args.rc)
        cmol = mol[clay][Cidx]
        tC = cKDTree(Cimg)
        nb = tC.query_ball_point(O, args.rc)
        npart = np.array([len(set(cmol[l])) for l in nb])
        cls = np.where(npart >= 2, 0, np.where(npart == 1, 1, 2))          # 0 between particles, 1 at one surface, 2 pore
        # hydrogen bonds (Zhou et al.: O-O < 3.5 A, H...O < 2.45 A, angle H-Od...Oa < 30 deg)
        nhb = np.zeros(len(O))
        lists = tO.query_ball_point(O, 3.5)
        cos30 = np.cos(np.radians(30))
        for i, l in enumerate(lists):
            l = np.array(l)
            if len(l) == 0:
                continue
            d = Oimg[l] - O[i]; dn = np.linalg.norm(d, axis=1)
            keep = dn > 0.5; l, d, dn = l[keep], d[keep], dn[keep]
            for h in vOH[i]:
                hn = np.linalg.norm(h)
                hb = (np.linalg.norm(d - h, axis=1) < 2.45) & ((d @ h) / (dn * hn) > cos30)
                if hb.any():
                    nhb[i] += hb.sum(); np.add.at(nhb, Oidx[l[hb]], 1)
        # g_OO(r) of each class, normalised with the mean water O density of the box
        dr = 0.1; edges = np.arange(0, args.rmax + dr, dr); rr = 0.5 * (edges[1:] + edges[:-1])
        rho = len(O) / abs(np.linalg.det(H))
        g = {}
        far = tO.query_ball_point(O, args.rmax)
        for ccl in (0, 1, 2):
            idx = np.nonzero(cls == ccl)[0]
            if len(idx) == 0:
                continue
            h = np.zeros(len(rr))
            for i in idx:
                dist = np.linalg.norm(Oimg[far[i]] - O[i], axis=1)
                h += np.histogram(dist[dist > 0.5], edges)[0]
            g[ccl] = h / (len(idx) * 4 * np.pi * rr ** 2 * dr * rho)
        out.append(dict(step=step, frame=k, nhb={ccl: nhb[cls == ccl].mean() for ccl in (0, 1, 2) if (cls == ccl).any()},
                        frac={ccl: (cls == ccl).mean() for ccl in (0, 1, 2)}, r=rr, g=g))
    return out


CLS = {0: 'between two particles', 1: 'at one surface', 2: 'pore'}


def fig9(res, fn, summary):
    runs = [(k, r) for k, r in res.items() if r.get('water')]
    if not runs:
        return
    fig, ax = plt.subplots(len(runs), 2, figsize=(12, 3.8 * len(runs)), squeeze=False)
    for row, (key, r) in enumerate(runs):
        W = r['water']
        x = np.arange(3); wbar = 0.8 / len(W)
        for n, fr in enumerate(W):
            gam = r['gamma_of'](fr['step'])
            ax[row, 0].bar(x + n * wbar, [fr['nhb'].get(c, np.nan) for c in (0, 1, 2)], wbar, label=f'gamma {gam:.2f}')
            summary.append(f'{key} gamma {gam:.2f}: H-bonds per water ' + ', '.join(f'{CLS[c]} {fr["nhb"][c]:.2f} ({fr["frac"][c]*100:.0f} %)' for c in fr['nhb']))
        ax[row, 0].set_xticks(x + 0.4 - wbar / 2); ax[row, 0].set_xticklabels([CLS[c] for c in (0, 1, 2)], fontsize=8)
        ax[row, 0].set_ylabel('H-bonds per water molecule'); ax[row, 0].legend(fontsize=7); ax[row, 0].set_title(key, fontsize=9)
        for n, fr in enumerate(W):
            for c, g in fr['g'].items():
                ax[row, 1].plot(fr['r'], g, ls=['-', '--', ':'][n % 3], label=f'{CLS[c]}, gamma {r["gamma_of"](fr["step"]):.2f}')
        ax[row, 1].set_xlabel('r (A)'); ax[row, 1].set_ylabel('g_OO(r)'); ax[row, 1].legend(fontsize=6)
    plt.tight_layout(); plt.savefig(fn, dpi=110); plt.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--sigma', nargs='*', default=None, help='normal stresses to look at (default: every shear_<s>.dat)')
    ap.add_argument('--gmin', type=float, default=None, help='steady state from this strain (default: half of the strain reached)')
    ap.add_argument('--hmax', type=float, default=12.0, help='A: largest water film between two stacked particles')
    ap.add_argument('--angle', type=float, default=25.0, help='deg: largest angle between the normals of stacked particles')
    ap.add_argument('--overlap', type=float, default=0.75, help='lateral offset of stacked particles < overlap * (R1 + R2)')
    ap.add_argument('--rc', type=float, default=6.0, help='A: a water O within rc of clay atoms of 2 particles is "between two particles"')
    ap.add_argument('--rmax', type=float, default=8.0, help='A: range of g_OO(r)')
    ap.add_argument('--frames', type=int, nargs='*', default=[], help='more trajectory frames for fig9 (first and last always)')
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    res, summary = {}, []
    for run in args.runs:
        run = run.rstrip('/')
        dt, opts = run_params(run)
        sigs = args.sigma or sorted({re.match(r'shear_(.+)\.dat$', os.path.basename(f)).group(1) for f in glob.glob(os.path.join(run, 'shear_*.dat'))}, key=float)
        for sig in sigs:
            fn = os.path.join(run, f'shear_{sig}.dat')
            if not os.path.exists(fn):
                print(f'{fn}: not there'); continue
            sh = read_shear(fn)
            if sh is None:
                print(f'{fn}: too few lines yet'); continue
            gmin = args.gmin if args.gmin is not None else 0.5 * sh['gamma'].max()
            m = sh['gamma'] >= gmin
            rate_dt = np.polyfit(sh['step'], sh['gamma'], 1)[0]
            key = f'{os.path.basename(run)} {sig} MPa'
            r = dict(sig=sig, tau=sh['tau'][m].mean(), sn=sh['sn'][m].mean(), lz=sh['lz'][m].mean(), gdot=rate_dt / dt * 1e15)
            r['V'] = rate_dt / dt * r['lz'] * V_MS
            r['gamma_of'] = lambda st, sh=sh: float(np.interp(st, sh['step'], sh['gamma']))
            print(f'{key}: gamma up to {sh["gamma"].max():.3f}, steady state from {gmin:.3f}: tau {r["tau"]:.3g} MPa, sigma_n {r["sn"]:.3g} MPa, '
                  f'mu {r["tau"]/r["sn"]:.3f}, gamma_dot {r["gdot"]:.3g} /s, V {r["V"]:.3g} m/s')
            summary.append(f'{key}: tau {r["tau"]:.4g} MPa  sigma_n {r["sn"]:.4g} MPa  mu {r["tau"]/r["sn"]:.4f}  gamma_dot {r["gdot"]:.4g} /s  V {r["V"]:.4g} m/s  Lz {r["lz"]:.1f} A')
            r['ifc'] = interfaces(run, sig, sh, dt, gmin, args)
            if r['ifc']:
                v = np.array([abs(p['v']) for p in r['ifc']['pairs']])
                summary.append(f'  {len(v)} stacked interfaces: slip velocity median {np.median(v):.3g} m/s, fastest {v.max():.3g} m/s (gouge V {r["V"]:.3g} m/s)')
            r['prof'] = profiles(run, sig, sh, gmin)
            try:
                r['water'] = water_structure(run, sig, args)
            except ImportError:
                print('  scipy missing: no fig9'); r['water'] = None
            res[key] = r
    if not res:
        sys.exit('nothing to analyse yet')
    o = lambda n: os.path.join(args.outdir, n)
    fig6(res, o('fig6_interface_slip.png')); fig7(res, o('fig7_velocity_profile.png'))
    fig8_10(res, o('fig8_tau_v.png'), o('fig10_viscosity.png'), summary); fig9(res, o('fig9_water_structure.png'), summary)
    open(o('hot_spots_summary.txt'), 'w').write('\n'.join(summary) + '\n')
    print('\n'.join(summary)); print(f'figures in {os.path.abspath(args.outdir)}')


if __name__ == '__main__':
    main()
