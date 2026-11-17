#!/usr/bin/env python3
"""Figures of the kind of Lin & Hatano (2024, JGR Solid Earth 129, e2024JB028850, Figs. 5-10) for a series of water contents (porosities):
the run folder of in.water (go_water.sh) with water_contents.dat and one folder w<m_w/m_c>/ per water content.

  lin_fig5_rate.png          friction coefficient and porosity against the strain rate (the velocity steps), per water content       (their Fig. 5)
  lin_fig6_profiles_<s>.png  velocity profile vx / (rate Lz) across the gouge, one line per rate, and the angle between the particle normals
                             and z (dots: particles at the end of the first rate; line: layer means; grey: start of the shear)          (their Fig. 6)
  lin_fig7_stacking.png      parallel radial distribution function g||(r) of the particles (along each particle's normal, within its radius)
                             at the start and at the end of the shear, and the stacked neighbours per particle against the strain      (their Fig. 7)
  lin_fig8_healing.png       stacked neighbours per particle against the rest time in the hold (05, constant normal stress),
                             fit a + b log10(t); the slope b against the porosity of the hold                                            (their Fig. 8)
  lin_fig9_a_b.png           velocity dependence (a - b) = d(mu_ss) / ln(V2 / V1) between the velocity steps, against the porosity       (their Fig. 9)
  lin_fig10_band.png         R^2 of the velocity profile (1 = uniform shear) against the strain; the share of the height that shears
                             at least half as fast as the mean (100 % = uniform, small = a shear band), against the porosity            (their Fig. 10)
  lin_porosity.png           summary against the porosity during the shear: mu, (a - b), R^2, stacked neighbours
  lin_points.csv, lin_heal.csv   the numbers
Porosity phi = 1 - rho_b / rho_s (rho_b = dry bulk density written by the shear, rho_s = 2.7 g/cm3); m_w/m_c from water_contents.dat.
A stacked pair: normals within 25 deg, water film < 12 A between the faces, lateral offset < 0.75 (R1 + R2) (as hot_spots.py).
The steady state of a rate is the second half of its strain.  Few particles (N = 22) make g||(r), the orientation profile and the band noisy.
  python3 lin_figs.py RUN_N22_water [--outdir RUN_N22_water/figs]
numpy, matplotlib (hot_spots.py next to this file).
"""
import argparse, glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hot_spots as hs
from hot_spots import plt, style, head, layout, logaxes, dots, csv, SLOT, INK2, MUTED, LW, MS, PX, HEAD, SURF

RHO_S = 2.7
SEQ = ['#8fbcef', '#6aa6ea', '#4b92e3', '#2a78d6', '#2165b8', '#195299', '#123f78', '#0c2d58']   # one hue, wet (light) to dry (dark)


def seq(n):
    return [SEQ[int(round(i))] for i in (np.linspace(1, 7, n) if n > 1 else [3])]


def water_table(root):
    """the water contents of the series: folder, target and reached m_w/m_c (water_contents.dat; else the folder names)"""
    out = []
    f = os.path.join(root, 'water_contents.dat')
    if os.path.exists(f):
        for l in open(f):
            p = l.split()
            if len(p) >= 7 and not l.startswith('#'):
                out.append(dict(folder=os.path.join(root, p[6]), target=float(p[1]), mwmc=float(p[2])))
    else:
        for d in sorted(glob.glob(os.path.join(root, 'w*/'))):
            m = re.match(r'w([0-9.eE+-]+)$', os.path.basename(d.rstrip('/')))
            if m:
                out.append(dict(folder=d.rstrip('/'), target=float(m.group(1)), mwmc=float(m.group(1))))
    return sorted(out, key=lambda w: -w['mwmc'])


def min_image(d, H, Hinv):
    return d - (np.round(d @ Hinv.T) @ H.T)


def frame_geometry(com, gyr):
    n, t, R = hs.particle_frame(gyr)
    return n, t, R


def stacked_per_particle(com, gyr, H, angle=25.0, hmax=12.0, overlap=0.75):
    n, t, R = frame_geometry(com, gyr)
    N = len(com); Hinv = np.linalg.inv(H); cnt = 0
    for i in range(N - 1):
        j = np.arange(i + 1, N)
        c = n[j] @ n[i]
        par = np.abs(c) >= np.cos(np.radians(angle))
        if not par.any():
            continue
        j = j[par]; c = c[par]
        d = min_image(com[j] - com[i], H, Hinv)
        nb = n[i][None, :] + np.sign(c)[:, None] * n[j]; nb /= np.linalg.norm(nb, axis=1)[:, None]
        hc = np.abs((d * nb).sum(1)); lat = np.linalg.norm(d - (d * nb).sum(1)[:, None] * nb, axis=1)
        cnt += int(((hc < 0.5 * (t[i] + t[j]) + hmax) & (lat < overlap * (R[i] + R[j]))).sum())
    return 2.0 * cnt / N


def pardf(frames, rmax=60.0, h=1.0):
    """g||(r) averaged over frames (com, gyr, H): pairs along the normal of particle i within the radius of i; 1 for no correlation"""
    edges = np.arange(0, rmax + h, h); acc = np.zeros(len(edges) - 1); k = 0
    for com, gyr, H in frames:
        n, t, R = frame_geometry(com, gyr)
        N = len(com); Hinv = np.linalg.inv(H); rad = np.median(R); V = abs(np.linalg.det(H))
        cnt = np.zeros_like(acc)
        for i in range(N):
            d = min_image(np.delete(com, i, 0) - com[i], H, Hinv)
            rp = d @ n[i]; rq = np.linalg.norm(d - rp[:, None] * n[i], axis=1)
            cnt += np.histogram(np.abs(rp[rq < rad]), edges)[0]
        acc += cnt * V / (N * N * 2 * np.pi * rad ** 2 * h); k += 1
    return 0.5 * (edges[1:] + edges[:-1]), (acc / max(k, 1))


def load(w, sig, dt):
    """one water content and normal stress: shear, phases, particles, profiles, hold"""
    fn = os.path.join(w['folder'], f'shear_{sig}.dat')
    if not os.path.exists(fn):
        return None
    sh = hs.read_shear(fn)
    if sh is None or not np.isfinite(sh['lx']).any():
        return None
    segs = hs.segments(sh, dt)
    X = hs.tilt_of(sh)
    ok = np.isfinite(sh['lx'])
    box = lambda s: np.array([[np.interp(s, sh['step'][ok], sh['lx'][ok]), 0, X(s)], [0, np.interp(s, sh['step'][ok], sh['ly'][ok]), 0],
                              [0, 0, np.interp(s, sh['step'], sh['lz'])]])
    r = dict(w=w, sig=sig, sh=sh, segs=segs, box=box, dt=dt)
    st, P = hs.read_avetime(hs.job_files(w['folder'], 'particles', sig))
    if st is not None:
        r['pst'], r['com'], r['gyr'] = st, P[:, :, 0:3], P[:, :, 6:12]
    st, P = hs.read_avetime(hs.job_files(w['folder'], 'prof', sig))
    if st is not None:
        r['fst'], r['prof'] = st, P
    # the hold (05): particles, porosity; its box is the one at the end of the first rate (the start of the hold)
    st, P = hs.read_avetime(hs.job_files(w['folder'], 'healparticles', sig))
    hf = os.path.join(w['folder'], f'heal_{sig}.dat')
    if st is not None and os.path.exists(hf) and segs:
        hd = np.array([[float(x) for x in l.split()] for l in open(hf) if l.strip() and not l.startswith('#')])
        s_end0 = sh['step'][segs[0]['rows'][-1]]
        H0 = box(s_end0)
        lz = np.interp(st, hd[:, 0], hd[:, 3]) if len(hd) else np.full(len(st), H0[2, 2])
        r['heal'] = dict(steps=st, com=P[:, :, 0:3], gyr=P[:, :, 6:12], H=[np.array([[H0[0, 0], 0, H0[0, 2]], [0, H0[1, 1], 0], [0, 0, z]]) for z in lz],
                         phi=1 - np.nanmean(hd[len(hd) // 2:, 4]) / RHO_S if len(hd) else np.nan)
    return r


def phase_numbers(r):
    """per rate: steady state mu, tau, sigma_n, porosity, slip velocity V = rate Lz (m/s), R^2 and the shearing share at its end, stacked neighbours"""
    sh, out = r['sh'], []
    for k, sg in enumerate(r['segs']):
        win = sg['win'] if len(sg['win']) else sg['rows']
        d = dict(k=k, rate=sg['rate'], mu=np.nanmean(sh['mu'][win]), tau=np.nanmean(sh['tau'][win]), sn=np.nanmean(sh['sn'][win]),
                 phi=1 - np.nanmean(sh['rho'][win]) / RHO_S, V=sg['rate'] * np.nanmean(sh['lz'][win]) * 1e-10,
                 s0=sh['step'][sg['rows'][0]], s1=sh['step'][sg['rows'][-1]], g1=sh['gamma'][sg['rows'][-1]])
        d.update(r2=np.nan, share=np.nan, stack=np.nan)
        if 'prof' in r:
            m = np.where((r['fst'] > d['s0']) & (r['fst'] <= d['s1']))[0]
            if len(m):
                d['r2'], d['share'] = r2_share(r, m[-1])
        if 'pst' in r:
            m = np.where((r['pst'] > d['s0']) & (r['pst'] <= d['s1']))[0]
            if len(m):
                i = m[-1]; d['stack'] = stacked_per_particle(r['com'][i], r['gyr'][i], r['box'](r['pst'][i]))
        out.append(d)
    for a, b in zip(out[:-1], out[1:]):
        a['a_b'] = (b['mu'] - a['mu']) / np.log(b['V'] / a['V']) if a['V'] > 0 and b['V'] > 0 else np.nan
    if out:
        out[-1]['a_b'] = np.nan
    return out


def r2_share(r, i):
    """R^2 of the velocity profile of the profile block i against the uniform shear, and the share of the height shearing at >= half the mean rate"""
    D = r['prof'][i]; zr, nw, vw, nc, vc = D[:, 0], D[:, 1], D[:, 2], D[:, 3], D[:, 4]
    v = (nw * vw + nc * vc) / np.maximum(nw + nc, 1e-9)
    s = r['fst'][i]; sh = r['sh']
    rate = sh['rate'][np.argmin(np.abs(sh['step'] - s))] * 1e-15            # 1/fs
    Lz = np.interp(s, sh['step'], sh['lz']); V = rate * Lz                  # A/fs
    f = V * zr; f = f + np.mean(v - f)
    r2 = 1 - np.sum((v - f) ** 2) / max(np.sum((v - v.mean()) ** 2), 1e-30)
    loc = np.gradient(v, zr * Lz)
    return r2, 100.0 * np.mean(loc >= 0.5 * rate)


def r2_series(r):
    if 'prof' not in r:
        return np.array([]), np.array([])
    g = np.interp(r['fst'], r['sh']['step'], r['sh']['gamma'])
    return g, np.array([r2_share(r, i)[0] for i in range(len(r['fst']))])


def stack_series(r):
    if 'pst' not in r:
        return np.array([]), np.array([])
    g = np.interp(r['pst'], r['sh']['step'], r['sh']['gamma'])
    return g, np.array([stacked_per_particle(r['com'][i], r['gyr'][i], r['box'](r['pst'][i])) for i in range(len(r['pst']))])


def wlab(w):
    return f"m_w/m_c {w['mwmc']:.2g}"


# ---------------- figures ----------------
def fig5(R, sigs, ws, cols, fn):
    fig, ax = plt.subplots(2, len(sigs), figsize=(max(6.4, 4.4 * len(sigs) + 0.4), 6.2 + HEAD), squeeze=False)
    for c, sig in enumerate(sigs):
        for w, col in zip(ws, cols):
            r = R.get((w['folder'], sig))
            if not r or not r['ph']:
                continue
            x = np.array([p['rate'] for p in r['ph']])
            for row, key in ((0, 'mu'), (1, 'phi')):
                y = np.array([p[key] for p in r['ph']])
                ax[row, c].plot(x, y, color=col, lw=LW, zorder=2); dots(ax[row, c], x, y, col)
        for row, lab in ((0, 'friction coefficient μ = τ / σn'), (1, 'porosity φ during the shear')):
            a = ax[row, c]; a.set_xscale('log'); a.margins(x=0.15, y=0.15)
            style(a, f'σn = {sig} MPa'); a.set_xlabel('strain rate (1/s)'); a.set_ylabel(lab)
    head(fig, 'Friction and porosity against the strain rate (cf. Lin & Hatano 2024, Fig. 5)', [(wlab(w), col, 'line') for w, col in zip(ws, cols)])
    layout(fig); fig.savefig(fn); plt.close(fig)


def fig6(R, sig, ws, fn):
    ws = [w for w in ws if (w['folder'], sig) in R]
    if not ws:
        return
    fig, ax = plt.subplots(len(ws), 2, figsize=(7.2, 2.3 * len(ws) + HEAD), squeeze=False)
    rates = sorted({round(p['rate'], -3) for w in ws for p in R[(w['folder'], sig)]['ph']})
    rcol = {rt: SLOT[i % len(SLOT)] for i, rt in enumerate(rates)}
    for k, w in enumerate(ws):
        r = R[(w['folder'], sig)]
        a = ax[k, 0]
        a.plot([0, 1], [0, 1], color=MUTED, lw=PX * 2, ls='-', zorder=1)
        if 'prof' in r:
            for p in r['ph']:
                m = (r['fst'] > p['s0'] + 0.5 * (p['s1'] - p['s0'])) & (r['fst'] <= p['s1'])
                if not m.any():
                    continue
                D = r['prof'][m]; nw, nc = D[:, :, 1], D[:, :, 3]
                v = ((D[:, :, 2] * nw + D[:, :, 4] * nc).sum(0)) / np.maximum((nw + nc).sum(0), 1e-9)
                Lz = np.interp(p['s1'], r['sh']['step'], r['sh']['lz']); V = p['rate'] * 1e-15 * Lz
                a.plot(v / V - np.mean(v / V - D[0, :, 0]), D[0, :, 0], color=rcol[round(p['rate'], -3)], lw=LW)
        style(a, f"{wlab(w)} · φ {r['ph'][0]['phi']:.2f}" if r['ph'] else wlab(w)); a.set_xlabel('vx / (γ̇ Lz)'); a.set_ylabel('z / Lz')
        a = ax[k, 1]
        if 'pst' in r and r['ph']:
            p0 = r['ph'][0]
            ie = np.where(r['pst'] <= p0['s1'])[0]; ie = ie[-1] if len(ie) else len(r['pst']) - 1
            for i, col, lw in ((0, MUTED, LW), (ie, SLOT[0], LW)):
                n, _, _ = frame_geometry(r['com'][i], r['gyr'][i]); H = r['box'](r['pst'][i])
                ang = np.degrees(np.arccos(np.clip(np.abs(n[:, 2]), 0, 1)))
                zf = (r['com'][i] @ np.linalg.inv(H).T)[:, 2] % 1.0
                if i == ie:
                    a.plot(ang, zf, 'o', ls='none', color=col, ms=MS * 0.7, mec=SURF, mew=PX * 2, zorder=3)
                b = np.linspace(0, 1, 6); mid = 0.5 * (b[1:] + b[:-1])
                mean = [np.mean(ang[(zf >= b[q]) & (zf < b[q + 1])]) if ((zf >= b[q]) & (zf < b[q + 1])).any() else np.nan for q in range(5)]
                a.plot(mean, mid, color=col, lw=lw, zorder=2)
        a.set_xlim(0, 90)
        style(a, 'particle normal to z'); a.set_xlabel('angle (deg)'); a.set_ylabel('z / Lz (periodic)')
    head(fig, f'Velocity and orientation across the gouge, σn = {sig} MPa (cf. Lin & Hatano 2024, Fig. 6)',
         [(hs.rl(rt), rcol[rt], 'line') for rt in rates] + [('start', MUTED, 'line'), ('end of the first rate', SLOT[0], 'dot')])
    layout(fig); fig.savefig(fn); plt.close(fig)


def fig7(R, sigs, ws, cols, fn):
    nw = len(ws)
    fig = plt.figure(figsize=(max(7.2, 2.0 * nw + 0.6), 2.4 * len(sigs) + 2.8 + HEAD))
    gs = fig.add_gridspec(len(sigs) + 1, max(nw, 1))
    for c, sig in enumerate(sigs):
        for k, (w, col) in enumerate(zip(ws, cols)):
            a = fig.add_subplot(gs[c, k]); r = R.get((w['folder'], sig))
            if r and 'pst' in r:
                for i, cc in ((0, MUTED), (len(r['pst']) - 1, col)):
                    rr, g = pardf([(r['com'][i], r['gyr'][i], r['box'](r['pst'][i]))])
                    a.plot(rr / 10, g, color=cc, lw=LW)
                a.axhline(1, color=MUTED, lw=PX)
            style(a, f'{wlab(w)} · {sig} MPa'); a.set_xlabel('r|| (nm)')
            if k == 0:
                a.set_ylabel('g||(r)')
    a = fig.add_subplot(gs[len(sigs), :])
    for c, sig in enumerate(sigs):
        for w, col in zip(ws, cols):
            r = R.get((w['folder'], sig))
            if r:
                g, s = stack_series(r)
                if len(g):
                    a.plot(g, s, color=col, lw=LW, ls='-' if c == 0 else (0, (4, 2)))
    style(a, 'stacked neighbours per particle during the shear' + (f'  (solid {sigs[0]} MPa, dashed {sigs[1]} MPa)' if len(sigs) > 1 else ''))
    a.set_xlabel('shear strain γ'); a.set_ylabel('stacked neighbours')
    head(fig, 'Stacking of the particles during the shear (cf. Lin & Hatano 2024, Fig. 7)',
         [('start of the shear', MUTED, 'line')] + [(wlab(w), col, 'line') for w, col in zip(ws, cols)])
    layout(fig); fig.savefig(fn); plt.close(fig)


def fig8(R, sigs, ws, cols, fn, rows):
    fig, ax = plt.subplots(2, len(sigs), figsize=(max(6.4, 4.4 * len(sigs) + 0.4), 6.0 + HEAD), squeeze=False)
    for c, sig in enumerate(sigs):
        pb = []
        for w, col in zip(ws, cols):
            r = R.get((w['folder'], sig))
            if not r or 'heal' not in r:
                continue
            h = r['heal']; t = h['steps'] * r['dt'] * 1e-3                        # ps
            s = np.array([stacked_per_particle(h['com'][i], h['gyr'][i], h['H'][i]) for i in range(len(t))])
            m = t > 0
            for tt, ss in zip(t, s):
                rows.append([w['folder'], w['mwmc'], sig, tt, ss])
            ax[0, c].plot(t[m], s[m], color=col, lw=LW)
            if m.sum() >= 3:
                b, a0 = np.polyfit(np.log10(t[m]), s[m], 1); pb.append((h['phi'], b, col))
        a = ax[0, c]; a.set_xscale('log'); style(a, f'σn = {sig} MPa · constant normal stress'); a.set_xlabel('rest time (ps)'); a.set_ylabel('stacked neighbours per particle')
        a = ax[1, c]
        if pb:
            pb.sort()
            a.plot([p[0] for p in pb], [p[1] for p in pb], color=MUTED, lw=LW, zorder=2)
            for x, y, col in pb:
                dots(a, [x], [y], col)
        a.axhline(0, color=MUTED, lw=PX); a.margins(x=0.2, y=0.2)
        style(a, 'restoration rate b of a + b log10(t)'); a.set_xlabel('porosity φ of the hold'); a.set_ylabel('b (neighbours per decade)')
    head(fig, 'Stacking during the hold after the shear (cf. Lin & Hatano 2024, Fig. 8)', [(wlab(w), col, 'line') for w, col in zip(ws, cols)])
    layout(fig); fig.savefig(fn); plt.close(fig)


def by_rate_pairs(R, sig, ws):
    pairs = {}
    for w in ws:
        r = R.get((w['folder'], sig))
        if not r:
            continue
        for p, q in zip(r['ph'][:-1], r['ph'][1:]):
            key = (round(p['rate'], -3), round(q['rate'], -3))
            pairs.setdefault(key, []).append((0.5 * (p['phi'] + q['phi']), p['a_b']))
    return pairs


def fig9(R, sigs, ws, fn):
    fig, ax = plt.subplots(1, len(sigs), figsize=(max(6.4, 4.4 * len(sigs) + 0.4), 3.6 + HEAD), squeeze=False)
    keys = sorted({k for sig in sigs for k in by_rate_pairs(R, sig, ws)})
    kcol = {k: SLOT[i % len(SLOT)] for i, k in enumerate(keys)}
    for c, sig in enumerate(sigs):
        a = ax[0, c]
        for k, pts in sorted(by_rate_pairs(R, sig, ws).items()):
            pts = sorted(p for p in pts if np.isfinite(p[1]))
            if pts:
                a.plot([p[0] for p in pts], [p[1] for p in pts], color=kcol[k], lw=LW, zorder=2); dots(a, [p[0] for p in pts], [p[1] for p in pts], kcol[k])
        a.axhline(0, color=MUTED, lw=PX); a.margins(x=0.2, y=0.2)
        style(a, f'σn = {sig} MPa'); a.set_xlabel('porosity φ during the shear'); a.set_ylabel('(a − b) = Δμss / ln(V2/V1)')
    head(fig, 'Velocity dependence of friction against the porosity (cf. Lin & Hatano 2024, Fig. 9)',
         [(f'{hs.rl(k[0]).replace(" /s", "")} → {hs.rl(k[1])}', kcol[k], 'line') for k in keys])
    layout(fig); fig.savefig(fn); plt.close(fig)


def fig10(R, sigs, ws, cols, fn):
    fig, ax = plt.subplots(2, len(sigs), figsize=(max(6.4, 4.4 * len(sigs) + 0.4), 6.2 + HEAD), squeeze=False)
    rates = sorted({round(p['rate'], -3) for r in R.values() for p in r['ph']})
    rcol = {rt: SLOT[i % len(SLOT)] for i, rt in enumerate(rates)}
    for c, sig in enumerate(sigs):
        a = ax[0, c]
        for w, col in zip(ws, cols):
            r = R.get((w['folder'], sig))
            if r:
                g, r2 = r2_series(r)
                if len(g):
                    a.plot(g, r2, color=col, lw=LW)
        a.axhline(0.9, color=MUTED, lw=PX)
        style(a, f'σn = {sig} MPa · R² of the velocity profile'); a.set_xlabel('shear strain γ'); a.set_ylabel('R² (1 = uniform shear)')
        a = ax[1, c]
        for rt in rates:
            pts = sorted((p['phi'], p['share']) for w in ws for p in R.get((w['folder'], sig), {'ph': []})['ph'] if round(p['rate'], -3) == rt and np.isfinite(p['share']))
            if pts:
                a.plot([p[0] for p in pts], [p[1] for p in pts], color=rcol[rt], lw=LW, zorder=2); dots(a, [p[0] for p in pts], [p[1] for p in pts], rcol[rt])
        a.set_ylim(0, 105); a.margins(x=0.2)
        style(a, 'height shearing at ≥ half the mean rate (end of each rate)'); a.set_xlabel('porosity φ during the shear'); a.set_ylabel('share of the height (%)')
    head(fig, 'Uniform shear or a shear band (cf. Lin & Hatano 2024, Fig. 10)',
         [(wlab(w), col, 'line') for w, col in zip(ws, cols)] + [(hs.rl(rt), rcol[rt], 'dot') for rt in rates])
    layout(fig); fig.savefig(fn); plt.close(fig)


def fig_summary(R, sigs, ws, fn):
    fig, ax = plt.subplots(2, 2, figsize=(7.6, 6.2 + HEAD), squeeze=False)
    scol = {sig: SLOT[i % len(SLOT)] for i, sig in enumerate(sigs)}
    panels = ((0, 0, 'mu', 'friction coefficient μ (first rate)'), (0, 1, 'a_b', '(a − b) (first velocity step)'),
              (1, 0, 'r2', 'R² of the velocity profile (end of the first rate)'), (1, 1, 'stack', 'stacked neighbours per particle (end of the first rate)'))
    for i, j, key, lab in panels:
        a = ax[i, j]
        for sig in sigs:
            pts = sorted((R[(w['folder'], sig)]['ph'][0]['phi'], R[(w['folder'], sig)]['ph'][0].get(key, np.nan)) for w in ws
                         if (w['folder'], sig) in R and R[(w['folder'], sig)]['ph'])
            pts = [p for p in pts if np.isfinite(p[1])]
            if pts:
                a.plot([p[0] for p in pts], [p[1] for p in pts], color=scol[sig], lw=LW, zorder=2); dots(a, [p[0] for p in pts], [p[1] for p in pts], scol[sig])
        a.margins(x=0.2, y=0.2); style(a, lab); a.set_xlabel('porosity φ during the shear'); a.set_ylabel(lab.split(' (')[0])
    head(fig, 'Shear behaviour against the porosity (water content) of the gouge', [(f'σn = {sig} MPa', scol[sig], 'line') for sig in sigs])
    layout(fig); fig.savefig(fn); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root', help='run folder of in.water (water_contents.dat, w*/)')
    ap.add_argument('--outdir', default=None, help='default: <root>/figs')
    args = ap.parse_args()
    root = args.root.rstrip('/'); out = args.outdir or os.path.join(root, 'figs'); os.makedirs(out, exist_ok=True)
    o = lambda n: os.path.join(out, n)
    dt, _ = hs.run_params(root)
    ws = water_table(root)
    if not ws:
        print(f'{root}: no water contents yet (water_contents.dat, w*/)'); return
    sigs = sorted({re.match(r'shear_(.+)\.dat$', os.path.basename(f)).group(1) for w in ws for f in glob.glob(os.path.join(w['folder'], 'shear_*.dat'))}, key=float)
    R = {}
    for w in ws:
        for sig in sigs:
            r = load(w, sig, dt)
            if r:
                r['ph'] = phase_numbers(r); R[(w['folder'], sig)] = r
    ws = [w for w in ws if any((w['folder'], s) in R for s in sigs)]
    if not R:
        print(f'{root}: no shear yet'); return
    cols = seq(len(ws))
    rows = [[os.path.basename(r['w']['folder']), r['w']['mwmc'], r['sig'], p['rate'], p['phi'], p['tau'], p['sn'], p['mu'], p['V'], p['a_b'], p['r2'], p['share'], p['stack'], p['g1']]
            for r in R.values() for p in r['ph']]
    csv(o('lin_points.csv'), ['folder', 'mw_mc', 'sigma_MPa', 'rate_1_s', 'porosity_shear', 'tau_MPa', 'sigma_n_MPa', 'mu', 'V_m_s', 'a_minus_b_to_next_rate',
                               'R2_end', 'shearing_share_pct_end', 'stacked_per_particle_end', 'strain_end'], rows)
    hrows = []
    for name, draw in (('fig5', lambda: fig5(R, sigs, ws, cols, o('lin_fig5_rate.png'))),
                    ('fig6', lambda: [fig6(R, s, ws, o(f'lin_fig6_profiles_{s}.png')) for s in sigs]),
                    ('fig7', lambda: fig7(R, sigs, ws, cols, o('lin_fig7_stacking.png'))),
                    ('fig8', lambda: fig8(R, sigs, ws, cols, o('lin_fig8_healing.png'), hrows)),
                    ('fig9', lambda: fig9(R, sigs, ws, o('lin_fig9_a_b.png'))),
                    ('fig10', lambda: fig10(R, sigs, ws, cols, o('lin_fig10_band.png'))),
                    ('summary', lambda: fig_summary(R, sigs, ws, o('lin_porosity.png')))):
        try:
            draw()
        except Exception as e:                       # one figure without its data must not stop the others
            print(f'{name}: not drawn ({type(e).__name__}: {e})')
    if hrows:
        csv(o('lin_heal.csv'), ['folder', 'mw_mc', 'sigma_MPa', 'rest_time_ps', 'stacked_per_particle'], hrows)
    for r in R.values():
        p = r['ph'][0]
        print(f"{os.path.basename(r['w']['folder'])} (m_w/m_c {r['w']['mwmc']:.3g}) {r['sig']} MPa: porosity {p['phi']:.3f}  mu {p['mu']:.4f}  (a-b) {p['a_b']:.4g}  R2 {p['r2']:.3f}  stacked {p['stack']:.2f}")
    print(f'figures and numbers in {out}')


if __name__ == '__main__':
    main()
