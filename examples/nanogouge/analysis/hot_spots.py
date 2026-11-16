#!/usr/bin/env python3
"""Slip between clay particles during the shear (stage 04): are particle-particle interfaces "hot spots" (Zhou et al., Am. Mineral. 108, 87, 2023)?
Figures of the same kind as theirs, from the files that in.04_shear writes in each run folder:

  fig6_interface_slip.png    relative in-plane displacement s(t) of stacked particle pairs (their Fig. 6: stationary, stepwise, steady sliding)
                             + how much faster the water films shear than the gouge: (v / h) / gamma_dot of every interface
  fig7_velocity_profile.png  vx(z) / (gamma_dot Lz) of the water and of the clay across the gouge (diagonal = uniform shear; a step = a shear band),
                             water share and T(z) (their Fig. 7)
  fig8_tau_v.png             shear stress against slip velocity, log-log, ln tau = k ln v + ln tau_i (their Fig. 8; one point per run folder)
  fig9_water_structure.png   hydrogen bonds per water molecule and g_OO(r) of water between two particles / at one surface / in pores
                             (their Fig. 9; needs the water in the trajectory: env DUMP=all ./go.sh ...)
  fig10_viscosity.png        effective viscosity tau / shear rate against the shear rate, of the gouge and of the sliding films (their Fig. 10)
  fig*.csv                   the numbers behind each figure (to plot them again elsewhere);  hot_spots_summary.txt  the main numbers

  One run is enough: in.04_shear shears at 1e8 /s, then steps the rate up to 1e10 /s (velocity steps: one point per rate in fig8, fig10), and writes
  the whole system with the water at the start and at the end of each phase (shear_<s>_all_<step>.dump, for fig9).  The jobs run this script
  after each job of the shear (RUN/figs); the analysed trajectory frames are kept in <outdir>/.cache.

  python3 hot_spots.py RUN_N67 --outdir figs                     one run
  python3 hot_spots.py RUN_N22 RUN_N22_r3e8 RUN_N22_r1e9 RUN_N22_r3e9 --sigma 5.6    a series of shear rates (env FROM=RUN_N22 RATE=... ./go.sh)
Options: --gmin  strain from which the steady state is averaged (default: half of the strain reached);  --outdir  where the figures go
numpy, matplotlib; scipy for fig9.
"""
import argparse, glob, json, mmap, os, re, sys, time
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter

A2NM = 0.1
V_MS = 1e5          # 1 A/fs = 1e5 m/s

# ---------------- look of the figures: one palette, thin marks, recessive grid ----------------
DPI = 150
PX = 72.0 / DPI                                       # one pixel in points
SURF, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
# categorical colours in a fixed order (checked for colour-blind readers: the first three are distinct in every pair); never cycled
SLOT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
C_GOUGE, C_IFC = SLOT[0], SLOT[1]                     # the same meaning in fig8 and fig10
LW, MS, RING = 2 * PX, 10 * PX, 2 * PX
plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF, 'figure.dpi': DPI, 'savefig.dpi': DPI,
                     'font.family': 'sans-serif', 'font.size': 8, 'axes.edgecolor': AXIS, 'axes.linewidth': PX, 'axes.labelcolor': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.axisbelow': True,
                     'lines.linewidth': LW, 'lines.solid_capstyle': 'round', 'lines.solid_joinstyle': 'round',
                     'legend.frameon': False, 'xtick.color': AXIS, 'ytick.color': AXIS})


def style(ax, title=None):
    ax.tick_params(which='both', colors=MUTED, labelsize=7, width=PX)
    ax.tick_params(which='minor', length=0)
    ax.grid(True, which='major', color=GRID, lw=PX, ls='-')
    ax.grid(False, which='minor')
    if title:
        ax.set_title(title, color=INK, loc='left', fontsize=8)


def handles(items):
    """legend keys: (label, colour, 'line' | 'dot' | 'bar')"""
    h = []
    for lab, c, kind in items:
        if kind == 'line':
            h.append(Line2D([], [], color=c, lw=2 * LW, label=lab))
        elif kind == 'dot':
            h.append(Line2D([], [], color=c, marker='o', ls='none', ms=MS, mec=SURF, mew=RING, label=lab))
        else:
            h.append(Patch(facecolor=c, edgecolor='none', label=lab))
    return h


HEAD = 0.62                                           # inches at the top of a figure for its title and legend


def head(fig, title, items):
    h = fig.get_figheight()
    fig.text(0.008, 1 - 0.1 / h, title, color=INK, fontsize=10, fontweight='bold', va='top')
    if items:
        leg = fig.legend(handles=handles(items), loc='upper left', bbox_to_anchor=(0.004, 1 - 0.3 / h), ncol=len(items), fontsize=8,
                         handlelength=1.6, columnspacing=1.8, borderaxespad=0)
        for t in leg.get_texts():
            t.set_color(INK2)


def layout(fig):
    fig.tight_layout(rect=(0, 0, 1, 1 - HEAD / fig.get_figheight()), w_pad=1.6, h_pad=1.4)


def end_labels(ax, xs, ys, labels):
    """names at the line ends, only if they do not collide (else the legend carries them)"""
    y0, y1 = ax.get_ylim(); ys = np.asarray(ys, float)
    if len(ys) > 1 and np.min(np.diff(np.sort(ys))) < 0.08 * (y1 - y0):
        return
    for x, y, l in zip(xs, ys, labels):
        ax.annotate(l, (x, y), xytext=(4, 0), textcoords='offset points', va='center', fontsize=7, color=INK2)


def _num(v, _=None):
    if v <= 0:
        return ''
    if 1e-2 <= v < 1e4:
        return f'{v:g}'
    e = int(np.floor(np.log10(v) + 1e-9)); m = v / 10 ** e
    return f'$10^{{{e}}}$' if abs(m - 1) < 1e-6 else f'${m:g}\\times10^{{{e}}}$'


def logaxes(ax, x, y):
    """log-log axes whose ticks are clean numbers, two or three per decade when the data span less than about a decade"""
    ax.set_xscale('log'); ax.set_yscale('log')
    for axis, d in ((ax.xaxis, x), (ax.yaxis, y)):
        d = np.asarray([v for v in np.ravel(d) if np.isfinite(v) and v > 0])
        dec = np.log10(d.max() / d.min()) if len(d) > 1 else 0
        axis.set_major_locator(LogLocator(base=10, subs=(1, 1.5, 2, 3, 4, 5, 6, 8) if dec < 0.4 else (1, 2, 3, 5) if dec < 0.8 else (1, 2, 5) if dec < 2 else (1.0,)))
        axis.set_major_formatter(FuncFormatter(_num)); axis.set_minor_formatter(NullFormatter())


def px_data(ax):
    bb = ax.get_window_extent(); (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    return (x1 - x0) / bb.width, (y1 - y0) / bb.height


def bar(ax, x, w, h, color):
    """a bar from the baseline 0 to h: 4 px rounded data end, square at the baseline (draw after the limits and the layout are fixed)"""
    if not np.isfinite(h) or h == 0:
        return
    dx, dy = px_data(ax)
    rx, ry = min(4 * dx, w / 2), min(4 * dy, abs(h)); s = 1 if h > 0 else -1; a, b = x - w / 2, x + w / 2
    v = [(a, 0), (a, h - s * ry), (a, h), (a + rx, h), (b - rx, h), (b, h), (b, h - s * ry), (b, 0), (a, 0)]
    c = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(v, c), facecolor=color, edgecolor='none', lw=0))


def dots(ax, x, y, c, **kw):
    ax.plot(x, y, 'o', ls='none', color=c, ms=MS, mec=SURF, mew=RING, zorder=3, **kw)


def csv(fn, header, rows):
    with open(fn, 'w') as f:
        f.write(','.join(header) + '\n')
        for r in rows:
            f.write(','.join(x if isinstance(x, str) else f'{x:.6g}' for x in r) + '\n')


def run_params(run):
    """dt (fs) and the -var options of the folder (config.env lines written by go.sh)."""
    opts = {}
    f = os.path.join(run, 'config.env')
    if os.path.exists(f):
        for m in re.finditer(r'-var\s+(\S+)\s+(\S+?)["\s]', open(f).read()):
            opts[m.group(1)] = m.group(2)
    return float(opts.get('dt', 1.0)), opts


def read_shear(fn):
    """shear_<s>.dat: step gamma tau sigma_n mu Lz rho T Lx Ly rate(1/s).  Lines of older inputs have fewer columns (the missing ones are nan)."""
    rows = []
    with open(fn) as f:
        for l in f:
            if l.startswith('#') or not l.strip():
                continue
            try:
                v = [float(x) for x in l.split()]
            except ValueError:
                continue
            rows.append((v + [np.nan] * 11)[:11])
    if len(rows) < 2:
        return None
    d = np.array(rows)
    sh = {c: d[:, i] for i, c in enumerate(['step', 'gamma', 'tau', 'sn', 'mu', 'lz', 'rho', 'T', 'lx', 'ly', 'rate'])}
    if np.nanmedian(sh['sn']) < 0:                                # an older input wrote the normal stress with the opposite sign
        sh['sn'] = -sh['sn']; sh['mu'] = sh['tau'] / sh['sn']
    return sh


def rl(rate):
    """1e8 -> '1e8 /s', 3.16e8 -> '3.2e8 /s'"""
    m, e = f'{rate:.1e}'.split('e')
    return f'{m.rstrip("0").rstrip(".")}e{int(e)} /s'


def segments(sh, dt, gmin=None):
    """the phases of the shear, one rate each (phase 0, then the velocity steps): rate (1/s), its rows, the rows of its steady state
    (the second half of its strain; phase 0 from gmin if given)"""
    st, g, rate = sh['step'], sh['gamma'], sh['rate'].copy()
    if not np.isfinite(rate).all():                              # an older input: one rate, from gamma = rate * dt * step
        rate[:] = np.polyfit(st, g, 1)[0] / dt * 1e15
    out, a = [], 0
    for b in range(1, len(st) + 1):
        if b == len(st) or abs(rate[b] / rate[a] - 1) > 1e-3:
            idx = np.arange(a, b)
            gm = gmin if (gmin is not None and not out) else g[a] + 0.5 * (g[b - 1] - g[a])
            out.append(dict(rate=rate[a], rows=idx, win=idx[g[idx] >= gm], s0=st[a]))
            a = b
    return out


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


def tilt_of(sh):
    """cumulative tilt X(step) = integral of Lz d(gamma) (the box flips do not matter for unwrapped positions; any sequence of rates)"""
    st0 = np.r_[0.0, sh['step']]; lz0 = np.r_[sh['lz'][0], sh['lz']]; g0 = np.r_[0.0, sh['gamma']]
    X = np.r_[0.0, np.cumsum(0.5 * (lz0[1:] + lz0[:-1]) * np.diff(g0))]
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


def interfaces(run, sig, sh, dt, args):
    """stacked particle pairs at the first frame and their relative in-plane displacement s(t) over the whole shear (all phases)"""
    steps, P = read_avetime(job_files(run, 'particles', sig))
    if steps is None or len(steps) < 3:
        return None
    if not np.isfinite(sh['lx']).any():
        print(f'  {run}: shear_{sig}.dat has no Lx/Ly columns (older in.04_shear): no interface analysis'); return None
    X = tilt_of(sh)
    ok = np.isfinite(sh['lx'])
    lx = np.interp(steps, sh['step'][ok], sh['lx'][ok]); ly = np.interp(steps, sh['step'][ok], sh['ly'][ok]); lz = np.interp(steps, sh['step'], sh['lz'])
    com, gyr = P[:, :, 0:3], P[:, :, 6:12]
    n0, t0, R0 = particle_frame(gyr[0])
    H = lambda k: np.array([[lx[k], 0, X(steps[k])], [0, ly[k], 0], [0, 0, lz[k]]])
    H0inv = np.linalg.inv(H(0))
    N = com.shape[1]
    pairs = []
    for i in range(N):
        for j in range(i + 1, N):
            if abs(n0[i] @ n0[j]) < np.cos(np.radians(args.angle)):
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
    out = []
    for (i, j, img, nb, d0, h) in pairs:
        d = np.array([com[k, j] - com[k, i] - H(k) @ img for k in range(len(steps))])
        dd = d - d0
        s = dd - np.outer(dd @ nb, nb)                           # in-plane relative displacement (A)
        # a uniform shear moves the two centres apart by gamma * dz in x: its in-plane part per unit strain (A)
        a1 = np.array([d0[2], 0.0, 0.0]); a1 = a1 - (a1 @ nb) * nb
        out.append(dict(i=i, j=j, h=h, s=s, smag=np.linalg.norm(s, axis=1), a1n=np.linalg.norm(a1)))
    return dict(t=steps * dt * 1e-6, steps=steps, pairs=out)


def slip(I, s0, s1, gdot):
    """slip velocity (m/s; slope of |s|) of every interface between the steps s0 and s1 (one phase), and what a uniform shear would give"""
    m = (I['steps'] >= s0) & (I['steps'] <= s1)
    if m.sum() < 3:
        return None
    return [dict(n=n, i=p['i'], j=p['j'], h=p['h'], v=np.polyfit(I['t'][m], p['smag'][m], 1)[0] * 1e-6 * V_MS, vaff=gdot * p['a1n'] * 1e-10)
            for n, p in enumerate(I['pairs'])]


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


def water_frames(run, sig, args):
    """(file, offset) of the frames with water: the all-atom snapshots shear_<s>_all_<step>.dump (start and end of each phase),
    else the first and last frame of shear_<s>.dump (env DUMP=all ./go.sh ...)"""
    snaps = glob.glob(os.path.join(run, f'shear_{sig}_all_*.dump'))
    if snaps:
        return [(f, 0) for f in sorted(snaps, key=lambda f: int(f.rsplit('_', 1)[1].split('.')[0]))]
    fn = os.path.join(run, f'shear_{sig}.dump')
    if not os.path.exists(fn):
        return []
    offs = dump_frames(fn)
    sel = sorted(set([0, len(offs) - 1] + [k for k in args.frames if 0 <= k < len(offs)])) if offs else []
    return [(fn, offs[k]) for k in sel]


def water_frame(fn, off, names, args):
    """one frame: class of each water molecule (between two particles / at one surface / pore), H-bonds per molecule, g_OO(r) of each class"""
    from scipy.spatial import cKDTree
    tow = [t for t, n in names.items() if n.lower() == 'ow'] or [20]
    thw = [t for t, n in names.items() if n.lower() == 'hw'] or [19]
    tna = [t for t, n in names.items() if n.lower() == 'na']
    step, lo, H, c = read_frame(fn, off)
    typ = c['type'].astype(int); mol = c['mol'].astype(int)
    r = np.c_[c['xu'], c['yu'], c['zu']]
    isO = np.isin(typ, tow); isH = np.isin(typ, thw)
    if isO.sum() == 0:
        return dict(step=step, nowater=True)
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
    nb = cKDTree(Cimg).query_ball_point(O, args.rc)
    npart = np.array([len(set(cmol[l])) for l in nb])
    cls = np.where(npart >= 2, 0, np.where(npart == 1, 1, 2))          # 0 between particles, 1 at one surface, 2 pore
    # hydrogen bonds (Zhou et al.: O-O < 3.5 A, H...O < 2.45 A, angle H-Od...Oa < 30 deg)
    nhb = np.zeros(len(O))
    cos30 = np.cos(np.radians(30))
    for i, l in enumerate(tO.query_ball_point(O, 3.5)):
        l = np.array(l)
        if len(l) == 0:
            continue
        d = Oimg[l] - O[i]; dn = np.linalg.norm(d, axis=1)
        keep = dn > 0.5; l, d, dn = l[keep], d[keep], dn[keep]
        for h in vOH[i]:
            hb = (np.linalg.norm(d - h, axis=1) < 2.45) & ((d @ h) / (dn * np.linalg.norm(h)) > cos30)
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
        hh = np.zeros(len(rr))
        for i in idx:
            dist = np.linalg.norm(Oimg[far[i]] - O[i], axis=1)
            hh += np.histogram(dist[dist > 0.5], edges)[0]
        g[ccl] = (hh / (len(idx) * 4 * np.pi * rr ** 2 * dr * rho)).tolist()
    return dict(step=int(step), nhb={ccl: float(nhb[cls == ccl].mean()) for ccl in (0, 1, 2) if (cls == ccl).any()},
                frac={ccl: float((cls == ccl).mean()) for ccl in (0, 1, 2)}, r=rr.tolist(), g=g)


def water_cached(fn, off, names, args, cache_dir, allowed):
    """water_frame, kept in cache_dir (a frame is analysed once: the figures after every job stay quick); allowed = False: only from the cache"""
    st = os.stat(fn)
    cf = os.path.join(cache_dir, re.sub(r'[^A-Za-z0-9._-]', '_', f'{os.path.abspath(fn)}_{off}')[-150:] + '.json')
    key = [st.st_size, int(st.st_mtime), args.rc, args.rmax]
    if os.path.exists(cf):
        try:
            d = json.load(open(cf))
            if d.get('key') == key:
                w = d['w']
                for k in ('nhb', 'frac', 'g'):
                    if k in w:
                        w[k] = {int(a): b for a, b in w[k].items()}
                return w
        except Exception:
            pass
    if not allowed:
        return 'later'
    try:
        w = water_frame(fn, off, names, args)
    except (ValueError, IndexError) as e:                          # a snapshot still being written
        print(f'  {fn}: not readable ({e}): skipped'); return None
    os.makedirs(cache_dir, exist_ok=True)
    json.dump(dict(key=key, w=w), open(cf, 'w'))
    return w


def profiles(run, sig, s0, s1, V):
    """profiles across the gouge averaged between the steps s0 and s1 (a steady state); the velocity as a fraction of the uniform shear:
    u = vx / (gamma_dot Lz), V = gamma_dot Lz in A/fs"""
    steps, P = read_avetime(job_files(run, 'prof', sig))
    if steps is None:
        return None
    m = (steps >= s0) & (steps <= s1)
    if m.sum() == 0:
        return None
    D = P[m]                                                     # coord (reduced), n_water, vx_water, n_clay, vx_clay, T
    nw, nc = D[:, :, 1], D[:, :, 3]
    vw = (D[:, :, 2] * nw).sum(0) / np.maximum(nw.sum(0), 1e-9)
    vc = (D[:, :, 4] * nc).sum(0) / np.maximum(nc.sum(0), 1e-9)
    return dict(zr=D[0, :, 0], uw=vw / V, uc=vc / V, share=nw.sum(0) / np.maximum((nw + nc).sum(0), 1e-9), T=D[:, :, 5].mean(0))


def fig6(combos, fn, out_csv):
    rows_ = []
    for cb in combos:
        seg = next((r for r in cb['segs'] if r.get('ifc')), None)      # phase 0 (the main rate) when it has enough frames
        if cb['I'] and seg:
            rows_.append((cb, seg))
    if not rows_:
        return
    fig, ax = plt.subplots(len(rows_), 4, figsize=(15, 2.6 * len(rows_) + HEAD), squeeze=False)
    comp = [('along the flow (x)', SLOT[0]), ('across the flow (y)', SLOT[1])]
    hist, rows = [], []
    for row, (cb, seg) in enumerate(rows_):
        I = cb['I']; key = f'{cb["run"]} {cb["sig"]} MPa'
        ps = sorted(seg['ifc'], key=lambda p: -p['v'])
        idx = sorted({0, len(ps) // 2, len(ps) - 1}); names = {0: 'fastest', len(ps) // 2: 'median', len(ps) - 1: 'slowest'}
        starts = [r['seg_start'] * cb['dt'] * 1e-6 for r in cb['segs'][1:]]   # the velocity steps begin (ns)
        for c in range(3):
            a = ax[row, c]
            if c >= len(idx):
                a.axis('off'); continue
            p = ps[idx[c]]; s = I['pairs'][p['n']]['s']
            ys = [s[:, 0] * A2NM, s[:, 1] * A2NM]
            for t0 in starts:
                a.axvline(t0, color=AXIS, lw=PX, zorder=1)
            for y, (lab, col) in zip(ys, comp):
                a.plot(I['t'], y, color=col)
            a.set_xlim(I['t'][0], I['t'][-1] + 0.12 * (I['t'][-1] - I['t'][0]))
            style(a, f'{key} · {names[idx[c]]} of {len(ps)} interfaces (particles {p["i"]+1}-{p["j"]+1})\n' + (
                     f'slip {p["v"]:.3g} m/s at {seg["rl"]}' if p['v'] > 0 else f'no steady slip at {seg["rl"]} (stuck)') + f', film {p["h"]:.1f} Å')
            a.set_xlabel('time (ns)' + ('; gray lines: velocity steps' if starts else '')); a.set_ylabel('in-plane slip (nm)')
            end_labels(a, [I['t'][-1]] * 2, [y[-1] for y in ys], ['x', 'y'])
        # how much faster than the gouge the film between the two particles is sheared
        f = np.array([abs(p['v']) / (p['h'] * 1e-10) / seg['gdot'] for p in ps])
        hist.append((ax[row, 3], np.log10(np.clip(f, 1e-3, None)), f'{key}, {seg["rl"]}', f))
        for r in cb['segs']:
            for p in (r['ifc'] or []):
                fi = abs(p['v']) / (p['h'] * 1e-10) / r['gdot']
                rows.append([cb['run'], cb['sig'], r['rate'], p['i'] + 1, p['j'] + 1, p['h'], p['v'], p['vaff'], fi, abs(p['v']) / p['vaff'] if p['vaff'] > 0 else np.nan])
    head(fig, 'Slip between stacked clay particles (cf. Zhou et al. 2023, Fig. 6)', [(l, c, 'line') for l, c in comp])
    for a, lf, key, f in hist:                                   # limits first, then the layout, then the bars (their rounding is in pixels)
        nb = int(min(24, max(6, len(lf) // 3)))
        cnt, edges = np.histogram(lf, bins=nb, range=(min(lf.min(), -0.5), max(lf.max(), 0.5)))
        a.set_xlim(edges[0] - 0.05 * (edges[-1] - edges[0]), edges[-1] + 0.05 * (edges[-1] - edges[0])); a.set_ylim(0, cnt.max() * 1.18)
        a._hs = (cnt, edges)
        style(a, f'{key} · film rate / gouge rate\n{np.mean(f > 1) * 100:.0f} % of the films shear faster than the gouge')
        a.set_xlabel('log10 [ (v / h) / γ̇ ]'); a.set_ylabel('interfaces')
    layout(fig)
    for a, lf, key, f in hist:
        cnt, edges = a._hs; dx, dy = px_data(a); bw = edges[1] - edges[0]; w = min(bw - 2 * dx, 24 * dx)
        for k, n in enumerate(cnt):
            bar(a, 0.5 * (edges[k] + edges[k + 1]), w, n, SLOT[0])
        a.axvline(0, color=INK2, lw=PX, zorder=4)
        a.annotate('gouge rate', (0, a.get_ylim()[1]), xytext=(3, -2), textcoords='offset points', va='top', fontsize=7, color=INK2)
    fig.savefig(fn); plt.close(fig)
    csv(out_csv, ['run', 'sigma_MPa', 'gouge_rate_1_s', 'particle_i', 'particle_j', 'film_h_A', 'slip_v_m_s', 'uniform_shear_v_m_s',
                  'film_rate_over_gouge_rate', 'v_over_uniform'], rows)


def fig7(res, entcol, fn, out_csv):
    keys = [k for k, r in res.items() if r.get('prof')]
    if not keys:
        return
    one_run = len({res[k]['run'] for k in keys}) == 1
    sigs = sorted({res[k]['sig'] for k in keys}, key=float)
    fig, ax = plt.subplots(len(sigs), 4, figsize=(14, 3.3 * len(sigs) + HEAD), sharey=True, squeeze=False)
    rows, seen = [], {}
    for row, sig in enumerate(sigs):
        for k in [k for k in keys if res[k]['sig'] == sig]:
            r = res[k]; p = r['prof']; c = entcol[r['ent']]
            seen[r['rl'] if one_run else f'{r["run"]} {r["rl"]}'] = c
            ax[row, 0].plot(p['uw'], p['zr'], color=c)
            ax[row, 1].plot(p['uc'], p['zr'], color=c)
            ax[row, 2].plot(p['share'], p['zr'], color=c)
            ax[row, 3].plot(p['T'], p['zr'], color=c)
            for z, uw, uc, s, T in zip(p['zr'], p['uw'], p['uc'], p['share'], p['T']):
                rows.append([r['run'], sig, r['rate'], z, uw, uc, s, T])
        for a in ax[row, :2]:
            a.plot([0, 1], [0, 1], color=AXIS, lw=1.5 * LW, zorder=1)
        style(ax[row, 0], f'σn = {sig} MPa · water velocity\ngray diagonal: uniform shear; a step: a shear band')
        style(ax[row, 1], f'σn = {sig} MPa · clay velocity')
        style(ax[row, 2], f'σn = {sig} MPa · where the water is')
        style(ax[row, 3], f'σn = {sig} MPa · temperature (flow removed)')
        ax[row, 0].set_xlabel('vx / (γ̇ Lz)'); ax[row, 1].set_xlabel('vx / (γ̇ Lz)')
        ax[row, 2].set_xlabel('water atoms / all atoms of the layer'); ax[row, 3].set_xlabel('T (K)')
        ax[row, 0].set_ylabel('z / Lz'); ax[row, 0].set_ylim(0, 1)
        for a in ax[row]:
            a.xaxis.set_major_locator(MaxNLocator(5))
    items = [(lab, c, 'line') for lab, c in seen.items()]
    title = 'Velocity profile across the gouge, steady state of each shear rate (cf. Zhou et al. 2023, Fig. 7)'
    head(fig, title + (f' · {res[keys[0]]["run"]}' if one_run else ''), items)
    layout(fig); fig.savefig(fn); plt.close(fig)
    csv(out_csv, ['run', 'sigma_MPa', 'gouge_rate_1_s', 'z_over_Lz', 'vx_water_over_gdotLz', 'vx_clay_over_gdotLz', 'water_share', 'T_K'], rows)


def fig8_10(res, fn8, fn10, out_csv, summary):
    by_sig = {}
    for key, r in res.items():
        by_sig.setdefault(r['sig'], []).append((key, r))
    sigs = sorted(by_sig, key=float)
    wid = max(6.4, 4.6 * len(sigs) + 0.4)
    fig8, ax8 = plt.subplots(1, len(sigs), figsize=(wid, 3.6 + HEAD), squeeze=False)
    fig10, ax10 = plt.subplots(1, len(sigs), figsize=(wid, 3.6 + HEAD), squeeze=False)
    rows = []
    for c, sig in enumerate(sigs):
        lst = by_sig[sig]
        V, tau, vh, gd, gdf = [], [], [], [], []
        for key, r in lst:
            V.append(r['V']); tau.append(r['tau']); gd.append(r['gdot'])
            if r.get('ifc'):
                ps = sorted(r['ifc'], key=lambda p: -abs(p['v']))
                k = max(1, len(ps) // 10)                           # the fastest 10 % (at least one): the hot spots
                vh.append(np.mean([abs(p['v']) for p in ps[:k]])); gdf.append(np.mean([abs(p['v']) / (p['h'] * 1e-10) for p in ps[:k]]))
            else:
                vh.append(np.nan); gdf.append(np.nan)
        V, tau, vh, gd, gdf = map(np.array, (V, tau, vh, gd, gdf))
        for (key, r), a, b, e, f, g in zip(lst, V, tau, vh, gd, gdf):
            rows.append([r['run'], sig, f, a, b, r['sn'], e, g, b * 1e6 / f, b * 1e6 / g if np.isfinite(g) else np.nan])
        a = ax8[0, c]; ks = []
        logaxes(a, np.r_[V, vh], tau)
        for x, col, lab in ((V, C_GOUGE, 'gouge'), (vh, C_IFC, 'hot interfaces')):
            m = np.isfinite(x) & (x > 0) & (tau > 0)
            if not m.any():
                continue
            dots(a, x[m], tau[m], col)
            if m.sum() >= 2:
                k, lti = np.polyfit(np.log(x[m]), np.log(tau[m]), 1)
                xx = np.logspace(np.log10(x[m].min()), np.log10(x[m].max()), 20)
                a.plot(xx, np.exp(lti) * xx ** k, color=col, lw=LW, alpha=0.6, zorder=2)
                ks.append(f'{lab} {k:.3f}')
                summary.append(f'sigma {sig} {lab}: ln tau = {k:.3f} ln v + {lti:.3f}  (k = {k:.3f}, tau_i = {np.exp(lti):.3g} MPa at 1 m/s)')
        a.margins(x=0.25, y=0.2)
        style(a, f'σn = {sig} MPa · lines: ln τ = k ln v + ln τi' + ('\nk: ' + ', '.join(ks) if ks else '')); a.set_xlabel('slip velocity v (m/s)'); a.set_ylabel('shear stress τ (MPa)')
        a = ax10[0, c]
        logaxes(a, np.r_[gd, gdf], np.r_[tau * 1e6 / gd * 1e3, tau * 1e6 / gdf * 1e3])
        dots(a, gd, tau * 1e6 / gd * 1e3, C_GOUGE)
        m = np.isfinite(gdf) & (gdf > 0)
        if m.any():
            dots(a, gdf[m], tau[m] * 1e6 / gdf[m] * 1e3, C_IFC)
        a.margins(x=0.25, y=0.2)
        style(a, f'σn = {sig} MPa'); a.set_xlabel('shear rate (1/s)'); a.set_ylabel('effective viscosity τ / rate (mPa s)')
    head(fig8, 'Shear stress against slip velocity (cf. Zhou et al. 2023, Fig. 8)',
         [('gouge: V = γ̇ Lz', C_GOUGE, 'dot'), ('hot interfaces: fastest 10 %', C_IFC, 'dot')])
    head(fig10, 'Effective viscosity against shear rate (cf. Zhou et al. 2023, Fig. 10)',
         [('gouge: τ / γ̇', C_GOUGE, 'dot'), ('sliding films: τ / (v / h)', C_IFC, 'dot')])
    layout(fig8); fig8.savefig(fn8); layout(fig10); fig10.savefig(fn10)
    plt.close(fig8); plt.close(fig10)
    csv(out_csv, ['run', 'sigma_MPa', 'gouge_rate_1_s', 'V_m_s', 'tau_MPa', 'sigma_n_MPa', 'hot_interface_v_m_s', 'film_rate_1_s',
                  'viscosity_gouge_Pa_s', 'viscosity_film_Pa_s'], rows)


CLS = {0: 'between two particles', 1: 'at one surface', 2: 'pore'}


def fig9(combos, fn, out_csv, summary):
    cbs = [cb for cb in combos if cb.get('water')]
    if not cbs:
        return
    ng = max(min(2, len(cb['water'])) for cb in cbs)
    fig, ax = plt.subplots(len(cbs), 1 + ng, figsize=(max(4.4 * (1 + ng), 4.4 + 0.9 * max(len(cb['water']) for cb in cbs)), 3.2 * len(cbs) + HEAD),
                           squeeze=False)
    rows, barsets = [], []
    for row, cb in enumerate(cbs):
        W = cb['water']; gof = cb['gamma_of']; key = f'{cb["run"]} {cb["sig"]} MPa'
        lab = lambda fr: 'start' if fr['step'] == 0 else f'γ {gof(fr["step"]):.2f}' + (f'\n{rl(cb["rate_of"](fr["step"]))}' if cb['rate_of'](fr['step']) > 0 else '')
        for fr in W:
            gam = gof(fr['step'])
            summary.append(f'{key} gamma {gam:.2f}: H-bonds per water ' + ', '.join(f'{CLS[c]} {fr["nhb"][c]:.2f} ({fr["frac"][c]*100:.0f} %)' for c in fr['nhb']))
            for c in (0, 1, 2):
                rows.append([cb['run'], cb['sig'], fr['step'], gam, cb['rate_of'](fr['step']), CLS[c], fr['frac'][c], fr['nhb'].get(c, np.nan)])
        # bars: one group per snapshot (strain; the rate of the phase it ends), one bar per class of water
        a = ax[row, 0]
        top = max(max(fr['nhb'].values()) for fr in W)
        a.set_xlim(-0.6, len(W) - 0.4); a.set_ylim(0, top * 1.18)
        a.set_xticks(range(len(W))); a.set_xticklabels([lab(fr) for fr in W], fontsize=6.5)
        style(a, f'{key} · hydrogen bonds per water molecule\nO–O < 3.5 Å, H···O < 2.45 Å, angle < 30°'); a.set_ylabel('H-bonds per molecule')
        barsets.append((a, W))
        show = [W[0]] + ([W[-1]] if len(W) > 1 else [])
        for n, fr in enumerate(show):
            a = ax[row, 1 + n]
            for c, g in fr['g'].items():
                a.plot(fr['r'], g, color=SLOT[c])
            a.axhline(1, color=AXIS, lw=PX, zorder=1)
            style(a, f'{key} · g_OO(r), {lab(fr)}'.replace('\n', ', ')); a.set_xlabel('r (Å)'); a.set_ylabel('g_OO(r)')
        for n in range(len(show), ng):
            ax[row, 1 + n].axis('off')
        if ng == 2 and len(show) == 2:                          # the same scale for the two frames
            y1 = max(ax[row, 1].get_ylim()[1], ax[row, 2].get_ylim()[1])
            ax[row, 1].set_ylim(0, y1); ax[row, 2].set_ylim(0, y1)
    head(fig, 'Water between and around the clay particles (cf. Zhou et al. 2023, Fig. 9)', [(CLS[c], SLOT[c], 'bar') for c in (0, 1, 2)])
    layout(fig)
    for a, W in barsets:
        dx, dy = px_data(a); w = min(0.26, 24 * dx); gap = 2 * dx
        for n, fr in enumerate(W):
            for c in (0, 1, 2):
                bar(a, n + (c - 1) * (w + gap), w, fr['nhb'].get(c, np.nan), SLOT[c])
    fig.savefig(fn); plt.close(fig)
    csv(out_csv, ['run', 'sigma_MPa', 'step', 'gamma', 'rate_1_s', 'water_class', 'fraction_of_water', 'hbonds_per_molecule'], rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--sigma', nargs='*', default=None, help='normal stresses to look at (default: every shear_<s>.dat)')
    ap.add_argument('--gmin', type=float, default=None, help='phase 0: steady state from this strain (default: the second half of each phase)')
    ap.add_argument('--hmax', type=float, default=12.0, help='A: largest water film between two stacked particles')
    ap.add_argument('--angle', type=float, default=25.0, help='deg: largest angle between the normals of stacked particles')
    ap.add_argument('--overlap', type=float, default=0.75, help='lateral offset of stacked particles < overlap * (R1 + R2)')
    ap.add_argument('--rc', type=float, default=6.0, help='A: a water O within rc of clay atoms of 2 particles is "between two particles"')
    ap.add_argument('--rmax', type=float, default=8.0, help='A: range of g_OO(r)')
    ap.add_argument('--frames', type=int, nargs='*', default=[], help='more frames of shear_<s>.dump for fig9 when there are no snapshots (first and last always)')
    ap.add_argument('--budget', type=float, default=0, help='s: no new trajectory frame for fig9 is started after half of this time (0: no limit; job.sh uses it)')
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()
    t_start = time.time()
    os.makedirs(args.outdir, exist_ok=True)
    o = lambda n: os.path.join(args.outdir, n)
    res, combos, summary, ents = {}, [], [], []
    for run in args.runs:
        run = run.rstrip('/') or '/'
        rn = os.path.basename(os.path.abspath(run))
        dt, opts = run_params(run)
        sigs = args.sigma or sorted({re.match(r'shear_(.+)\.dat$', os.path.basename(f)).group(1) for f in glob.glob(os.path.join(run, 'shear_*.dat'))}, key=float)
        for sig in sigs:
            fn = os.path.join(run, f'shear_{sig}.dat')
            if not os.path.exists(fn):
                print(f'{fn}: not there'); continue
            sh = read_shear(fn)
            if sh is None:
                print(f'{fn}: too few lines yet'); continue
            segs = segments(sh, dt, args.gmin)
            I = interfaces(run, sig, sh, dt, args)
            st0, g0, r0 = np.r_[0.0, sh['step']], np.r_[0.0, sh['gamma']], np.r_[np.nan, [s['rate'] for s in segs for _ in s['rows']]]
            cb = dict(run=rn, path=run, sig=sig, dt=dt, I=I, segs=[],
                      gamma_of=lambda st, a=st0, b=g0: float(np.interp(st, a, b)),
                      rate_of=lambda st, a=st0, b=r0: float(np.nan_to_num(b[max(0, np.searchsorted(a, st, 'right') - 1)])))
            print(f'{rn} {sig} MPa: strain up to {sh["gamma"].max():.3f}, {len(segs)} rate(s): ' + ', '.join(rl(s['rate']) for s in segs))
            for k, sg in enumerate(segs):
                w = sg['win']
                if len(w) < 2:
                    print(f'  {rl(sg["rate"])}: too short yet for a steady state'); continue
                gdot = sg['rate']; lz = sh['lz'][w].mean()
                r = dict(run=rn, sig=sig, k=k, rate=gdot, rl=rl(gdot), tau=sh['tau'][w].mean(), sn=sh['sn'][w].mean(), lz=lz, gdot=gdot,
                         V=gdot * lz * 1e-10, s0=sh['step'][w[0]], s1=sh['step'][w[-1]], seg_start=sg['s0'])
                r['ent'] = (rn, r['rl'])
                if r['ent'] not in ents:
                    ents.append(r['ent'])
                r['ifc'] = slip(I, r['s0'], r['s1'], gdot) if I else None
                r['prof'] = profiles(run, sig, r['s0'], r['s1'], gdot * 1e-15 * lz)
                key = f'{rn} {sig} MPa {r["rl"]}'
                summary.append(f'{key}: tau {r["tau"]:.4g} MPa  sigma_n {r["sn"]:.4g} MPa  mu {r["tau"]/r["sn"]:.4f}  V {r["V"]:.4g} m/s  Lz {lz:.1f} A'
                               f'  (strain {sh["gamma"][w[0]]:.3f}-{sh["gamma"][w[-1]]:.3f})')
                if r['ifc']:
                    v = np.array([abs(p['v']) for p in r['ifc']])
                    summary.append(f'  {len(v)} stacked interfaces: slip velocity median {np.median(v):.3g} m/s, fastest {v.max():.3g} m/s (gouge V {r["V"]:.3g} m/s)')
                res[key] = r; cb['segs'].append(r)
            combos.append(cb)
    if not res:
        sys.exit('nothing to analyse yet')
    # one colour per run and rate, in the order met, the same in every figure; past eight: gray
    entcol = {e: (SLOT[k] if k < len(SLOT) else MUTED) for k, e in enumerate(ents)}
    fig6(combos, o('fig6_interface_slip.png'), o('fig6_interfaces.csv'))
    fig7(res, entcol, o('fig7_velocity_profile.png'), o('fig7_profiles.csv'))
    fig8_10(res, o('fig8_tau_v.png'), o('fig10_viscosity.png'), o('fig8_fig10_points.csv'), summary)
    open(o('hot_spots_summary.txt'), 'w').write('\n'.join(summary) + '\n')
    # fig9 last: the water of the snapshots takes the longest (a frame is analysed once and kept in <outdir>/.cache)
    try:
        import scipy  # noqa: F401
        later = 0
        for cb in combos:
            names = type_names(cb['path']); W = []
            for f, off in water_frames(cb['path'], cb['sig'], args):
                w = water_cached(f, off, names, args, o('.cache'), not args.budget or time.time() - t_start < 0.5 * args.budget)
                if w == 'later':
                    later += 1
                elif w and w.get('nowater'):
                    print(f'  {f}: no water in this trajectory (the snapshots shear_<s>_all_<step>.dump have it)'); break
                elif w:
                    W.append(w)
            cb['water'] = W
        if later:
            print(f'  {later} trajectory frame(s) for fig9 left for the next time (time budget)')
        fig9(combos, o('fig9_water_structure.png'), o('fig9_water.csv'), summary)
    except ImportError:
        print('  scipy missing: no fig9')
    open(o('hot_spots_summary.txt'), 'w').write('\n'.join(summary) + '\n')
    print('\n'.join(summary)); print(f'figures and their numbers (.csv) in {os.path.abspath(args.outdir)}  ({time.time() - t_start:.0f} s)')


if __name__ == '__main__':
    main()
