"""Analyses of a soapff-cui MD run, computed by established tools (no analysis algorithm of soapff-cui itself):

  msd       LAMMPS `compute msd` (in in.lammps, written every 0.1 ps to msd_<group>.dat) of water O and each ion type;
            self-diffusion coefficient from the Einstein relation, D = slope / (2 d), slope of a least-squares line
            over the production part: total (d = 3), in plane (xy, d = 2), normal (z, d = 1).  No finite-size correction.
            Ionic conductivity (Nernst-Einstein, sigma = e^2 / (V kB T) sum_i N_i z_i^2 D_i, ion-ion correlations neglected).
  PDF, RDF  SOVA (sovapy, soapff_cui/sova_structure.py): pair distribution functions g(r) (atom_pair_hist,
            pair_dist_func) and radial distribution functions N(r) = 4 pi r^2 rho g(r) (radial_dist_fun) of water O -
            water O, ion - water O, ion - mineral O and every element pair of the whole system; coordination numbers
            (N(r) integrated to the first minimum); X-ray structure factor S(Q) of the whole system
            (partial_structure_factor, xcoeff, structure_factor_xray)
  compute   the output files of the averaging fixes in in.lammps (fix ave/time, ave/chunk, ave/correlate, ave/histo with
            `file`, i.e. the computes added from the catalog): scalar columns against time and their mean over the
            production part, vector / chunk / histogram blocks averaged over the production part
  swelling  the swelling curve written by the soapff swelling block (swelling.dat; Cygan et al. 2004): for every water
            content the mean and standard deviation of d001 = lz and of the cell over the averaging part of its NPT, against
            the water content in g H2O / g clay (clay = every atom that is not water); swelling and deswelling separately;
            an experimental curve can be drawn with it from swelling_exp.csv (water_g_per_g, d001_A[, label])
  XRD       Debyer (Debye scattering equation; https://github.com/wojdyr/debyer): X-ray powder pattern I(2 theta) averaged
            over the last frames of traj.lammpstrj.  cluster mode (default; crystals, minerals): the cell repeated to a
            crystallite with edges >= 50 A, without periodicity; periodic mode (liquids, glasses; orthorhombic box only,
            as Debyer's periodic boundaries are orthorhombic): cut-off half the shortest box length with sinc damping and
            the number density of the box

  python -m soapff_cui analyze OUTDIR [--skip-ps X] [--sova-python PY] [--sova-home DIR] [--what msd,sova,compute,xrd]
                                      [--debyer PATH] [--wavelength 1.5406] [--two-theta 3,80,0.02] [--xrd-frames 3]
                                      [--xrd-mode cluster|periodic] [--xrd-size 50]
      -> OUTDIR/analysis/{summary.txt, summary.json, analysis.png, msd_fit.csv, sova_pdf.csv, sova_rdf.csv, sova_sq.csv,
                          sova_summary.json, sova_sq.png, compute_<file>.png, compute_<file>_average.csv, xrd.png, xrd.csv}

SOVA: used directly when sovapy is importable; otherwise sova_structure.py is run with the Python of a SOVA installation
(--sova-python or $SOVA_PYTHON; --sova-home or $SOVA_HOME is put on its PYTHONPATH for a source checkout).
Debyer: --debyer, $DEBYER, debyer on PATH, or a source build ~/debyer/debyer/debyer or ~/*/debyer/debyer/debyer.
"""
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np


A2PS_TO_1E9M2S = 10.0          # 1 A^2/ps = 1e-8 m^2/s = 10 x 1e-9 m^2/s
WHAT = ("msd", "sova", "xrd")      # offered analyses; "compute" (fix ave/*) and "swelling" are kept in the code, not offered
AVE_FIXES = ("ave/time", "ave/chunk", "ave/correlate", "ave/histo")


def msd(out_dir, skip_ps, dt_fs):
    """D from the LAMMPS compute msd output (Einstein relation, least-squares slope after skip_ps)."""
    res, curves = {}, {}
    for path in sorted(glob.glob(os.path.join(out_dir, "msd_*.dat"))):
        g = os.path.basename(path)[4:-4]
        a = np.loadtxt(path, comments="#", ndmin=2)
        if len(a) < 3:
            continue
        t = a[:, 0] * dt_fs / 1000.0
        m = t >= skip_ps
        if m.sum() < 3:
            continue

        def D(y, dim):
            return float(np.polyfit(t[m], y[m], 1)[0] / (2 * dim) * A2PS_TO_1E9M2S)
        res[g] = dict(D_total=D(a[:, 4], 3), D_xy=D(a[:, 1] + a[:, 2], 2), D_z=D(a[:, 3], 1),
                      fit_from_ps=float(t[m][0]), fit_to_ps=float(t[m][-1]), msd_at_end_A2=float(a[-1, 4]))
        curves[g] = (t, a[:, 1:5], (float(t[m][0]), float(t[m][-1])))
    return res, curves


def sova(out_dir, dest, skip_ps, dt_fs, sova_python=None, sova_home=None, log=print, traj="traj.lammpstrj"):
    """run sova_structure.py (SOVA PDF, RDF, S(Q)) in this Python if sovapy is importable, else with the given SOVA Python."""
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sova_structure.py")
    args = ["--traj", os.path.join(out_dir, traj), "--charges", os.path.join(out_dir, "charges.csv"),
            "--out", dest, "--skip-ps", str(skip_ps), "--dt-fs", str(dt_fs)]
    sova_python = sova_python or os.environ.get("SOVA_PYTHON")
    sova_home = sova_home or os.environ.get("SOVA_HOME")
    if sova_python:
        py = sova_python
    else:
        try:
            import sovapy  # noqa: F401
            py = sys.executable
        except ImportError:
            raise RuntimeError("SOVA (sovapy) not found: pip install sovapy, or give the Python of a SOVA installation "
                               "(--sova-python / $SOVA_PYTHON, and --sova-home / $SOVA_HOME for a source checkout)")
    env = dict(os.environ)
    if sova_home:
        env["PYTHONPATH"] = sova_home + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("MPLCONFIGDIR", os.path.join(dest, ".mplconfig"))
    log(f"SOVA: {py} {os.path.basename(script)}")
    p = subprocess.run([py, script] + args, env=env, capture_output=True, text=True, cwd=dest)
    if p.returncode != 0:
        raise RuntimeError("SOVA failed: " + (p.stderr or p.stdout)[-1500:])
    return json.load(open(os.path.join(dest, "sova_summary.json")))



# --------------------------------------------------------------------------------------------- ionic conductivity
E, KB, A2M = 1.602176634e-19, 1.380649e-23, 1.0e-10          # C, J/K, m per A


def ionic_conductivity(out_dir, D, temp=None):
    """Nernst-Einstein conductivity (S/m) of the ions with MSD results D {group: dict(D_total, D_xy, D_z)} (1e-9 m2/s);
    charges and numbers from charges.csv, T the given temperature, V the box of system.data."""
    rows = list(csv.DictReader(open(os.path.join(out_dir, "charges.csv"))))
    ions = {}
    for r in rows:
        if r.get("role") == "ion":
            ions.setdefault(r["element"], []).append(float(r.get("charge", r.get("q", 0.0)) or 0.0))
    T, V = temp, _data_volume(out_dir)
    if not ions or not T or not V:
        return None
    pref = E ** 2 / (V * A2M ** 3 * KB * T)
    res = dict(T=T, V=V, ions={}, sigma_total=0.0, sigma_xy=0.0, sigma_z=0.0)
    for el, qs in ions.items():
        if el not in D:
            continue
        z2 = float(np.mean(np.square(qs)))                                    # e^2 (the SOAP charge of the ion)
        n = len(qs)
        part = {k: pref * n * z2 * D[el][f"D_{k}"] * 1e-9 for k in ("total", "xy", "z")}
        res["ions"][el] = dict(N=n, q=float(np.mean(qs)), **{f"sigma_{k}": v for k, v in part.items()})
        res["sigma_total"] += part["total"]
        res["sigma_xy"] += part["xy"]
        res["sigma_z"] += part["z"]
    return res if res["ions"] else None


def _data_volume(out_dir):
    try:
        L = open(os.path.join(out_dir, "system.data")).read().splitlines()
        b = {}
        for l in L[:40]:
            w = l.split()
            if len(w) >= 4 and w[2].endswith("lo") and w[3].endswith("hi"):
                b[w[2][0]] = float(w[1]) - float(w[0])
        return b["x"] * b["y"] * b["z"]
    except (OSError, KeyError, ValueError):
        return None


def ionic_lines(ionic):
    L = [f"イオン伝導度（ネルンスト–アインシュタイン：σ = e²/(V k_B T) Σ N z² D、イオン間の相関は含まない；T = {ionic['T']:.1f} K、"
         f"V = {ionic['V']:.0f} Å³）：全体 {ionic['sigma_total']:.4g} S/m（面内 xy {ionic['sigma_xy']:.4g}、垂直 z {ionic['sigma_z']:.4g}）"]
    for el, d in ionic["ions"].items():
        L.append(f"  {el}（{d['N']} 個、電荷 {d['q']:+.3f} e）：{d['sigma_total']:.4g} S/m（xy {d['sigma_xy']:.4g}）")
    return L

# --------------------------------------------------------------------------------------------- compute (fix ave/*)
def ave_files(out_dir):
    """[(fix ID, style, file)] of the averaging fixes with `file` in in.lammps whose file exists (msd_*.dat excluded)."""
    out = []
    try:
        lines = open(os.path.join(out_dir, "in.lammps")).read().splitlines()
    except OSError:
        return out
    for line in lines:
        w = line.split("#")[0].split()
        if len(w) > 4 and w[0] == "fix" and w[3] in AVE_FIXES and "file" in w[:-1]:
            f = w[w.index("file") + 1]
            if not os.path.basename(f).startswith("msd_") and os.path.isfile(os.path.join(out_dir, f)):
                out.append((w[1], w[3], f))
    return out


def read_ave(path):
    """a LAMMPS fix ave/* file -> ("scalar", column names, array) or ("blocks", row column names, [(step, array)]).
    Scalar files have 2 comment lines (title, TimeStep + values); block files (ave/time vector, ave/chunk, ave/histo,
    ave/correlate) have 3 (title, TimeStep + block header, row columns) and per output a header line and its rows."""
    heads, body = [], []
    for line in open(path):
        if line.startswith("#"):
            if not body:
                heads.append(line[1:].split())
        elif line.strip():
            body.append(line.split())
    names = heads[-1] if heads else []
    if len(heads) >= 3 and heads[-2] and heads[-2][0].lower() == "timestep":
        blocks, i = [], 0
        while i < len(body):
            step, n = float(body[i][0]), int(float(body[i][1]))
            rows = body[i + 1:i + 1 + n]
            if len(rows) < n:                                  # last block still being written
                break
            blocks.append((step, np.array(rows, dtype=float)))
            i += 1 + n
        return "blocks", names, blocks
    width = len(body[0]) if body else len(names)
    return "scalar", names, np.array([r for r in body if len(r) == width], dtype=float).reshape(-1, width)


def computes(out_dir, dest, skip_ps, dt_fs):
    """figures (one panel per column, as the columns have different units) and tables of the averaging fixes' files;
    a summary per file."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = []
    for fid, style, f in ave_files(out_dir):
        base = os.path.splitext(os.path.basename(f))[0]
        item = dict(fix=fid, style=style, file=f)
        try:
            kind, names, data = read_ave(os.path.join(out_dir, f))
        except Exception as e:
            item["error"] = str(e)
            res.append(item)
            continue
        if kind == "scalar" and len(data):
            x = data[:, 0] * dt_fs / 1000.0
            m = x >= skip_ps
            cols = [(k, names[k] if k < len(names) else f"col{k}") for k in range(1, data.shape[1])]
            ys, xlab = [data[:, k] for k, _ in cols], "time (ps)"
            item["mean_over_ps"] = [float(x[m][0]), float(x[-1])] if m.any() else None
            item["means"] = {n: float(data[m, k].mean()) for k, n in cols} if m.any() else {}
            item["rows"] = int(len(data))
            np.savetxt(os.path.join(dest, f"compute_{base}.csv"), np.c_[x, data[:, 1:]], delimiter=",",
                       header=",".join(["time_ps"] + [n for _, n in cols]), comments="", fmt="%.6g")
        elif kind == "blocks" and data:
            use = [r for s, r in data if s * dt_fs / 1000.0 >= skip_ps] or [data[-1][1]]
            use = [r for r in use if r.shape == use[-1].shape]
            avg = np.mean(np.array(use), axis=0)
            # x: Coord1 (ave/chunk), Coord (ave/histo), TimeDelta (ave/correlate), or the first column of an ave/time
            # vector when it increases (e.g. compute rdf: r); else the row number
            xk = 1 if len(names) > 1 and (names[1] in ("Coord1", "Coord", "TimeDelta") or
                                          (names[0] == "Row" and avg.shape[0] > 2 and np.all(np.diff(avg[:, 1]) > 0))) else 0
            cols = [(k, names[k] if k < len(names) else f"col{k}") for k in range(xk + 1, avg.shape[1])]
            cols = [(k, n) for k, n in cols if n not in ("Ncount", "OrigID")]
            x, ys, xlab = avg[:, xk], [avg[:, k] for k, _ in cols], names[xk] if xk < len(names) else "row"
            np.savetxt(os.path.join(dest, f"compute_{base}_average.csv"), avg, delimiter=",",
                       header=",".join(names[:avg.shape[1]]), comments="", fmt="%.6g")
            item["blocks_averaged"] = len(use)
            item["blocks_total"] = len(data)
        else:
            item["error"] = "no data yet"
            res.append(item)
            continue
        npan = min(len(cols), 6)                    # columns beyond the 6th go into the last panel
        fig, axs = plt.subplots(npan, 1, figsize=(8, 1.2 + 1.9 * npan), sharex=True, squeeze=False)
        for i, ((k, n), y) in enumerate(zip(cols, ys)):
            ax = axs[min(i, npan - 1), 0]
            ax.plot(x, y, lw=1.2, label=n, color=None if i >= npan - 1 and len(cols) > npan else "C0")
            ax.set_ylabel(n if len(cols) <= npan or i < npan - 1 else "", fontsize=9)
            if kind == "scalar" and m.any():
                ax.axvspan(x[m][0], x[-1], color="#999", alpha=0.07, lw=0)
                ax.axhline(item["means"][n], color="#d62728", lw=0.8, ls="--")
        if len(cols) > npan:
            axs[-1, 0].legend(fontsize=7)
        axs[-1, 0].set_xlabel(xlab)
        axs[0, 0].set_title(f"fix {fid} {style} -> {f}" + (f" (mean of {item['blocks_averaged']} outputs)" if "blocks_averaged" in item
                                                         else " (dashed: mean over the shaded part)"), loc="left", fontsize=10)
        fig.tight_layout()
        item["png"] = f"compute_{base}.png"
        fig.savefig(os.path.join(dest, item["png"]), dpi=100)
        plt.close(fig)
        res.append(item)
    return res


# ------------------------------------------------------------------------------------------------------ XRD (Debyer)
def find_debyer(path=None):
    """the debyer executable: the given path, $DEBYER, debyer on PATH, or a source build ~/debyer/debyer/debyer or
    ~/*/debyer/debyer/debyer; None when not found."""
    for c in (path, os.environ.get("DEBYER"), shutil.which("debyer")):
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    home = os.path.expanduser("~")
    for c in [os.path.join(home, "debyer", "debyer", "debyer")] + sorted(glob.glob(os.path.join(home, "*", "debyer", "debyer", "debyer"))):
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None


def xrd(out_dir, dest, skip_ps, dt_fs, debyer=None, wavelength=1.5406, two_theta=(3.0, 80.0, 0.02), frames=3,
        mode="cluster", size_A=50.0, log=print, traj="traj.lammpstrj"):
    """X-ray powder pattern I(2 theta) by Debyer, averaged over the last `frames` frames of traj.lammpstrj.
    mode "cluster" (crystals, layered minerals): the cell repeated n_a x n_b x n_c times so that every edge is at least
        size_A, without periodicity (Debye formula of a crystallite; peaks broadened by its size, small-angle scattering
        of the particle below a few degrees);
    mode "periodic" (liquids, glasses; orthorhombic box only): the box with periodic boundaries, cut-off half the
        shortest box length, sinc damping and the number density of the box (no reflection with d > half the box)."""
    from .sova_structure import frames as read_frames
    exe = find_debyer(debyer)
    if not exe:
        raise RuntimeError("debyer not found: build https://github.com/wojdyr/debyer and give its path "
                           "(--debyer / $DEBYER, or the GUI settings)")
    el = [r["element"] for r in csv.DictReader(open(os.path.join(out_dir, "charges.csv")))]
    fr = [f for f in read_frames(os.path.join(out_dir, traj)) if f[0] * dt_fs / 1000.0 >= skip_ps - 1e-9]
    if not fr:
        raise RuntimeError(f"no frames after the equilibration in {traj}")
    fr = fr[-max(1, int(frames)):]
    ortho = all(abs(c[i][j]) < 1e-6 for _, c, _ in fr for i in range(3) for j in range(3) if i != j)
    if mode == "periodic" and not ortho:
        raise RuntimeError("the periodic mode of Debyer needs an orthorhombic box (this box is tilted): use the cluster mode")
    reps = [max(1, int(np.ceil(size_A / np.linalg.norm(fr[-1][1][i]) - 1e-9))) for i in range(3)] if mode == "cluster" else [1, 1, 1]
    t0, t1, step = two_theta
    x, pats = None, []
    tmp = tempfile.mkdtemp(prefix=".debyer-", dir=dest)
    try:
        for k, (st, cell, pos) in enumerate(fr):
            if mode == "periodic":
                L = np.diag(cell)
                opts = ["-a", f"{L[0]:.6f}", "-b", f"{L[1]:.6f}", "-c", f"{L[2]:.6f}", "-r", f"{0.5 * L.min():.4f}",
                        "--sinc", "--ro", f"{len(pos) / float(np.prod(L)):.6f}"]
                syms, xyz = el, pos
            else:
                shift = np.array([[i, j, m] for i in range(reps[0]) for j in range(reps[1]) for m in range(reps[2])], float) @ cell
                xyz = (pos[None, :, :] + shift[:, None, :]).reshape(-1, 3)
                syms, opts = el * len(shift), []
            fxyz, fout = os.path.join(tmp, f"frame{k}.xyz"), os.path.join(tmp, f"frame{k}.dat")
            with open(fxyz, "w") as f:
                f.write(f"{len(xyz)}\nstep {st}\n")
                f.writelines(f"{s} {a:.5f} {b:.5f} {c:.5f}\n" for s, (a, b, c) in zip(syms, xyz))
            cmd = [exe, "-x", "-l", f"{wavelength}", "-f", f"{t0}", "-t", f"{t1}", "-s", f"{step}"] + opts + ["-o", fout, fxyz]
            log(f"debyer: frame {k + 1}/{len(fr)} (step {st}, {len(xyz)} atoms"
                + (f", cell x {reps[0]}x{reps[1]}x{reps[2]})" if mode == "cluster" else ", periodic)"))
            p = subprocess.run(cmd, capture_output=True, text=True)
            if p.returncode != 0 or not os.path.exists(fout):
                raise RuntimeError("debyer failed: " + (p.stderr or p.stdout)[-800:])
            d = np.loadtxt(fout, comments="#", ndmin=2)
            x = d[:, 0]
            pats.append(d[:, 1])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    inten = np.mean(pats, axis=0)
    size = [float(np.linalg.norm(fr[-1][1][i]) * reps[i]) for i in range(3)]
    # a maximum with d above a third of the crystallite (fewer than 3 lattice planes across it) or above the cut-off of
    # the periodic mode is not a Bragg reflection but an oscillation of the crystallite's size (or of the cut-off)
    d_lim = max(size) / 3.0 if mode == "cluster" else 0.5 * min(np.diag(fr[-1][1]))
    dsp = wavelength / (2 * np.sin(np.radians(x / 2)))
    loc = [k for k in range(1, len(inten) - 1) if inten[k] >= inten[k - 1] and inten[k] > inten[k + 1]]
    bragg = [k for k in loc if dsp[k] <= d_lim]
    top = max([inten[k] for k in bragg] or [inten.max()]) or 1.0       # strongest reflection = 100 (not the small-angle edge)
    with open(os.path.join(dest, "xrd.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["two_theta_deg", "intensity", "relative_intensity"])
        for a, b in zip(x, inten):
            w.writerow([f"{a:.4f}", f"{b:.6g}", f"{100 * b / top:.3f}"])
    peaks = [dict(two_theta=float(x[k]), d_A=float(dsp[k]), relative=float(100 * inten[k] / top))
             for k in sorted(sorted(bragg, key=lambda k: -inten[k])[:8])]
    size_osc = [dict(two_theta=float(x[k]), d_A=float(dsp[k]), relative=float(100 * inten[k] / top))
                for k in loc if dsp[k] > d_lim]
    return dict(debyer=exe, wavelength_A=float(wavelength), two_theta=[float(t0), float(t1), float(step)], frames=len(fr),
                steps=[int(f[0]) for f in fr], mode=mode, periodic=mode == "periodic", repeat=reps, size_A=size,
                atoms=len(el) * int(np.prod(reps)), d_limit_A=float(d_lim), peaks=peaks, size_oscillations=size_osc)


# ------------------------------------------------------------------------------------------------------------ main
def _time_step_fs(out_dir, default):
    """the time step of the in.lammps that was run (its last `timestep` command), else the default."""
    dt = default
    try:
        for line in open(os.path.join(out_dir, "in.lammps")):
            w = line.split("#")[0].split()
            if len(w) == 2 and w[0] == "timestep":
                try:
                    dt = float(w[1])
                except ValueError:
                    pass
    except OSError:
        pass
    return dt


def trajectories(out_dir):
    """the LAMMPS dumps of an output folder (traj.lammpstrj first, then the others)."""
    fs = sorted(os.path.basename(f) for f in glob.glob(os.path.join(out_dir, "*.lammpstrj"))
                if os.path.basename(f) != "swelling_states.lammpstrj")      # one frame per water content: not an MD
    return sorted(fs, key=lambda f: f != "traj.lammpstrj")


def analyze(out_dir, skip_ps=None, sova_python=None, sova_home=None, log=print, what=WHAT, debyer=None, xrd_opts=None,
            traj=None):
    """traj: the dump used by SOVA and XRD (default traj.lammpstrj, else the first *.lammpstrj of the folder)."""
    out_dir = os.path.abspath(out_dir)
    rep = json.load(open(os.path.join(out_dir, "report.json"))) if os.path.exists(os.path.join(out_dir, "report.json")) else {}
    dt_fs = _time_step_fs(out_dir, float(rep.get("md", {}).get("timestep_fs", 1.0)))
    dest = os.path.join(out_dir, "analysis")
    os.makedirs(dest, exist_ok=True)
    for old in ("thermo.csv", "density_z.csv", "msd.csv", "rdf.csv", "sova_partial_gr.csv"):
        if os.path.exists(os.path.join(dest, old)):                 # outputs of earlier versions
            os.remove(os.path.join(dest, old))
    traj = traj or (trajectories(out_dir) or ["traj.lammpstrj"])[0]
    has_traj = os.path.exists(os.path.join(out_dir, traj))
    msd_files = glob.glob(os.path.join(out_dir, "msd_*.dat"))
    t_end = 0.0
    for f in msd_files:
        a = np.loadtxt(f, comments="#", ndmin=2)
        if len(a):
            t_end = max(t_end, a[-1, 0] * dt_fs / 1000.0)
    if not msd_files and has_traj:
        with open(os.path.join(out_dir, traj)) as f:                 # time of the last frame
            for line in f:
                if line.startswith("ITEM: TIMESTEP"):
                    t_end = int(next(f)) * dt_fs / 1000.0
    if skip_ps is None:
        skip_ps = round(0.2 * t_end, 3)
    log(f"MD 0-{t_end:.1f} ps; equilibration skipped: first {skip_ps} ps")
    summary = dict(out_dir=out_dir, t_end_ps=float(t_end), skip_ps=float(skip_ps), timestep_fs=dt_fs, what=list(what), traj=traj,
                   has_mineral=any(r.get("role") == "mineral" for r in csv.DictReader(open(os.path.join(out_dir, "charges.csv")))),
                   tools=dict(msd="LAMMPS compute msd (Einstein relation)", pdf_rdf="SOVA (sovapy) pair_dist_func, radial_dist_fun",
                              sq="SOVA (sovapy) partial_structure_factor, structure_factor_xray",
                              compute="LAMMPS fix ave/time, ave/chunk, ave/correlate, ave/histo output files",
                              xrd="Debyer (Debye scattering equation)",
                              swelling="LAMMPS fix npt; d001 = lz of the cell (swelling.dat of the soapff swelling block)"),
                   figures=[])
    dres, dcurves = {}, {}
    if "msd" in what:
        dres, dcurves = msd(out_dir, skip_ps, dt_fs)
        summary["diffusion_1e-9_m2_s"] = dres
        if not msd_files:
            summary["diffusion_note"] = "no LAMMPS msd_*.dat in this folder: add the MSD of soapff (catalog) before run"
        with open(os.path.join(dest, "msd_fit.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["group", "D_total_1e-9_m2_s", "D_xy", "D_z", "fit_from_ps", "fit_to_ps"])
            for g, d in dres.items():
                w.writerow([g, d["D_total"], d["D_xy"], d["D_z"], d["fit_from_ps"], d["fit_to_ps"]])

    if "sova" in what:
        try:
            if not has_traj:
                raise RuntimeError(f"no {traj} in this folder")
            summary["sova"] = sova(out_dir, dest, skip_ps, dt_fs, sova_python, sova_home, log, traj=traj)
        except Exception as e:
            summary["sova"] = dict(error=str(e))
            log(f"WARNING {e}")
    if "msd" in what or "sova" in what:
        plot(dest, dcurves, dres, sova="sova" in what and "error" not in summary["sova"])
        summary["figures"].append(["MSD・PDF・RDF", "analysis.png"])
        if "sova" in what and summary["sova"].get("sq_xray"):
            plot_sq(dest, summary["sova"].get("sq_rmax"))
            summary["figures"].append(["S(Q)", "sova_sq.png"])
    if "compute" in what:
        summary["computes"] = computes(out_dir, dest, skip_ps, dt_fs)
        summary["figures"] += [[f"compute：{c['file']}", c["png"]] for c in summary["computes"] if c.get("png")]
    if "xrd" in what:
        try:
            if not has_traj:
                raise RuntimeError(f"no {traj} in this folder")
            summary["xrd"] = xrd(out_dir, dest, skip_ps, dt_fs, debyer, log=log, traj=traj, **(xrd_opts or {}))
            plot_xrd(dest, summary["xrd"])
            summary["figures"].append(["XRD", "xrd.png"])
        except Exception as e:
            summary["xrd"] = dict(error=str(e))
            log(f"WARNING {e}")
    if "swelling" in what:
        try:
            summary["swelling"] = swelling(out_dir, dest)
            summary["figures"].append(["膨潤曲線", "swelling.png"])
        except FileNotFoundError:
            summary["swelling"] = dict(error="no swelling.dat in this folder")
        except Exception as e:
            summary["swelling"] = dict(error=str(e))
            log(f"WARNING swelling: {e}")
    json.dump(summary, open(os.path.join(dest, "summary.json"), "w"), indent=1, ensure_ascii=False)
    txt = summary_text(summary)
    open(os.path.join(dest, "summary.txt"), "w").write(txt)
    log(txt)
    return summary


def summary_text(s):
    what = s.get("what", ("msd", "sova"))
    L = [f"解析：MD 0–{s['t_end_ps']:.1f} ps、最初の {s['skip_ps']:.1f} ps を平衡化として除外"
         + (f"（時間刻み {s['timestep_fs']} fs）" if "timestep_fs" in s else "")
         + (f"。SOVA と XRD の軌道：{s['traj']}" if s.get("traj") and ("sova" in s.get("what", ()) or "xrd" in s.get("what", ())) else "")]
    if "msd" in what:
        if s.get("diffusion_1e-9_m2_s"):
            L.append("拡散係数（10⁻⁹ m²/s、LAMMPS compute msd の傾き、アインシュタインの式、有限サイズ補正なし）")
            for g, d in s["diffusion_1e-9_m2_s"].items():
                L.append(f"  {('水（O）' if g == 'Ow' else g):8s} D = {d['D_total']:.3f}（面内 xy {d['D_xy']:.3f}、垂直 z {d['D_z']:.3f}）  "
                         f"{d['fit_from_ps']:.1f}–{d['fit_to_ps']:.1f} ps で直線近似")
            if min(d["fit_to_ps"] - d["fit_from_ps"] for d in s["diffusion_1e-9_m2_s"].values()) < 5.0:
                L.append("  ⚠ 直線近似の範囲が 5 ps に満たないので、拡散係数は信頼できません。MD を長く（目安 50 ps 以上）してください")
            if s.get("has_mineral"):
                L.append("  層間や表面の近くでは z 方向の動きが制限されるので、拡散係数は面内 xy の値を使う")
            if s.get("ionic_conductivity"):
                L += ionic_lines(s["ionic_conductivity"])
        elif s.get("diffusion_note"):
            L.append("拡散係数：msd_*.dat がありません。カタログの「拡散係数の MSD」（soapff）を run の前に入れてください")
    if "sova" in what:
        st = s.get("sova", {})
        if st.get("pairs"):
            L.append(f"PDF g(r) と RDF N(r) = 4πr²ρg(r)（SOVA {st.get('sovapy_version', '')}、r ≤ {st['rmax']:.1f} Å、"
                     f"{st['frames']} フレーム）：sova_pdf.csv、sova_rdf.csv")
            for p, d in st.get("first_peaks", {}).items():
                cn = st.get("coordination", {}).get(p)
                L.append(f"  {p.replace('Ow', '水O').replace('O_min', '鉱物O'):12s} 第 1 ピーク {d['r_peak']:.2f} Å（g = {d['g_peak']:.2f}）"
                         + (f"、配位数 {cn['cn']:.2f}（{cn['r_min']:.2f} Å の極小まで）" if cn else ""))
            n_all = sum(1 for p in st["pairs"] if p.startswith("all:"))
            if n_all:
                L.append(f"  ほかに全原子の元素の組 {n_all} 通り（all:）")
            if st.get("sq_xray"):
                L.append(f"X 線の構造因子 S(Q)（SOVA、全原子、g(r) を {st.get('sq_rmax', 0):.1f} Å まで使用、"
                         f"Q ≤ {st.get('sq_qmax', 0):.0f} Å⁻¹）：sova_sq.csv")
                if st.get("sq_rmax"):
                    L.append(f"  Q < 2π/{st['sq_rmax']:.1f} Å = {2 * np.pi / st['sq_rmax']:.2f} Å⁻¹ は g(r) の打ち切りで決まらない"
                             + ("。層状・界面の系では g(r) が 1 に近づかないので、低 Q の値は使わない" if s.get("has_mineral") else ""))
        elif st.get("error"):
            L.append("PDF と RDF（SOVA）：実行できませんでした。" + st["error"].splitlines()[0])
    if "swelling" in what and s.get("swelling", {}).get("states"):
        sw = s["swelling"]
        L.append(f"膨潤曲線（d001 = lz、各状態の NPT の後半 {sw['avg_ps']:.1f} ps の平均と標準偏差；粘土 {sw['clay_mass_amu']:.0f} g/mol）：swelling.csv")
        for r in sw["states"]:
            if r["direction"] == "up" and (r["waters"] <= 3 or r["waters"] % 10 == 0
                                           or r["waters"] == max(u["waters"] for u in sw["states"])):
                L.append(f"  水 {r['waters']:3d} 個（{r['water_g_per_g']:.3f} g/g）：d001 = {r['d001']:.2f} ± {r['d001_sd']:.2f} Å")
        if any(r["direction"] == "down" for r in sw["states"]):
            dd = [abs(r["d001"] - u["d001"]) for r in sw["states"] if r["direction"] == "down"
                  for u in sw["states"] if u["direction"] == "up" and u["waters"] == r["waters"]]
            if dd:
                L.append(f"  脱水（水を減らす向き）と膨潤の差：平均 {np.mean(dd):.2f} Å、最大 {np.max(dd):.2f} Å（同じ水の量の {len(dd)} 状態）")
    elif "swelling" in what and s.get("swelling", {}).get("error") and "no swelling.dat" not in s["swelling"]["error"]:
        L.append("膨潤曲線：解析できませんでした。" + s["swelling"]["error"])
    if "compute" in what:
        cs = s.get("computes") or []
        if not cs:
            L.append("compute の結果：in.lammps に file 付きの fix ave/time・ave/chunk・ave/correlate・ave/histo がありません")
        else:
            L.append(f"compute の結果（in.lammps の fix ave/* が書いたファイル {len(cs)} 個）")
        for c in cs:
            head = f"  {c['file']}（fix {c['fix']} {c['style']}）"
            if c.get("error"):
                L.append(head + f"：読めませんでした（{c['error']}）")
            elif "means" in c:
                L.append(head + (f"：{c['mean_over_ps'][0]:.1f}–{c['mean_over_ps'][1]:.1f} ps の平均 " if c.get("mean_over_ps") else "：")
                         + "、".join(f"{k} = {v:.5g}" for k, v in c["means"].items()))
            else:
                L.append(head + f"：{c['blocks_total']} 回の出力のうち {c['blocks_averaged']} 回を平均（compute_*_average.csv）")
    if "xrd" in what:
        x = s.get("xrd", {})
        if x.get("peaks") is not None:
            L.append(f"XRD（debyer、λ = {x['wavelength_A']} Å、2θ = {x['two_theta'][0]:g}–{x['two_theta'][1]:g}°、"
                     f"最後の {x['frames']} フレームの平均、"
                     + ("周期境界（カットオフ = 箱の半分、sinc 減衰）" if x["periodic"] else
                        "セルを {}×{}×{} 並べた結晶子（{:.0f}×{:.0f}×{:.0f} Å、{} 原子）".format(*x["repeat"], *x["size_A"], x["atoms"]))
                     + "）：xrd.csv。強いピーク（最も強いピーク = 100）")
            for pk in x["peaks"]:
                L.append(f"  2θ = {pk['two_theta']:6.2f}°   d = {pk['d_A']:6.3f} Å   相対強度 {pk['relative']:5.1f}")
            if x.get("size_oscillations"):
                L.append(f"  d > {x['d_limit_A']:.1f} Å の極大は" + ("結晶子の大きさ" if not x["periodic"] else "カットオフ")
                         + "による振動（ブラッグ反射ではない）：" + "、".join(f"{p['two_theta']:.2f}°" for p in x["size_oscillations"]))
            if not x["periodic"]:
                L.append("  結晶子が小さいほどピークは広がり、低角のピーク（001 など）は低角側へずれる。--xrd-size（GUI では結晶子の大きさ）で大きくできる（計算時間は原子数の 2 乗）")
        elif x.get("error"):
            L.append("XRD（debyer）：実行できませんでした。" + x["error"].splitlines()[0])
    return "\n".join(L) + "\n"


def plot(dest, dcurves, dres, sova=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    cols = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b", "#7f7f7f"]
    a = ax[0]
    for (g, (tt, m, (f0, f1))), c in zip(dcurves.items(), cols):
        lab = "water O" if g == "Ow" else g
        a.plot(tt, m[:, 3], color=c, lw=1.5, label=f"{lab}: D = {dres[g]['D_total']:.2f}")
        a.plot(tt, m[:, 0] + m[:, 1], color=c, lw=1.0, ls="--", label=f"{lab} xy: D = {dres[g]['D_xy']:.2f}")
        a.axvspan(f0, f1, color=c, alpha=0.05, lw=0)
    if dcurves:
        a.legend(fontsize=8, title="D (1e-9 m2/s)", title_fontsize=8)
    a.set_xlabel("time (ps)"); a.set_ylabel("MSD (A$^2$)"); a.set_title("(a) LAMMPS compute msd", loc="left", fontsize=10)
    for a, fname, ylab, title in ((ax[1], "sova_pdf.csv", "g(r)", "(b) SOVA PDF g(r)"),
                                  (ax[2], "sova_rdf.csv", "N(r) = 4$\\pi$r$^2\\rho$g(r) (A$^{-1}$)", "(c) SOVA RDF")):
        p = os.path.join(dest, fname)
        if sova and os.path.exists(p):
            rows = list(csv.reader(open(p)))
            head, data = rows[0], np.array(rows[1:], dtype=float)
            for n, c in zip([h for h in head[1:] if not h.startswith("all:")], cols):
                a.plot(data[:, 0], data[:, head.index(n)], color=c, lw=1.3, label=n.replace("O_min", "O(mineral)"))
            if fname == "sova_pdf.csv":
                a.axhline(1, color="#999", lw=0.8, ls=":")
            a.legend(fontsize=8)
        a.set_xlabel("r (A)"); a.set_ylabel(ylab); a.set_title(title, loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(dest, "analysis.png"), dpi=100); plt.close(fig)


def plot_sq(dest, rmax=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = np.loadtxt(os.path.join(dest, "sova_sq.csv"), delimiter=",", skiprows=1, ndmin=2)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(d[:, 0], d[:, 1], lw=1.3, color="#1f77b4")
    ax.axhline(1, color="#999", lw=0.8, ls=":")
    if rmax:                    # Q < 2 pi / r_max is not resolved by g(r) up to r_max (truncation, inhomogeneous box)
        q0 = 2 * np.pi / rmax
        ax.axvspan(0, q0, color="#999", alpha=0.15, lw=0)
        ok = d[:, 0] >= q0
        if ok.any():
            lo, hi = d[ok, 1].min(), d[ok, 1].max()
            ax.set_ylim(lo - 0.15 * (hi - lo), hi + 0.15 * (hi - lo))
        ax.text(q0, 0.97, f" Q < 2π/{rmax:.1f} Å: not resolved", transform=ax.get_xaxis_transform(), fontsize=8,
                va="top", color="#666")
    ax.set_xlabel("Q (A$^{-1}$)"); ax.set_ylabel("S(Q)")
    ax.set_title("SOVA X-ray structure factor S(Q) of the whole system", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(dest, "sova_sq.png"), dpi=100); plt.close(fig)


def swelling(out_dir, dest):
    """swelling.dat (step added_waters phase lx ly lz xy xz yz pe temp, every 100 fs) -> one row per state (the
    consecutive lines with the same number of added waters; phase 1 = averaged): water molecules, g H2O / g clay, d001
    (lz) mean and standard deviation, cell a b c alpha beta gamma, energy; 'up' until the largest water content, then
    'down'.  Writes swelling.csv and swelling.png (with swelling_exp.csv if present)."""
    a = np.loadtxt(os.path.join(out_dir, "swelling.dat"), comments="#", ndmin=2)
    if not len(a):
        raise ValueError("swelling.dat is empty")
    rows = list(csv.DictReader(open(os.path.join(out_dir, "charges.csv"))))
    from ._clayff.chem import MASS
    clay_mass = sum(MASS[r["element"]] for r in rows if r["clayff_type"] not in ("ow", "hw"))
    w0 = sum(1 for r in rows if r["clayff_type"] == "ow")
    m_w = 2 * MASS["H"] + MASS["O"]
    states, start = [], 0
    nmax = a[:, 1].max()
    reached = False
    for k in range(1, len(a) + 1):
        if k == len(a) or a[k, 1] != a[start, 1]:
            blk = a[start:k]
            prod = blk[blk[:, 2] == 1]
            if len(prod):
                lx, ly, lz, xy, xz, yz = (prod[:, j] for j in range(3, 9))
                A = lx; B = np.sqrt(ly ** 2 + xy ** 2); C = np.sqrt(lz ** 2 + xz ** 2 + yz ** 2)
                al = np.degrees(np.arccos((xy * xz + ly * yz) / (B * C)))
                be = np.degrees(np.arccos(xz / C)); ga = np.degrees(np.arccos(xy / B))
                n = int(round(blk[0, 1])) + w0
                states.append(dict(waters=n, water_g_per_g=n * m_w / clay_mass, direction="down" if reached else "up",
                                   d001=float(lz.mean()), d001_sd=float(lz.std()), a=float(A.mean()), b=float(B.mean()),
                                   c=float(C.mean()), alpha=float(al.mean()), beta=float(be.mean()), gamma=float(ga.mean()),
                                   pe=float(prod[:, 9].mean()), T=float(prod[:, 10].mean()), frames=int(len(prod))))
            if a[start, 1] == nmax:
                reached = True
            start = k
    dt = np.median(np.diff(a[:, 0])) if len(a) > 1 else 1
    avg_steps = np.median([s_["frames"] for s_ in states]) * dt
    with open(os.path.join(dest, "swelling.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(states[0]))
        w.writeheader(); w.writerows(states)
    exp = []
    ef = os.path.join(out_dir, "swelling_exp.csv")
    if os.path.exists(ef):
        exp = list(csv.DictReader(open(ef)))
    plot_swelling(dest, states, exp, clay_mass)
    dt_fs = _time_step_fs(out_dir, 0.5)
    return dict(states=states, clay_mass_amu=clay_mass, water0=w0, avg_ps=float(avg_steps * dt_fs / 1000.0),
                experiment=os.path.basename(ef) if exp else None)


def plot_swelling(dest, states, exp, clay_mass):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    up = [r for r in states if r["direction"] == "up"]
    dn = [r for r in states if r["direction"] == "down"]
    ax.plot([r["water_g_per_g"] for r in up], [r["d001"] for r in up], "o-", color="#c0392b", ms=4, lw=1, label="MD (swelling)")
    if dn:
        ax.plot([r["water_g_per_g"] for r in dn], [r["d001"] for r in dn], "s--", color="#2c6fbb", mfc="none", ms=4, lw=1,
                label="MD (deswelling)")
    if exp:
        labs = sorted({e.get("label") or "experiment" for e in exp})
        for lab in labs:
            pts = [e for e in exp if (e.get("label") or "experiment") == lab]
            ax.plot([float(e["water_g_per_g"]) for e in pts], [float(e["d001_A"]) for e in pts], "^", color="#2a9d4b", ms=6,
                    label=lab)
    ax.set_xlabel("water content (g H$_2$O / g clay)")
    ax.set_ylabel("d$_{001}$ (Å)")
    sec = ax.secondary_xaxis("top", functions=(lambda x: x * clay_mass / 18.0153, lambda n: n * 18.0153 / clay_mass))
    sec.set_xlabel("H$_2$O per cell")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(dest, "swelling.png"), dpi=200)
    plt.close(fig)


def plot_xrd(dest, x):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = np.loadtxt(os.path.join(dest, "xrd.csv"), delimiter=",", skiprows=1, ndmin=2)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(d[:, 0], np.minimum(d[:, 2], 111), lw=1.0, color="#1f77b4")
    for pk in x["peaks"]:
        ax.annotate(f"{pk['d_A']:.2f} A", (pk["two_theta"], min(pk["relative"], 104)), xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=8)
    if x.get("d_limit_A") and x["d_limit_A"] * 2 > x["wavelength_A"]:
        tl = 2 * np.degrees(np.arcsin(x["wavelength_A"] / (2 * x["d_limit_A"])))
        if tl > d[0, 0]:
            ax.axvspan(d[0, 0], tl, color="#999", alpha=0.15, lw=0)
            ax.text(tl, 0.97, f" d > {x['d_limit_A']:.0f} A: size oscillation", transform=ax.get_xaxis_transform(),
                    fontsize=8, va="top", color="#666")
    ax.set_xlim(d[0, 0], d[-1, 0]); ax.set_ylim(0, 112)
    ax.set_xlabel(f"2$\\theta$ (deg), $\\lambda$ = {x['wavelength_A']} A"); ax.set_ylabel("relative intensity")
    ax.set_title(f"Debyer X-ray powder pattern: mean of the last {x['frames']} frames, "
                 + ("periodic box" if x["periodic"] else "crystallite {}x{}x{} cells".format(*x["repeat"])), loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(dest, "xrd.png"), dpi=100); plt.close(fig)
