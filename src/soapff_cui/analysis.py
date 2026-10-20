"""Analyses of a soapff-cui MD run, computed by established tools (no analysis algorithm of soapff-cui itself):

  msd       LAMMPS `compute msd` (in in.lammps, written every 0.1 ps to msd_<group>.dat) of water O and each ion type;
            self-diffusion coefficient from the Einstein relation, D = slope / (2 d), slope of a least-squares line
            over the production part: total (d = 3), in plane (xy, d = 2), normal (z, d = 1).  No finite-size correction.
  PDF, RDF  SOVA (sovapy, soapff_cui/sova_structure.py): pair distribution functions g(r) (atom_pair_hist,
            pair_dist_func) and radial distribution functions N(r) = 4 pi r^2 rho g(r) (radial_dist_fun) of water O -
            water O, ion - water O, ion - mineral O and every element pair of the whole system

  python -m soapff_cui analyze OUTDIR [--skip-ps X] [--sova-python PY] [--sova-home DIR]
      -> OUTDIR/analysis/{summary.txt, summary.json, analysis.png, msd_fit.csv, sova_pdf.csv, sova_rdf.csv,
                          sova_summary.json}

SOVA: used directly when sovapy is importable; otherwise sova_structure.py is run with the Python of a SOVA installation
(--sova-python or $SOVA_PYTHON; --sova-home or $SOVA_HOME is put on its PYTHONPATH for a source checkout).
"""
import csv
import glob
import json
import os
import subprocess
import sys

import numpy as np

A2PS_TO_1E9M2S = 10.0          # 1 A^2/ps = 1e-8 m^2/s = 10 x 1e-9 m^2/s


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


def sova(out_dir, dest, skip_ps, dt_fs, sova_python=None, sova_home=None, log=print):
    """run sova_structure.py (SOVA PDF and RDF) in this Python if sovapy is importable, else with the given SOVA Python."""
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sova_structure.py")
    args = ["--traj", os.path.join(out_dir, "traj.lammpstrj"), "--charges", os.path.join(out_dir, "charges.csv"),
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


def analyze(out_dir, skip_ps=None, sova_python=None, sova_home=None, log=print):
    out_dir = os.path.abspath(out_dir)
    rep = json.load(open(os.path.join(out_dir, "report.json"))) if os.path.exists(os.path.join(out_dir, "report.json")) else {}
    dt_fs = float(rep.get("md", {}).get("timestep_fs", 1.0))
    dest = os.path.join(out_dir, "analysis")
    os.makedirs(dest, exist_ok=True)
    for old in ("thermo.csv", "density_z.csv", "msd.csv", "rdf.csv", "sova_sq.csv", "sova_partial_gr.csv"):
        if os.path.exists(os.path.join(dest, old)):                 # outputs of earlier versions
            os.remove(os.path.join(dest, old))
    msd_files = glob.glob(os.path.join(out_dir, "msd_*.dat"))
    t_end = 0.0
    for f in msd_files:
        a = np.loadtxt(f, comments="#", ndmin=2)
        if len(a):
            t_end = max(t_end, a[-1, 0] * dt_fs / 1000.0)
    if not msd_files:
        with open(os.path.join(out_dir, "traj.lammpstrj")) as f:     # time of the last frame
            for line in f:
                if line.startswith("ITEM: TIMESTEP"):
                    t_end = int(next(f)) * dt_fs / 1000.0
    if skip_ps is None:
        skip_ps = round(0.2 * t_end, 3)
    log(f"MD 0-{t_end:.1f} ps; equilibration skipped: first {skip_ps} ps")
    summary = dict(out_dir=out_dir, t_end_ps=float(t_end), skip_ps=float(skip_ps),
                   has_mineral=any(r.get("role") == "mineral" for r in csv.DictReader(open(os.path.join(out_dir, "charges.csv")))),
                   tools=dict(msd="LAMMPS compute msd (Einstein relation)", pdf_rdf="SOVA (sovapy) pair_dist_func, radial_dist_fun"))
    dres, dcurves = msd(out_dir, skip_ps, dt_fs)
    summary["diffusion_1e-9_m2_s"] = dres
    if not msd_files:
        summary["diffusion_note"] = "no LAMMPS msd_*.dat in this folder (input written by an older soapff-cui): rerun make"
    with open(os.path.join(dest, "msd_fit.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["group", "D_total_1e-9_m2_s", "D_xy", "D_z", "fit_from_ps", "fit_to_ps"])
        for g, d in dres.items():
            w.writerow([g, d["D_total"], d["D_xy"], d["D_z"], d["fit_from_ps"], d["fit_to_ps"]])
    try:
        summary["sova"] = sova(out_dir, dest, skip_ps, dt_fs, sova_python, sova_home, log)
    except Exception as e:
        summary["sova"] = dict(error=str(e))
        log(f"WARNING {e}")
    json.dump(summary, open(os.path.join(dest, "summary.json"), "w"), indent=1)
    txt = summary_text(summary)
    open(os.path.join(dest, "summary.txt"), "w").write(txt)
    plot(dest, dcurves, dres)
    log(txt)
    return summary


def summary_text(s):
    L = [f"解析：MD 0–{s['t_end_ps']:.1f} ps、最初の {s['skip_ps']:.1f} ps を平衡化として除外"]
    if s["diffusion_1e-9_m2_s"]:
        L.append("拡散係数（10⁻⁹ m²/s、LAMMPS compute msd の傾き、アインシュタインの式、有限サイズ補正なし）")
        for g, d in s["diffusion_1e-9_m2_s"].items():
            L.append(f"  {('水（O）' if g == 'Ow' else g):8s} D = {d['D_total']:.3f}（面内 xy {d['D_xy']:.3f}、垂直 z {d['D_z']:.3f}）  "
                     f"{d['fit_from_ps']:.1f}–{d['fit_to_ps']:.1f} ps で直線近似")
        if min(d["fit_to_ps"] - d["fit_from_ps"] for d in s["diffusion_1e-9_m2_s"].values()) < 5.0:
            L.append("  ⚠ 直線近似の範囲が 5 ps に満たないので、拡散係数は信頼できません。MD を長く（目安 50 ps 以上）してください")
        if s.get("has_mineral"):
            L.append("  層間や表面の近くでは z 方向の動きが制限されるので、拡散係数は面内 xy の値を使う")
    elif s.get("diffusion_note"):
        L.append("拡散係数：この出力には LAMMPS の MSD がありません（古い版の入力）。make からやり直してください")
    st = s.get("sova", {})
    if st.get("pairs"):
        L.append(f"PDF g(r) と RDF N(r) = 4πr²ρg(r)（SOVA {st.get('sovapy_version', '')}、r ≤ {st['rmax']:.1f} Å、"
                 f"{st['frames']} フレーム）：sova_pdf.csv、sova_rdf.csv")
        for p, d in st.get("first_peaks", {}).items():
            L.append(f"  {p.replace('Ow', '水O').replace('O_min', '鉱物O'):12s} 第 1 ピーク {d['r_peak']:.2f} Å（g = {d['g_peak']:.2f}）")
        n_all = sum(1 for p in st["pairs"] if p.startswith("all:"))
        if n_all:
            L.append(f"  ほかに全原子の元素の組 {n_all} 通り（all:）")
    elif st.get("error"):
        L.append("PDF と RDF（SOVA）：実行できませんでした。" + st["error"].splitlines()[0])
    return "\n".join(L) + "\n"


def plot(dest, dcurves, dres):
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
        if os.path.exists(p):
            rows = list(csv.reader(open(p)))
            head, data = rows[0], np.array(rows[1:], dtype=float)
            for n, c in zip([h for h in head[1:] if not h.startswith("all:")], cols):
                a.plot(data[:, 0], data[:, head.index(n)], color=c, lw=1.3, label=n.replace("O_min", "O(mineral)"))
            if fname == "sova_pdf.csv":
                a.axhline(1, color="#999", lw=0.8, ls=":")
            a.legend(fontsize=8)
        a.set_xlabel("r (A)"); a.set_ylabel(ylab); a.set_title(title, loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(dest, "analysis.png"), dpi=100); plt.close(fig)
