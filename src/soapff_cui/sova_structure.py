#!/usr/bin/env python3
"""PDF and RDF of a soapff-cui MD run with SOVA (sovapy: Structural Order Visualization and Analysis, M. Shiga and
H. Morita; https://github.com/MotokiShiga/sova-cui).  Needs only numpy and sovapy, so that it can run with the Python of
a separate SOVA installation (soapff_cui.analysis calls it that way when sovapy is not importable).

For frames of traj.lammpstrj after the equilibration, for each subset of atoms a sovapy Atoms object (TriclinicVolume of
the LAMMPS cell, periodic) is built and, per frame,
    sovapy.computation.structure_factor.atom_pair_hist   -> pair-distance histogram
    sovapy.computation.structure_factor.pair_dist_func   -> partial pair distribution functions (PDF) g(r)
are evaluated and averaged over the frames; the radial distribution functions (RDF) are
    sovapy.computation.structure_factor.radial_dist_fun  -> N(r) = 4 pi r^2 rho_B g(r)   (rho_B = N_B / V)
SOVA pairs atoms by element, so water O and mineral O are separated by taking subsets:
    water O                 -> water O - water O
    each ion + water O      -> ion - water O
    each ion + mineral O    -> ion - mineral O
    all atoms               -> every element pair (all:<A>-<B>)

    python sova_structure.py --traj traj.lammpstrj --charges charges.csv --out DIR [--skip-ps 10] [--dt-fs 1]
        -> DIR/sova_pdf.csv, DIR/sova_rdf.csv, DIR/sova_summary.json
"""
import argparse
import csv
import json
import os

import numpy as np

ATOMIC_NUMBER = {"H": 1, "O": 8, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "Cl": 17, "K": 19, "Ca": 20, "Ti": 22, "Fe": 26, "Cs": 55}


def frames(path):
    """(step, cell rows, positions wrapped into the cell) per frame of a LAMMPS custom dump, sorted by id."""
    with open(path) as f:
        while True:
            l = f.readline()
            if not l:
                return
            if not l.startswith("ITEM: TIMESTEP"):
                continue
            try:
                step = int(f.readline()); f.readline(); n = int(f.readline())
                box = [f.readline() for _ in range(4)]
                cols = f.readline().split()[2:]
                lines = [f.readline() for _ in range(n)]
                if not lines or not lines[-1].endswith("\n"):
                    return
            except ValueError:
                return
            b = [list(map(float, x.split())) for x in box[1:4]]
            if "xy" in box[0]:
                xy, xz, yz = b[0][2], b[1][2], b[2][2]
                xlo = b[0][0] - min(0.0, xy, xz, xy + xz); xhi = b[0][1] - max(0.0, xy, xz, xy + xz)
                ylo = b[1][0] - min(0.0, yz); yhi = b[1][1] - max(0.0, yz)
            else:
                xy = xz = yz = 0.0
                xlo, xhi, ylo, yhi = b[0][0], b[0][1], b[1][0], b[1][1]
            cell = np.array([[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, b[2][1] - b[2][0]]])
            keys = ("x", "y", "z") if "x" in cols else ("xu", "yu", "zu")
            a = np.loadtxt(lines, usecols=[cols.index("id")] + [cols.index(c) for c in keys], ndmin=2)
            a = a[np.argsort(a[:, 0])]
            r = a[:, 1:4] - np.array([xlo, ylo, b[2][0]])
            yield step, cell, np.mod(r @ np.linalg.inv(cell), 1.0) @ cell


def sova_atoms(pos, symbols, cell):
    from sovapy.core import data
    from sovapy.core.volumes import TriclinicVolume
    order = np.argsort([ATOMIC_NUMBER.get(s, 99) for s in symbols], kind="stable")
    vol = TriclinicVolume(*cell[0], *cell[1], *cell[2])
    vol.origin = np.zeros(3)
    vol.periodic_boundary = [0.0, 1.0]
    return data.Atoms(pos[order], None, [symbols[i] for i in order], vol)


def first_peak(r, g):
    """position and height of the first peak of g(r) between 1.5 and 4 A (curve smoothed over 0.1 A for the search)."""
    dr = r[1] - r[0]
    w = max(1, int(round(0.1 / dr)))
    gs = np.convolve(g, np.ones(2 * w + 1) / (2 * w + 1), mode="same")
    lo, hi = int(1.5 / dr), min(len(gs), int(4.0 / dr))
    if hi <= lo or gs[lo:hi].max() <= 0:
        return None
    i = lo + int(np.argmax(gs[lo:hi]))
    return dict(r_peak=float(r[i]), g_peak=float(g[i]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True); ap.add_argument("--charges", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--skip-ps", type=float, default=0.0); ap.add_argument("--dt-fs", type=float, default=1.0)
    ap.add_argument("--max-frames", type=int, default=40); ap.add_argument("--dr", type=float, default=0.02)
    ap.add_argument("--rmax", type=float, default=8.0)
    a = ap.parse_args()
    import sovapy
    from sovapy.computation.structure_factor import atom_pair_hist, pair_dist_func, radial_dist_fun

    rows = list(csv.DictReader(open(a.charges)))
    el = np.array([r["element"] for r in rows]); role = np.array([r["role"] for r in rows])
    typ = np.array([r["clayff_type"] for r in rows])
    Ow = np.nonzero(typ == "ow")[0]
    ions = {e: np.nonzero((role == "ion") & (el == e))[0] for e in sorted(set(el[role == "ion"]))}
    Omin = np.nonzero((role == "mineral") & (el == "O"))[0]
    subsets = []                       # (label, atom indices, pairs to report or None = every element pair)
    if len(Ow) > 1:
        subsets.append(("Ow-Ow", Ow, [("O", "O")]))
    for e, idx in ions.items():
        if len(Ow):
            subsets.append((f"{e}-Ow", np.r_[idx, Ow], [(e, "O")]))
        if len(Omin):
            subsets.append((f"{e}-O_min", np.r_[idx, Omin], [(e, "O")]))
    subsets.append(("all", np.arange(len(el)), None))

    allf = [(s, c, p) for s, c, p in frames(a.traj) if s * a.dt_fs / 1000.0 >= a.skip_ps - 1e-9]
    if not allf:
        raise SystemExit("no frames after the equilibration")
    sel = [allf[i] for i in np.unique(np.linspace(0, len(allf) - 1, min(a.max_frames, len(allf))).astype(int))]
    widths = []
    for _, c, _ in sel:
        V = abs(np.linalg.det(c))
        widths.append(min(V / np.linalg.norm(np.cross(c[(k + 1) % 3], c[(k + 2) % 3])) for k in range(3)))
    rmax = min(a.rmax, 0.5 * min(widths))
    nr = int(rmax / a.dr) + 1
    vol = float(np.mean([abs(np.linalg.det(c)) for _, c, _ in sel]))

    pdf, rdf, peaks = {}, {}, {}
    r = None
    for label, idx, want in subsets:
        g_sum, atoms = None, None
        for _, cell, pos in sel:
            atoms = sova_atoms(pos[idx], [el[i] for i in idx], cell)
            rr, hist = atom_pair_hist(atoms, a.dr)
            g = pair_dist_func(atoms, rr, hist)[:nr]
            g_sum = g if g_sum is None else g_sum + g
        r, g = rr[:nr], g_sum / len(sel)
        pairs = [tuple(p) for p in atoms.pairs]
        for pa in (want or pairs):
            k = pairs.index(pa) if pa in pairs else pairs.index(pa[::-1])
            name = label if want else f"all:{pa[0]}-{pa[1]}"
            nb = int(np.sum(el[idx] == pa[1]))
            pdf[name] = g[:, k]
            rdf[name] = radial_dist_fun(r, g[:, k], nb / vol)            # N(r) = 4 pi r^2 rho_B g(r)
            if want:
                fp = first_peak(r, g[:, k])
                if fp:
                    peaks[name] = fp
    names = list(pdf)
    os.makedirs(a.out, exist_ok=True)
    for fname, table in (("sova_pdf.csv", pdf), ("sova_rdf.csv", rdf)):
        with open(os.path.join(a.out, fname), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["r_A"] + names)
            for i in range(len(r)):
                w.writerow([f"{r[i]:.3f}"] + [f"{table[n][i]:.6g}" for n in names])
    result = dict(sovapy_version=getattr(sovapy, "__version__", ""), frames=len(sel), rmax=float(rmax), dr=a.dr,
                  pairs=names, first_peaks=peaks,
                  note="PDF g(r) by SOVA pair_dist_func (normalised by the density over the whole box: in a slab or "
                       "interlayer it does not tend to 1); RDF N(r) = 4 pi r^2 rho_B g(r) by SOVA radial_dist_fun")
    json.dump(result, open(os.path.join(a.out, "sova_summary.json"), "w"), indent=1)
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
