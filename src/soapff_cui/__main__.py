import argparse

from .core import THRESHOLD, build_library, check_library, make, open_in_ovito


def main():
    ap = argparse.ArgumentParser(prog="soapff-cui", description="structure file -> LAMMPS ClayFF force field with SOAP charges")
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build-library", help="reference data file (+ dump frames) -> SOAP/charge library")
    b.add_argument("data", nargs="+", help="LAMMPS data file(s) with charges (atom_style full, type names in the Masses comments)")
    b.add_argument("--frames", nargs="*", default=[], help="MD dump frames of a reference: N:dump (N-th data file, "
                                                            "1-based) or dump (= first); thermal environments")
    b.add_argument("-o", "--out", default="library.npz")
    b.add_argument("--every", type=int, default=1, help="use every n-th dump frame")
    b.add_argument("--prune", type=float, default=None, help="drop entries at least this similar (e.g. 0.999) to a kept "
                                                             "entry of the same element, charge and type (smaller library)")
    b.add_argument("--water-in-soap", action="store_true", help="keep the water in the SOAP of the references (old libraries); "
                                                                "default: water left out, typed by rule in make")
    b.add_argument("--protect", type=int, default=0, help="with --prune: the first N reference files are kept whole")
    b.add_argument("--r-cut", type=float, default=5.0)
    b.add_argument("--l-max", type=int, default=3)
    b.add_argument("--weight-r0", type=float, default=None, help="radial weighting of SOAP, w(r) = 1/(1 + (r/r0)^m) "
                                                                  "(dscribe 'pow'): r0 in A; none by default")
    b.add_argument("--weight-m", type=float, default=4, help="exponent m of the radial weighting")
    b.add_argument("--species-weight", nargs="*", default=[], metavar="EL=W",
                   help="scale the SOAP components of every species pair (Z1, Z2) by w(Z1) w(Z2), e.g. Al=2 Si=2 Mg=2")
    sf = sub.add_parser("make-surfaces", help="bulk reference data file -> slabs in vacuum of every face (h k l) of its unit "
                                              "cell (surface references for build-library)")
    sf.add_argument("data", help="LAMMPS data file with charges and type names (water is left out)")
    sf.add_argument("-o", "--out", default="surfaces")
    sf.add_argument("--max-index", type=int, default=1, help="faces with |h|, |k|, |l| <= this")
    sf.add_argument("--thickness", type=float, default=16.0, help="minimum slab thickness (A)")
    sf.add_argument("--vacuum", type=float, default=16.0, help="vacuum along the normal (A)")
    c = sub.add_parser("check-library", help="developer check: nearest other entry has the same charge (one to one)")
    c.add_argument("library")
    m = sub.add_parser("make", help="structure -> system.data, system.settings, in.lammps, charges.csv, report.json")
    m.add_argument("structure")
    m.add_argument("-o", "--out", default="soapff_out")
    m.add_argument("--library", required=True)
    m.add_argument("--threshold", type=float, default=THRESHOLD, help="below this similarity: ClayFF type charge")
    m.add_argument("--cutoff", type=float, default=10.0)
    m.add_argument("--kspace", default="pppm 1.0e-5")
    m.add_argument("--temp", type=float, default=300.0)
    m.add_argument("--nvt-ps", type=float, default=50.0)
    m.add_argument("--timestep", type=float, default=0.5, help="time step (fs); 0.5 for the flexible SPC water of ClayFF")
    m.add_argument("--minimize-only", action="store_true", help="write in.lammps only up to the energy minimisation "
                                                                 "(the MD is added by the user)")
    m.add_argument("--relax-ps", type=float, default=10.0, help="NVT relaxation after the minimisation, before the "
                                                                "production (not recorded; 0: none)")
    m.add_argument("--water-rule", action="store_true", help="water molecules (O with two H, no framework cation) are "
                                                             "typed by rule (ow/hw, SPC), not by SOAP")
    m.add_argument("--no-check", action="store_true", help="write the input even if the molecule/charge consistency "
                                                           "check fails (problems are reported)")
    m.add_argument("--no-fallback", action="store_true", help="atoms below the threshold keep the nearest library entry "
                                                              "(reported) instead of the geometric ClayFF type/charge")
    m.add_argument("--dump-ps", type=float, default=1.0, help="trajectory output interval (ps); 0.1-0.5 for diffusion")
    m.add_argument("--npt-ps", type=float, default=0.0, help="NPT after the NVT (0: no NPT)")
    m.add_argument("--press", type=float, default=1.0, help="NPT pressure (atm)")
    m.add_argument("--baro", default="aniso", choices=["iso", "aniso", "tri", "z"],
                   help="NPT coupling: aniso/tri for bulk crystals, z for slabs and interfaces")
    m.add_argument("--charges", default="clayff2004", choices=["clayff2004", "underwood-bourg", "library"],
                   help="charge scheme: clayff2004 (Cygan et al. 2004), underwood-bourg (Mg 1.361 and neighbour O, smectite "
                        "edges), library (the charges stored in the library)")
    m.add_argument("--water", default="spc", choices=["spc", "spce"], help="water model: flexible SPC (ClayFF) or flexible SPC/E")
    m.add_argument("--boundary", default="p p p",
                   help="LAMMPS boundary, e.g. 'p p p' (default) or 'p p f' (slab with vacuum along z; NVT only)")
    o = sub.add_parser("open", help="open an output directory or file (traj.lammpstrj, data file) in OVITO")
    o.add_argument("path")
    an = sub.add_parser("analyze", help="analyses of an MD run: LAMMPS compute msd, fix ave/* files, "
                                        "SOVA (g(r), N(r), coordination, S(Q)) and Debyer (XRD)")
    an.add_argument("outdir")
    an.add_argument("--skip-ps", type=float, default=None, help="equilibration to skip (ps); default 20%% of the MD")
    an.add_argument("--what", default="msd,sova,xrd",
                    help="analyses, comma separated: msd (diffusion), sova (PDF, RDF, S(Q)), xrd (Debyer)")
    an.add_argument("--sova-dr", type=float, default=0.02, help="SOVA: bin width of g(r) (A)")
    an.add_argument("--sova-rmax", type=float, default=8.0, help="SOVA: largest r (A)")
    an.add_argument("--sova-qmin", type=float, default=0.05); an.add_argument("--sova-qmax", type=float, default=12.0)
    an.add_argument("--sova-dq", type=float, default=0.05, help="SOVA: Q range and step of S(Q) (1/A)")
    an.add_argument("--sova-radiation", default="xray", choices=["xray", "neutron"])
    an.add_argument("--sova-frames", type=int, default=40, help="SOVA: number of frames averaged")
    an.add_argument("--traj", default=None, help="dump for SOVA and XRD (default traj.lammpstrj, else the first *.lammpstrj)")
    an.add_argument("--debyer", default=None, help="debyer executable (default $DEBYER, PATH, ~/*/debyer/debyer/debyer)")
    an.add_argument("--wavelength", type=float, default=1.5406, help="X-ray wavelength in A (default Cu K-alpha 1.5406)")
    an.add_argument("--two-theta", default="3,80,0.02", help="2 theta from,to,step in degrees (default 3,80,0.02)")
    an.add_argument("--xrd-frames", type=int, default=3, help="number of last frames averaged for XRD (default 3)")
    an.add_argument("--xrd-mode", default="cluster", choices=["cluster", "periodic"],
                    help="cluster: the cell repeated to a crystallite of --xrd-size A (crystals, minerals; default); "
                         "periodic: periodic box with cut-off half the box (liquids, glasses; orthorhombic box only)")
    an.add_argument("--xrd-size", type=float, default=50.0, help="minimum crystallite edge in A for the cluster mode (default 50)")
    a = ap.parse_args()
    if a.cmd == "analyze":
        from .analysis import analyze
        analyze(a.outdir, skip_ps=a.skip_ps, sova_opts=dict(dr=a.sova_dr, rmax=a.sova_rmax, qmin=a.sova_qmin, qmax=a.sova_qmax, dq=a.sova_dq, radiation=a.sova_radiation, max_frames=a.sova_frames),
                what=tuple(w.strip() for w in a.what.split(",") if w.strip()), debyer=a.debyer, traj=a.traj,
                xrd_opts=dict(wavelength=a.wavelength, two_theta=tuple(float(v) for v in a.two_theta.split(",")),
                              frames=a.xrd_frames, mode=a.xrd_mode, size_A=a.xrd_size))
    elif a.cmd == "open":
        print("opening", open_in_ovito(a.path))
    elif a.cmd == "check-library":
        import json
        r = check_library(a.library)
        json.dump(r, open(a.library.replace(".npz", "_check.json"), "w"), indent=1)
    elif a.cmd == "make-surfaces":
        from .surfaces import write_surfaces
        write_surfaces(a.data, a.out, max_index=a.max_index, min_thickness=a.thickness, vacuum=a.vacuum)
    elif a.cmd == "build-library":
        build_library(a.data, a.frames, a.out, every=a.every, soap=dict(r_cut=a.r_cut, l_max=a.l_max, species_weights={x.split('=')[0]: float(x.split('=')[1]) for x in a.species_weight} or None, weighting=dict(function='pow', r0=a.weight_r0, m=a.weight_m, c=1, d=1) if a.weight_r0 else None), prune=a.prune,
                      water_in_soap=a.water_in_soap, protect=a.protect)
    else:
        make(a.structure, a.out, a.library, threshold=a.threshold, cutoff=a.cutoff, kspace=a.kspace, temp=a.temp,
             nvt_ps=a.nvt_ps, npt_ps=a.npt_ps, press=a.press, baro=a.baro, boundary=a.boundary, dump_ps=a.dump_ps, fallback=not a.no_fallback, check=not a.no_check, water_rule=a.water_rule, relax_ps=a.relax_ps, timestep_fs=a.timestep, md=not a.minimize_only,
             charges=a.charges, water=a.water)


if __name__ == "__main__":
    main()
