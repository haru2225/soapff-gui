"""Physical properties of a soapff-cui system: the LAMMPS commands that compute them (added to in.lammps after the energy
minimisation) and the analysis of the files LAMMPS writes.  Every method is the one of the LAMMPS manual and examples:

  structure  dump of the equilibrium MD -> SOVA (PDF, RDF, coordination numbers, S(Q)) and Debyer (XRD), soapff_cui.analysis
  diffusion  compute msd of water O and each ion type -> self-diffusion coefficients (Einstein) and the ionic conductivity
             (Nernst-Einstein: sigma = e^2 / (V kB T) sum_i N_i z_i^2 D_i, ion-ion correlations neglected)
  viscosity  Green-Kubo (Howto viscosity): fix ave/correlate of pxy, pxz, pyz; eta = V / (kB T) int <P_ab(0) P_ab(t)> dt
             (a liquid property: for a solid or a clay with interlayer water the integral does not converge)
  thermal    Green-Kubo (compute heat/flux, Howto kappa): fix ave/correlate of the heat flux J;
             kappa = 1 / (V kB T^2) int <J(0) J(t)> dt.  Per-atom stress: compute centroid/stress/atom for pair, bond and
             angle (recommended for molecules with angles, but it has no kspace part) plus compute stress/atom for the
             PPPM part; the heat flux is linear in the stress, so J = J(centroid) + J(kspace) - J(no stress)
  elastic0   0 K elastic constants (examples/ELASTIC; strain 0.005 by default: with PPPM the minimisations stop at a
             residual stress of ~50-100 atm, which with a strain of 0.001 makes C_ij uncertain by 5-10 GPa; 0.005 and
             0.01 give the same C_ij within 1-2 GPa): box relaxed to zero stress (fix box/relax, minimize), then for
             each strain direction (xx yy zz yz xz xy) +-strain applied by change_box (remap), atoms relaxed (minimize),
             C_ij = -(p_i(+) - p_i(-)) / (2 strain)
  elasticT   elastic constants at temperature T (examples/ELASTIC_T/BORN_MATRIX): compute born/matrix numdiff (Born
             term), the covariance of the virial stress and the kinetic term:
             C_ij = <C^B_ij> - V/(kB T) Cov(s_i, s_j) + N kB T / V (2 for C11 C22 C33, 1 for C44 C55 C66)
  tensile    uniaxial tension / compression: fix deform erate along x, y or z, the other two at constant pressure by
             fix npt; stress = -pressure; Young's modulus = initial slope of stress-strain, strength = peak stress
  shear      fix deform of a tilt factor (xy, xz, yz): solid (fix nvt, remap x) -> shear modulus and strength;
             flow (fix nvt/sllod with compute temp/deform, remap v) -> viscosity = shear stress / shear rate

Order in in.lammps: elastic0 (right after the minimisation, 0 K) -> equilibration at T (NVT, then NPT if asked; the
box is then set to its NPT average) -> one equilibrium MD collecting structure, diffusion, viscosity, thermal and
elasticT together (NVE when a Green-Kubo property is asked and elasticT is not, else NVT) -> each tensile / shear test,
each starting from the equilibrated state (read_dump of soapff_eq.dump).  Each block starts with a marker
'# soapff-property: <name> key=value ...' that analyze() reads.
"""
import csv
import glob
import os

import numpy as np

ALL_PROPS = (("structure", "構造：SOVA の PDF・RDF・配位数・S(Q)、debyer の XRD"),
             ("diffusion", "拡散係数とイオン伝導度（MSD、ネルンスト–アインシュタイン）"),
             ("viscosity", "粘度（グリーン–久保：応力の自己相関。液体のみ）"),
             ("thermal", "熱伝導率（グリーン–久保：熱流束の自己相関）"),
             ("elastic0", "弾性定数 0 K（ELASTIC：± 微小ひずみと最小化）"),
             ("elasticT", "弾性定数 有限温度（ELASTIC_T：born/matrix とゆらぎの式）"),
             ("tensile", "引張・圧縮（fix deform ＋ 横方向は fix npt）"),
             ("shear", "せん断（fix deform の傾き。流れは fix nvt/sllod）"))
# the properties offered in the GUI: those that do not depend on the mineral (any periodic structure); the others
# (structure, diffusion, viscosity, thermal) are still generated and analysed when in.lammps has them
PROPS = tuple(x for x in ALL_PROPS if x[0] in ("elastic0", "elasticT", "tensile", "shear"))
EQUILIBRIUM = ("structure", "diffusion", "viscosity", "thermal", "elasticT")
E0_MIN = "0.0 1.0e-8 20000 200000"          # minimisation of the reference and strained states (etol ftol maxiter maxeval)
MECHANICAL = ("tensile", "shear")
DEFAULTS = dict(temp=300.0, press=1.0, timestep=0.5, relax_ps=10.0, npt_ps=20.0, baro="aniso", prod_ps=100.0, dump_ps=1.0,
                msd_ps=0.1, gk_sample=4, gk_corr_ps=5.0, e0_strain=5.0e-3, eT_stress_every=10, eT_born_every=100,
                eT_delta=1.0e-4, tensile_dirs="x", tensile_mode="tension", tensile_rate=1.0e9, tensile_max=0.1,
                shear_planes="xy", shear_mode="solid", shear_rate=1.0e9, shear_max=0.2, mech_dump_ps=1.0, seed=4928459)

# constants (SI) and LAMMPS real units
KB, E, ATM2PA, A2M, FS2S = 1.380649e-23, 1.602176634e-19, 101325.0, 1.0e-10, 1.0e-15
KCALMOL2J = 4184.0 / 6.02214076e23
BOLTZ, NKTV2P, ATM2GPA = 0.0019872067, 68568.415, 1.01325e-4      # kcal/mol/K; atm per kcal/mol/A^3; GPa per atm
VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")
PCOMP = {"xx": "pxx", "yy": "pyy", "zz": "pzz", "yz": "pyz", "xz": "pxz", "xy": "pxy"}


def needs_triclinic(props):
    """elastic0 (shear strains) and shear change tilt factors: the box must be triclinic, and 'change_box all triclinic'
    must come before any dump (LAMMPS rule)."""
    return "elastic0" in props or "shear" in props


def triclinic_lines(out_dir):
    """'change_box all triclinic' and the kspace_style of system.settings again (LAMMPS: the kspace style must be
    redefined after changing to a triclinic box); placed right after the minimisation, before any dump."""
    out = ["change_box all triclinic"]
    try:
        out += [l.strip() for l in open(os.path.join(out_dir, "system.settings")) if l.split()[:1] == ["kspace_style"]]
    except OSError:
        pass
    return out


def lines(props, p, type_ids, msd_lines=None, existing=()):
    """LAMMPS commands for the properties `props` with parameters p (DEFAULTS); type_ids {ClayFF type name: type number}
    (the element names of the dump and the MSD groups); msd_lines(type_ids, every): the MSD block of soapff_cui.core;
    existing: fix and compute IDs already in in.lammps (an MSD block that is already there is not repeated)."""
    p = dict(DEFAULTS, **{k: v for k, v in p.items() if v is not None})
    props = [k for k, _ in ALL_PROPS if k in props]
    T, P, dt = float(p["temp"]), float(p["press"]), float(p["timestep"])
    steps = lambda ps: max(1, int(round(float(ps) * 1000.0 / dt)))
    tdamp, pdamp = 100.0, 1000.0                                   # fs (100 and 1000 time steps of 1 fs; ClayFF studies)
    names = " ".join(t for t, _ in sorted(type_ids.items(), key=lambda kv: kv[1]))
    out = ["", "# ===== soapff：物性の計算（" + "、".join(lab.split("（")[0] for k, lab in ALL_PROPS if k in props) + "）====="]
    if needs_triclinic(props):
        out += ["# (change_box all triclinic is needed for the tilt factors: it is placed before any dump)"]

    # --- 0 K elastic constants
    if "elastic0" in props:
        e = float(p["e0_strain"])
        out += ["", f"# soapff-property: elastic0 strain={e:g} file=elastic0.dat",
                "# 0 K elastic constants (LAMMPS examples/ELASTIC): zero-stress box, then +-strain per direction, atoms relaxed",
                "velocity all set 0.0 0.0 0.0",
                "# zero-stress box: box/relax and fixed-box minimisations alternated 4 times (box/relax alone leaves the",
                "# atoms unrelaxed); the last fixed-box minimisation is the reference, relaxed as the strained states"]
        for _ in range(4):
            out += ["fix soapff_br all box/relax tri 0.0 vmax 0.001 nreset 20", f"minimize {E0_MIN}", "unfix soapff_br",
                    f"minimize {E0_MIN}"]
        out += [
                'print "# dir sign strain pxx pyy pzz pyz pxz pxy (atm) lx ly lz xy xz yz" file elastic0.dat screen no',
                'print "0 0 0 $(pxx) $(pyy) $(pzz) $(pyz) $(pxz) $(pxy) $(lx) $(ly) $(lz) $(xy) $(xz) $(yz)" append elastic0.dat screen no',
                "reset_timestep 0",
                "write_dump all custom soapff_e0.dump id x y z ix iy iz modify sort id"]
        for d in range(1, 7):
            for sgn in (-1, 1):
                u = sgn * e
                box = {1: f"x delta 0 $({u:g}*lx) xy delta $({u:g}*xy) xz delta $({u:g}*xz)",
                       2: f"y delta 0 $({u:g}*ly) yz delta $({u:g}*yz)",
                       3: f"z delta 0 $({u:g}*lz)",
                       4: f"yz delta $({u:g}*lz)", 5: f"xz delta $({u:g}*lz)", 6: f"xy delta $({u:g}*ly)"}[d]
                out += ["read_dump soapff_e0.dump 0 x y z ix iy iz box yes replace yes",
                        f"change_box all {box} remap units box",
                        f"minimize {E0_MIN}",
                        f'print "{d} {sgn} {u:g} $(pxx) $(pyy) $(pzz) $(pyz) $(pxz) $(pxy) $(lx) $(ly) $(lz) $(xy) $(xz) $(yz)" '
                        "append elastic0.dat screen no"]
        out += ["read_dump soapff_e0.dump 0 x y z ix iy iz box yes replace yes"]

    finite = [k for k in props if k != "elastic0"]
    if not finite:
        return out

    # --- equilibration at T
    out += ["", f"# soapff: equilibration at {T:g} K (NVT {p['relax_ps']:g} ps" +
            (f", NPT {p['npt_ps']:g} ps at {P:g} atm, {p['baro']}; box set to its NPT average)" if float(p["npt_ps"]) > 0 else ")"),
            f"velocity all create {T:g} {int(p['seed'])} mom yes rot no dist gaussian",
            f"timestep {dt:g}",
            f"fix soapff_eq all nvt temp {T:g} {T:g} {tdamp:g}",
            f"run {steps(p['relax_ps'])}",
            "unfix soapff_eq"]
    if float(p["npt_ps"]) > 0:
        n = steps(p["npt_ps"])
        n = max(2000, int(round(n / 2000.0)) * 2000)               # two halves, each a multiple of the averaging output
        out += ["reset_timestep 0",                                # the box average is output on multiples of 1000
                f"fix soapff_eq all npt temp {T:g} {T:g} {tdamp:g} {p['baro']} {P:g} {P:g} {pdamp:g}",
                f"run {n // 2}",
                "# the box averaged over the second half of the NPT only (the first half still relaxes the volume)",
                "variable soapff_lx equal lx", "variable soapff_ly equal ly", "variable soapff_lz equal lz",
                "variable soapff_xy equal xy", "variable soapff_xz equal xz", "variable soapff_yz equal yz",
                "fix soapff_box all ave/time 10 100 1000 v_soapff_lx v_soapff_ly v_soapff_lz v_soapff_xy v_soapff_xz "
                "v_soapff_yz ave running",
                f"run {n // 2}",
                "unfix soapff_eq",
                "change_box all x final $(xlo) $(xlo+f_soapff_box[1]) y final $(ylo) $(ylo+f_soapff_box[2]) "
                "z final $(zlo) $(zlo+f_soapff_box[3]) remap units box"]
        out += ["change_box all xy final $(f_soapff_box[4]) xz final $(f_soapff_box[5]) yz final $(f_soapff_box[6]) "
                "remap units box"] if needs_triclinic(props) or p["baro"] == "tri" else []
        out += ["unfix soapff_box"]
    out += ["reset_timestep 0"]
    if any(k in props for k in MECHANICAL):
        out += ["write_dump all custom soapff_eq.dump id x y z vx vy vz ix iy iz modify sort id"]

    # --- equilibrium MD
    eq = [k for k in props if k in EQUILIBRIUM]
    if eq:
        gk = [k for k in eq if k in ("viscosity", "thermal")]
        nve = bool(gk) and "elasticT" not in eq
        N = steps(p["prod_ps"])
        s = max(1, int(p["gk_sample"]))
        nc = max(2, int(round(float(p["gk_corr_ps"]) * 1000.0 / dt / s)))
        out += ["", f"# soapff: equilibrium MD {p['prod_ps']:g} ps ({'NVE' if nve else 'NVT'}) for: " + ", ".join(eq),
                "fix soapff_md all nve" if nve else f"fix soapff_md all nvt temp {T:g} {T:g} {tdamp:g}",
                "variable soapff_vol equal vol",
                "fix soapff_state all ave/time 10 10 100 c_thermo_temp c_thermo_press v_soapff_vol c_thermo_pe "
                "file soapff_state.dat"]
        undo = ["unfix soapff_md", "unfix soapff_state"]
        if "structure" in eq:
            out += [f"# soapff-property: structure file=traj.lammpstrj dump_ps={p['dump_ps']:g}",
                    f"dump soapff_d all custom {steps(p['dump_ps'])} traj.lammpstrj id type element x y z ix iy iz q",
                    f"dump_modify soapff_d sort id element {names}"]
            undo += ["undump soapff_d"]
        if "diffusion" in eq:
            ml = msd_lines(type_ids, steps(p["msd_ps"])) if msd_lines else []
            if any(l.split()[1] in existing for l in ml if l.startswith(("fix ", "compute "))):
                ml = ["# (the MSD of soapff is already in in.lammps above: it is not repeated)"]
            out += ["# soapff-property: diffusion files=msd_*.dat"] + ml
            undo += [f"unfix {l.split()[1]}" for l in ml if l.startswith("fix ")]
            undo += [f"uncompute {l.split()[1]}" for l in ml if l.startswith("compute ")]
        if "viscosity" in eq:
            out += [f"# soapff-property: viscosity file=gk_viscosity.dat sample={s} corr={nc}",
                    "# Green-Kubo viscosity (LAMMPS Howto viscosity): autocorrelation of the off-diagonal pressure",
                    "# (c_thermo_press[4 5 6] = pxy pxz pyz; the default pressure compute, whatever the thermo_style)",
                    f"fix soapff_gkv all ave/correlate {s} {nc} {s * nc} c_thermo_press[4] c_thermo_press[5] c_thermo_press[6] "
                    "type auto file gk_viscosity.dat ave running"]
            undo += ["unfix soapff_gkv"]
        if "thermal" in eq:
            out += [f"# soapff-property: thermal file=gk_heatflux.dat sample={s} corr={nc}",
                    "# Green-Kubo thermal conductivity (compute heat/flux): centroid stress for pair, bond, angle; per-atom",
                    "# stress of PPPM by stress/atom (centroid/stress/atom has no kspace part); J = J_c + J_k - J_0",
                    "compute soapff_ke all ke/atom", "compute soapff_pe all pe/atom",
                    "compute soapff_sc all centroid/stress/atom NULL pair bond angle",
                    "compute soapff_sk all stress/atom NULL kspace",
                    "compute soapff_s0 all stress/atom NULL fix",
                    "compute soapff_jc all heat/flux soapff_ke soapff_pe soapff_sc",
                    "compute soapff_jk all heat/flux soapff_ke soapff_pe soapff_sk",
                    "compute soapff_j0 all heat/flux soapff_ke soapff_pe soapff_s0"]
            for i, a in enumerate("xyz", 1):
                out += [f"variable soapff_j{a} equal c_soapff_jc[{i}]+c_soapff_jk[{i}]-c_soapff_j0[{i}]"]
            out += [f"fix soapff_gkk all ave/correlate {s} {nc} {s * nc} v_soapff_jx v_soapff_jy v_soapff_jz "
                    "type auto file gk_heatflux.dat ave running"]
            undo += ["unfix soapff_gkk"] + [f"uncompute soapff_{c}" for c in ("jc", "jk", "j0", "sc", "sk", "s0", "ke", "pe")]
        if "elasticT" in eq:
            ns, nb = max(1, int(p["eT_stress_every"])), max(1, int(p["eT_born_every"]))
            out += [f"# soapff-property: elasticT born=elasticT_born.dat virial=elasticT_virial.dat delta={float(p['eT_delta']):g}",
                    "# elastic constants at T (LAMMPS examples/ELASTIC_T/BORN_MATRIX): Born term, virial stress, kinetic term",
                    "compute soapff_vir all pressure NULL virial",
                    f"compute soapff_born all born/matrix numdiff {float(p['eT_delta']):g} soapff_vir",
                    f"fix soapff_bornf all ave/time {nb} 1 {nb} c_soapff_born[*] file elasticT_born.dat",
                    f"fix soapff_virf all ave/time {ns} 1 {ns} c_soapff_vir[*] c_thermo_temp v_soapff_vol file elasticT_virial.dat"]
            undo += ["unfix soapff_bornf", "unfix soapff_virf", "uncompute soapff_born", "uncompute soapff_vir"]
        out += [f"run {N}"] + undo

    # --- mechanical tests
    def restore():
        return ["read_dump soapff_eq.dump 0 x y z vx vy vz ix iy iz box yes replace yes", "reset_timestep 0"]
    stress_vars = ["variable soapff_lx equal lx", "variable soapff_ly equal ly", "variable soapff_lz equal lz",
                   "variable soapff_xy equal xy", "variable soapff_xz equal xz", "variable soapff_yz equal yz"]
    # pressure tensor of the default compute thermo_press (pxx pyy pzz pxy pxz pyz), written as pxx pyy pzz pyz pxz pxy
    rec = ("c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] c_thermo_press[6] c_thermo_press[5] c_thermo_press[4] "
           "v_soapff_lx v_soapff_ly v_soapff_lz v_soapff_xy v_soapff_xz v_soapff_yz c_thermo_temp")
    if "tensile" in props:
        dirs = [d for d in str(p["tensile_dirs"]).replace(",", " ").split() if d in ("x", "y", "z")] or ["x"]
        modes = {"tension": [1], "compression": [-1], "both": [1, -1]}.get(p["tensile_mode"], [1])
        rate = float(p["tensile_rate"]) * FS2S                     # 1/s -> 1/fs (real units)
        n = steps(float(p["tensile_max"]) / (rate * 1.0e3))         # time = strain / rate (fs -> ps)
        for d in dirs:
            for sg in modes:
                name = f"{'tensile' if sg > 0 else 'compress'}_{d}"
                lat = [a for a in "xyz" if a != d]
                out += ["", f"# soapff-property: tensile dir={d} sign={sg} rate_s={float(p['tensile_rate']):g} "
                        f"max_strain={float(p['tensile_max']):g} file=mech_{name}.dat",
                        f"# uniaxial {'tension' if sg > 0 else 'compression'} along {d}: fix deform erate, "
                        f"{lat[0]} and {lat[1]} at {P:g} atm (fix npt); stress = -p{d}{d}"] + restore() + stress_vars + [
                        f"variable soapff_L0 equal $(l{d})",
                        f"variable soapff_strain equal (l{d}-v_soapff_L0)/v_soapff_L0",
                        f"fix soapff_t all npt temp {T:g} {T:g} {tdamp:g} {lat[0]} {P:g} {P:g} {pdamp:g} "
                        f"{lat[1]} {P:g} {P:g} {pdamp:g} couple none",
                        f"fix soapff_def all deform 1 {d} erate {sg * rate:.6g} units box remap x",
                        f"fix soapff_mech all ave/time 10 10 100 v_soapff_strain {rec} file mech_{name}.dat",
                        f"dump soapff_dm all custom {steps(p['mech_dump_ps'])} traj_{name}.lammpstrj id type element x y z ix iy iz q",
                        f"dump_modify soapff_dm sort id element {names}",
                        f"run {n}",
                        "unfix soapff_mech", "unfix soapff_def", "unfix soapff_t", "undump soapff_dm"]
    if "shear" in props:
        planes = [a for a in str(p["shear_planes"]).replace(",", " ").split() if a in ("xy", "xz", "yz")] or ["xy"]
        flow = p["shear_mode"] == "flow"
        rate = float(p["shear_rate"]) * FS2S
        n = steps(float(p["shear_max"]) / (rate * 1.0e3))
        perp = {"xy": "ly", "xz": "lz", "yz": "lz"}
        for a in planes:
            name = f"shear_{a}"
            out += ["", f"# soapff-property: shear plane={a} mode={'flow' if flow else 'solid'} rate_s={float(p['shear_rate']):g} "
                    f"max_strain={float(p['shear_max']):g} file=mech_{name}.dat",
                    f"# shear {a}: fix deform of the tilt factor {a} (engineering strain rate, flip no); "
                    + ("flow: fix nvt/sllod (temp/deform), remap v" if flow else "solid: fix nvt, remap x")
                    + f"; shear stress = -p{a}"] + restore() + stress_vars + [
                    f"variable soapff_T0 equal $({a})",
                    f"variable soapff_strain equal ({a}-v_soapff_T0)/{perp[a]}"]
            if flow:                                                 # nvt/sllod makes its own compute temp/deform
                out += [f"fix soapff_t all nvt/sllod temp {T:g} {T:g} {tdamp:g}",
                        f"fix soapff_def all deform 1 {a} erate {rate:.6g} remap v flip no"]
            else:
                out += [f"fix soapff_t all nvt temp {T:g} {T:g} {tdamp:g}",
                        f"fix soapff_def all deform 1 {a} erate {rate:.6g} remap x flip no"]
            out += [f"fix soapff_mech all ave/time 10 10 100 v_soapff_strain {rec} file mech_{name}.dat",
                    f"dump soapff_dm all custom {steps(p['mech_dump_ps'])} traj_{name}.lammpstrj id type element x y z ix iy iz q",
                    f"dump_modify soapff_dm sort id element {names}",
                    f"run {n}",
                    "unfix soapff_mech", "unfix soapff_def", "unfix soapff_t", "undump soapff_dm"]
    return out


# ================================================================================================== analysis
def markers(out_dir):
    """[(name, {key: value})] of the '# soapff-property:' lines of in.lammps."""
    res = []
    try:
        for line in open(os.path.join(out_dir, "in.lammps")):
            if line.startswith("# soapff-property:"):
                w = line.split(":", 1)[1].split()
                res.append((w[0], dict(x.split("=", 1) for x in w[1:] if "=" in x)))
    except OSError:
        pass
    return res


def _rows(path):
    """numeric rows of a LAMMPS fix ave/time (scalar mode) or print file."""
    return np.array([list(map(float, l.split())) for l in open(path) if l.strip() and not l.startswith("#")])


def state(out_dir):
    """mean temperature (K) and volume (A^3) of the equilibrium MD (soapff_state.dat), or None."""
    p = os.path.join(out_dir, "soapff_state.dat")
    if not os.path.exists(p):
        return None
    a = _rows(p)
    return (float(a[:, 1].mean()), float(a[:, 3].mean())) if len(a) else None


def _last_block(path):
    """the last (running-average) block of a fix ave/correlate file: (TimeDelta in steps, values)."""
    blocks, cur, n = [], [], None
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        w = line.split()
        if n is None:
            n, cur = int(w[1]), []
            continue
        cur.append(list(map(float, w)))
        if len(cur) == n:
            blocks.append(np.array(cur))
            n = None
    if not blocks:
        raise ValueError("no complete output yet")
    b = blocks[-1]
    return b[:, 1], b[:, 3:], len(blocks)


def green_kubo(out_dir, m, dt_fs, kind):
    """viscosity (Pa s) or thermal conductivity (W/m/K) from the running-average autocorrelation."""
    st = state(out_dir)
    if not st:
        raise ValueError("soapff_state.dat (temperature, volume) is missing")
    T, V = st
    lag, c, nb = _last_block(os.path.join(out_dir, m["file"]))
    t = lag * dt_fs                                                            # fs
    run = np.array([np.concatenate([[0.0], np.cumsum(0.5 * (c[1:, k] + c[:-1, k]) * np.diff(t))]) for k in range(c.shape[1])]).T
    if kind == "viscosity":                                                    # atm^2 fs -> Pa^2 s, A^3 -> m^3
        run *= ATM2PA ** 2 * FS2S * (V * A2M ** 3) / (KB * T)
    else:                                                                      # (kcal/mol A/fs)^2 fs / A^3 -> W/m/K
        run *= KCALMOL2J ** 2 / FS2S / A2M / (V * KB * T ** 2)
    mean = run.mean(axis=1)
    half = t >= 0.5 * t[-1]
    return dict(T=T, V=V, t_fs=t, acf=c, running=run, mean=mean, outputs=nb,
                value=float(mean[half].mean()), value_end=float(mean[-1]),
                components=[float(run[half, k].mean()) for k in range(run.shape[1])])


def _moduli(C):
    """Voigt, Reuss, Hill bulk and shear moduli (GPa), and the Young's moduli and Poisson ratios of the compliance."""
    K_V = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * (C[0, 1] + C[0, 2] + C[1, 2])) / 9.0
    G_V = (C[0, 0] + C[1, 1] + C[2, 2] - C[0, 1] - C[0, 2] - C[1, 2] + 3 * (C[3, 3] + C[4, 4] + C[5, 5])) / 15.0
    r = dict(K_Voigt=float(K_V), G_Voigt=float(G_V))
    try:
        S = np.linalg.inv(C)
        K_R = 1.0 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * (S[0, 1] + S[0, 2] + S[1, 2]))
        G_R = 15.0 / (4 * (S[0, 0] + S[1, 1] + S[2, 2]) - 4 * (S[0, 1] + S[0, 2] + S[1, 2]) + 3 * (S[3, 3] + S[4, 4] + S[5, 5]))
        r.update(K_Reuss=float(K_R), G_Reuss=float(G_R), K_Hill=float((K_V + K_R) / 2), G_Hill=float((G_V + G_R) / 2),
                 E_x=float(1 / S[0, 0]), E_y=float(1 / S[1, 1]), E_z=float(1 / S[2, 2]),
                 nu_xy=float(-S[0, 1] / S[0, 0]), nu_xz=float(-S[0, 2] / S[0, 0]), nu_yz=float(-S[1, 2] / S[1, 1]))
    except np.linalg.LinAlgError:
        pass
    return r


def elastic0(out_dir, m):
    a = _rows(os.path.join(out_dir, m["file"]))
    e = float(m["strain"])
    C = np.full((6, 6), np.nan)
    for d in range(1, 7):
        neg, pos = a[(a[:, 0] == d) & (a[:, 1] == -1)], a[(a[:, 0] == d) & (a[:, 1] == 1)]
        if len(neg) and len(pos):
            C[:, d - 1] = -(pos[0, 3:9] - neg[0, 3:9]) / (2 * e) * ATM2GPA
    done = ~np.isnan(C).any(axis=0)
    Cs = 0.5 * (C + C.T)
    return dict(C=Cs, raw=C, directions=int(done.sum()), ref_stress_GPa=(-a[a[:, 0] == 0][0, 3:9] * ATM2GPA).tolist()
                if (a[:, 0] == 0).any() else None, moduli=_moduli(Cs) if done.all() else {})


def elasticT(out_dir, m, n_atoms):
    b = _rows(os.path.join(out_dir, m["born"]))[:, 1:22]                      # Born term, energy (kcal/mol)
    v = _rows(os.path.join(out_dir, m["virial"]))                             # step, pxx pyy pzz pxy pxz pyz, T, V
    T, V = float(v[:, 7].mean()), float(v[:, 8].mean())
    s = v[:, [1, 2, 3, 6, 5, 4]]                                               # Voigt order xx yy zz yz xz xy (atm)
    born = b.mean(axis=0)
    B = np.zeros((6, 6))
    k = 6
    for i in range(6):
        B[i, i] = born[i]
    for i in range(6):
        for j in range(i + 1, 6):
            B[i, j] = B[j, i] = born[k]
            k += 1
    B *= NKTV2P / V                                                           # -> atm
    F = -np.cov(s.T, bias=True) * V / (BOLTZ * T) / NKTV2P                    # atm
    K = np.diag([2, 2, 2, 1, 1, 1]) * n_atoms * BOLTZ * T / V * NKTV2P        # atm
    C = (B + F + K) * ATM2GPA
    return dict(C=C, born=B * ATM2GPA, fluct=F * ATM2GPA, kinetic=K * ATM2GPA, T=T, V=V, samples_stress=len(s),
                samples_born=len(b), moduli=_moduli(C))


def mechanical(out_dir, m, dt_fs, fit_max=None):
    """stress-strain curve (fix ave/time: step, strain, pxx pyy pzz pyz pxz pxy, lx ly lz xy xz yz, T)."""
    a = _rows(os.path.join(out_dir, m["file"]))
    strain = a[:, 1]
    if "dir" in m:
        d = m["dir"]
        sgn = int(m.get("sign", 1))
        stress = -a[:, 2 + "xyz".index(d)] * ATM2GPA                         # sigma = -p
        lat = [k for k in "xyz" if k != d]
        L = {k: a[:, 8 + "xyz".index(k)] for k in "xyz"}
        lat_strain = {k: (L[k] - L[k][0]) / L[k][0] for k in lat}
    else:
        comp = {"yz": 5, "xz": 6, "xy": 7}[m["plane"]]
        stress = -a[:, comp] * ATM2GPA
        sgn, lat_strain = 1, {}
    x, y = sgn * strain, sgn * stress                                         # compression: positive numbers
    fit_max = fit_max or min(0.02, 0.25 * float(m.get("max_strain", 0.1)))
    sel = (x >= 0) & (x <= fit_max)
    res = dict(strain=strain, stress=stress, fit_max=fit_max, points=int(len(a)))
    if sel.sum() >= 3:
        res["modulus"] = float(np.polyfit(x[sel], y[sel], 1)[0])
        for k, ls in lat_strain.items():
            res[f"poisson_{k}"] = float(-np.polyfit(x[sel], sgn * ls[sel], 1)[0])
    i = int(np.argmax(y))
    res.update(peak=float(y[i]), peak_strain=float(x[i]), initial_stress=float(np.mean(y[:3])))
    if m.get("mode") == "flow":                                               # NEMD viscosity: tau / rate (second half)
        half = x >= 0.5 * x.max()
        res["viscosity_Pa_s"] = float(np.mean(y[half]) * 1e9 / float(m["rate_s"]))
    return res


def ionic_conductivity(out_dir, D, temp=None):
    """Nernst-Einstein conductivity (S/m) of the ions with MSD results D {group: dict(D_total, D_xy, D_z)} (1e-9 m2/s);
    charges and numbers from charges.csv, T and V from soapff_state.dat (or the given temperature and system.data)."""
    rows = list(csv.DictReader(open(os.path.join(out_dir, "charges.csv"))))
    ions = {}
    for r in rows:
        if r.get("role") == "ion":
            ions.setdefault(r["element"], []).append(float(r.get("charge", r.get("q", 0.0)) or 0.0))
    st = state(out_dir)
    if st:
        T, V = st
    else:
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


def analyze(out_dir, dest, dt_fs, log=print):
    """every '# soapff-property:' block of in.lammps with its output: figures in dest, a summary per block."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res, figs = [], []
    n_atoms = sum(1 for _ in csv.DictReader(open(os.path.join(out_dir, "charges.csv"))))
    for name, m in markers(out_dir):
        if name in ("structure", "diffusion"):                               # analysed by soapff_cui.analysis
            continue
        item = dict(property=name, **m)
        try:
            if name in ("viscosity", "thermal"):
                g = green_kubo(out_dir, m, dt_fs, name)
                unit, scale = ("mPa·s", 1e3) if name == "viscosity" else ("W/m/K", 1.0)
                comps = ("xy", "xz", "yz") if name == "viscosity" else ("x", "y", "z")
                item.update(value=g["value"] * scale, value_end=g["value_end"] * scale, unit=unit, T=g["T"], V=g["V"],
                            components={c: v * scale for c, v in zip(comps, g["components"])}, outputs=g["outputs"],
                            corr_ps=float(g["t_fs"][-1] / 1000.0))
                fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
                for k, c in enumerate(comps):
                    ax[0].plot(g["t_fs"] / 1000, g["acf"][:, k] / g["acf"][0, k], lw=1, label=c)
                    ax[1].plot(g["t_fs"] / 1000, g["running"][:, k] * scale, lw=1, label=c)
                ax[1].plot(g["t_fs"] / 1000, g["mean"] * scale, lw=2, color="k", label="mean")
                ax[1].axhline(item["value"], color="#d62728", ls="--", lw=0.8)
                ax[1].axvspan(g["t_fs"][-1] / 2000, g["t_fs"][-1] / 1000, color="#999", alpha=0.1, lw=0)
                ax[0].axhline(0, color="#999", lw=0.6)
                ax[0].set_xlabel("t (ps)"); ax[0].set_ylabel("normalised autocorrelation"); ax[0].legend(fontsize=8)
                ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel(f"running integral ({unit})"); ax[1].legend(fontsize=8)
                ax[0].set_title(f"Green-Kubo {'viscosity: <P_ab(0)P_ab(t)>' if name == 'viscosity' else 'thermal conductivity: <J(0)J(t)>'}",
                                loc="left", fontsize=10)
                ax[1].set_title(f"{item['value']:.4g} {unit} (mean over the shaded half)", loc="left", fontsize=10)
                png = f"prop_{name}.png"
            elif name in ("elastic0", "elasticT"):
                r = elastic0(out_dir, m) if name == "elastic0" else elasticT(out_dir, m, n_atoms)
                item.update(C=np.round(r["C"], 3).tolist(), moduli=r["moduli"],
                            **({"directions": r["directions"], "ref_stress_GPa": r["ref_stress_GPa"]} if name == "elastic0"
                               else {"T": r["T"], "samples_stress": r["samples_stress"], "samples_born": r["samples_born"],
                                     "born": np.round(r["born"], 3).tolist(), "fluct": np.round(r["fluct"], 3).tolist(),
                                     "kinetic": np.round(r["kinetic"], 3).tolist()}))
                fig, ax = plt.subplots(figsize=(6.4, 5.6))
                C = r["C"]
                lim = np.nanmax(np.abs(C)) if np.isfinite(C).any() else 1
                ax.imshow(np.nan_to_num(C), cmap="RdBu_r", vmin=-lim, vmax=lim)
                for i in range(6):
                    for j in range(6):
                        ax.text(j, i, "" if np.isnan(C[i, j]) else f"{C[i, j]:.1f}", ha="center", va="center", fontsize=9,
                                color="white" if abs(np.nan_to_num(C[i, j])) > 0.6 * lim else "black")
                lab = [f"{k + 1}\n{v}" for k, v in enumerate(VOIGT)]
                ax.set_xticks(range(6)); ax.set_xticklabels(lab, fontsize=8); ax.set_yticks(range(6)); ax.set_yticklabels(lab, fontsize=8)
                cond = "0 K (ELASTIC)" if name == "elastic0" else "{:.0f} K (born/matrix)".format(r["T"])
                ax.set_title(f"C_ij (GPa), {cond}", loc="left", fontsize=10)
                png = f"prop_{name}.png"
            elif name in ("tensile", "shear"):
                r = mechanical(out_dir, m, dt_fs)
                label = (f"{'tension' if int(m.get('sign', 1)) > 0 else 'compression'} {m['dir']}" if name == "tensile"
                         else f"shear {m['plane']} ({m.get('mode', 'solid')})")
                item.update(label=label, **{k: v for k, v in r.items() if k not in ("strain", "stress")})
                sg = int(m.get("sign", 1)) if name == "tensile" else 1
                x, y = sg * r["strain"], sg * r["stress"]
                fig, ax = plt.subplots(figsize=(7.5, 4.8))
                ax.plot(x, y, lw=1.2)
                if "modulus" in r:
                    xx = np.array([0, r["fit_max"]])
                    ax.plot(xx, r["modulus"] * xx + (np.polyfit(x[(x >= 0) & (x <= r["fit_max"])], y[(x >= 0) & (x <= r["fit_max"])], 1)[1]),
                            color="#d62728", ls="--", lw=1, label=f"slope {r['modulus']:.1f} GPa (strain ≤ {r['fit_max']:g})")
                ax.plot([r["peak_strain"]], [r["peak"]], "o", color="#d62728", label=f"peak {r['peak']:.2f} GPa at {r['peak_strain']:.3f}")
                ax.set_xlabel("engineering shear strain" if name == "shear" else ("strain" if sg > 0 else "compressive strain"))
                ax.set_ylabel(("shear stress" if name == "shear" else ("stress" if sg > 0 else "compressive stress")) + " (GPa)")
                ax.legend(fontsize=8); ax.set_title(f"{label}: {m['file']}", loc="left", fontsize=10)
                png = f"prop_{os.path.splitext(m['file'])[0]}.png"
            else:
                continue
            fig.tight_layout(); fig.savefig(os.path.join(dest, png), dpi=100); plt.close(fig)
            item["png"] = png
            figs.append(png)
        except Exception as e:
            item["error"] = f"{type(e).__name__}: {e}"
            log(f"WARNING {name}: {item['error']}")
        res.append(item)
    return res


def summary_lines(items, ionic=None):
    L = []
    if ionic:
        L.append(f"イオン伝導度（ネルンスト–アインシュタイン：σ = e²/(V k_B T) Σ N z² D、イオン間の相関は含まない；T = {ionic['T']:.1f} K、"
                 f"V = {ionic['V']:.0f} Å³）：全体 {ionic['sigma_total']:.4g} S/m（面内 xy {ionic['sigma_xy']:.4g}、垂直 z {ionic['sigma_z']:.4g}）")
        for el, d in ionic["ions"].items():
            L.append(f"  {el}（{d['N']} 個、電荷 {d['q']:+.3f} e）：{d['sigma_total']:.4g} S/m（xy {d['sigma_xy']:.4g}）")
    for it in items:
        n = it["property"]
        if it.get("error"):
            L.append(f"{dict(ALL_PROPS).get(n, n).split('（')[0]}：解析できませんでした（{it['error']}）")
            continue
        if n in ("viscosity", "thermal"):
            what = "粘度" if n == "viscosity" else "熱伝導率"
            L.append(f"{what}（グリーン–久保、{it['outputs']} 回目の累積平均、相関 {it['corr_ps']:.1f} ps、T = {it['T']:.1f} K）："
                     f"{it['value']:.4g} {it['unit']}（積分の後半の平均。終点 {it['value_end']:.4g}）")
            L.append("  成分：" + "、".join(f"{k} {v:.4g}" for k, v in it["components"].items())
                     + ("。固体や層間水を含む系では粘度の積分は収束しません（液体のみ）" if n == "viscosity" else
                        "。積分が平らになっていることを図で確かめてください（長い MD が要ります）"))
            if n == "thermal":
                L.append("  水やイオンを含む多成分系では、熱流束の対流項（原子のエネルギー × 速度）に拡散の寄与が入るので、"
                         "値は参考です（固体の結晶で確かめた方法）")
        elif n in ("elastic0", "elasticT"):
            head = "弾性定数 0 K（ELASTIC）" if n == "elastic0" else f"弾性定数 {it['T']:.0f} K（born/matrix、応力 {it['samples_stress']} 点、Born {it['samples_born']} 点）"
            L.append(head + "：C_ij（GPa、Voigt：1 xx 2 yy 3 zz 4 yz 5 xz 6 xy）")
            for row in it["C"]:
                L.append("  " + " ".join(f"{v:8.1f}" if v == v else "     nan" for v in row))
            diag = [it["C"][k][k] for k in range(6) if it["C"][k][k] == it["C"][k][k]]
            if n == "elasticT" and any(d <= 0 for d in diag):
                L.append("  ⚠ 対角成分が負：ゆらぎの式はサンプル不足です（Born 項と応力のゆらぎがほぼ打ち消し合うので、"
                         "長い平衡 MD（数百 ps〜ns）が要ります）")
            mo = it.get("moduli") or {}
            if "K_Hill" in mo:
                L.append(f"  体積弾性率 K（Hill）{mo['K_Hill']:.1f} GPa、せん断弾性率 G（Hill）{mo['G_Hill']:.1f} GPa、"
                         f"ヤング率 E_x {mo['E_x']:.1f}、E_y {mo['E_y']:.1f}、E_z {mo['E_z']:.1f} GPa、ポアソン比 ν_xy {mo['nu_xy']:.3f}")
        elif n in ("tensile", "shear"):
            if n == "tensile":
                s = f"{it['label']}：ヤング率 {it['modulus']:.1f} GPa（ひずみ {it['fit_max']:g} までの傾き）、" if "modulus" in it else f"{it['label']}："
                s += f"{'引張' if int(it.get('sign', 1)) > 0 else '圧縮'}強度 {it['peak']:.2f} GPa（ひずみ {it['peak_strain']:.3f}）"
                s += "".join(f"、ポアソン比 ν_{k[-1]} {v:.3f}" for k, v in it.items() if k.startswith("poisson_"))
            else:
                s = f"{it['label']}：" + (f"せん断弾性率 {it['modulus']:.1f} GPa、" if "modulus" in it else "")
                s += f"せん断強度 {it['peak']:.2f} GPa（ひずみ {it['peak_strain']:.3f}）"
                if "viscosity_Pa_s" in it:
                    s += f"、粘度（NEMD、後半の平均）{it['viscosity_Pa_s'] * 1e3:.4g} mPa·s"
            s += f"  ひずみ速度 {float(it['rate_s']):.2g} /s、初期応力 {it['initial_stress']:.2f} GPa"
            if abs(it["initial_stress"]) > 0.1:
                s += "（⚠ 0 から離れています：平衡化（NPT）を長くしてください）"
            L.append(s)
    return L
