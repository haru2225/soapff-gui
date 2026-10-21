"""The fix / compute styles of LAMMPS (and a few other commands) that soapff-GUI offers to insert after the basic input
(force field + minimisation): title, description, a template line, the names of its arguments and the manual page.

Placeholders in templates: {ID} (made unique in the input), {T} temperature (K), {P} pressure (atm), {DT} time step (fs),
{N} steps.  The manual: https://docs.lammps.org/<doc>."""

DOC = "https://docs.lammps.org/"
BASE = ["ID", "Group", "Style"]

# kind, style, title, description, template, labels after BASE (None: fixed keyword, not a field), doc page
_FIX = [
    # --- time integration and thermostats / barostats
    ("Integration, temperature, pressure", "nvt", "NVT (Nosé-Hoover, constant temperature)",
     "Constant-temperature MD. A damping time of about 100 fs (100 times the time step) is a good start.",
     "fix {ID} all nvt temp {T} {T} 100.0", [None, "Start temperature (K)", "End temperature (K)", "Damping (fs)"], "fix_nh.html"),
    ("Integration, temperature, pressure", "npt", "NPT (isotropic)", "Constant temperature and pressure. The three box edges change by the same ratio.",
     "fix {ID} all npt temp {T} {T} 100.0 iso {P} {P} 1000.0",
     [None, "Start temperature (K)", "End temperature (K)", "Damping (fs)", "Direction", "Start pressure (atm)", "End pressure (atm)", "Damping (fs)"], "fix_nh.html"),
    ("Integration, temperature, pressure", "npt", "NPT (normal to the surface only: interlayers, interfaces)",
     "Changes only the z box length. For clay interlayers and interfaces.",
     "fix {ID} all npt temp {T} {T} 100.0 z {P} {P} 1000.0",
     [None, "Start temperature (K)", "End temperature (K)", "Damping (fs)", "Direction", "Start pressure (atm)", "End pressure (atm)", "Damping (fs)"], "fix_nh.html"),
    ("Integration, temperature, pressure", "npt", "NPT (three edges independently: bulk crystals)", "Changes x, y and z independently.",
     "fix {ID} all npt temp {T} {T} 100.0 aniso {P} {P} 1000.0",
     [None, "Start temperature (K)", "End temperature (K)", "Damping (fs)", "Direction", "Start pressure (atm)", "End pressure (atm)", "Damping (fs)"], "fix_nh.html"),
    ("Integration, temperature, pressure", "npt", "NPT (all components, including box tilt)",
     "Changes every component of a tilted box (a rectangular box first needs change_box all triclinic).",
     "fix {ID} all npt temp {T} {T} 100.0 tri {P} {P} 1000.0",
     [None, "Start temperature (K)", "End temperature (K)", "Damping (fs)", "Direction", "Start pressure (atm)", "End pressure (atm)", "Damping (fs)"], "fix_nh.html"),
    ("Integration, temperature, pressure", "nphug", "NPHug (Hugoniostat: Hugoniot state of shock compression)",
     "Drives the system to a state on the shock Hugoniot curve (pressure, volume and energy satisfy the Hugoniot condition). "
     "The pressure is the shock pressure (atm); the start and end temperatures are not used, only the damping matters. "
     "The reference state (E0, V0, P0) is the value at the first step (it can also be set with keywords). "
     "z is a uniaxial shock (only the z box length changes).",
     "fix {ID} all nphug temp {T} {T} 100.0 z 10000.0 10000.0 1000.0",
     [None, "Start temperature (unused)", "End temperature (unused)", "Damping (fs)", "Direction", "Start pressure (atm)", "End pressure (atm)", "Damping (fs)"], "fix_nphug.html"),
    ("Integration, temperature, pressure", "nph", "NPH (constant pressure, no temperature control)", "",
     "fix {ID} all nph iso {P} {P} 1000.0", ["Direction", "Start pressure (atm)", "End pressure (atm)", "Damping (fs)"], "fix_nh.html"),
    ("Integration, temperature, pressure", "nve", "NVE (constant energy)", "Time integration only. Sometimes combined with a thermostat fix (langevin etc.).",
     "fix {ID} all nve", [], "fix_nve.html"),
    ("Integration, temperature, pressure", "langevin", "Langevin thermostat", "Use together with fix nve. The random seed is a positive integer.",
     "fix {ID} all langevin {T} {T} 100.0 48279", ["Start temperature (K)", "End temperature (K)", "Damping (fs)", "Random seed"], "fix_langevin.html"),
    ("Integration, temperature, pressure", "temp/csvr", "CSVR thermostat (Bussi)", "Use together with fix nve.",
     "fix {ID} all temp/csvr {T} {T} 100.0 54324", ["Start temperature (K)", "End temperature (K)", "Damping (fs)", "Random seed"], "fix_temp_csvr.html"),
    ("Integration, temperature, pressure", "temp/berendsen", "Berendsen thermostat", "Use together with fix nve. Suited to equilibration.",
     "fix {ID} all temp/berendsen {T} {T} 100.0", ["Start temperature (K)", "End temperature (K)", "Damping (fs)"], "fix_temp_berendsen.html"),
    ("Integration, temperature, pressure", "temp/rescale", "Velocity rescaling", "Resets the temperature to the target every N steps. Use together with fix nve.",
     "fix {ID} all temp/rescale 100 {T} {T} 10.0 1.0",
     ["Interval (steps)", "Start temperature (K)", "End temperature (K)", "Window (K)", "Fraction"], "fix_temp_rescale.html"),
    ("Integration, temperature, pressure", "press/berendsen", "Berendsen barostat", "Use together with a thermostat fix.",
     "fix {ID} all press/berendsen iso {P} {P} 1000.0", ["Direction", "Start pressure (atm)", "End pressure (atm)", "Damping (fs)"],
     "fix_press_berendsen.html"),
    ("Integration, temperature, pressure", "box/relax", "Box minimization (place before minimize)", "Relaxes the box size too during energy minimization.",
     "fix {ID} all box/relax iso 0.0 vmax 0.001", ["Direction", "Target pressure (atm)", None, "Max change per iteration"], "fix_box_relax.html"),
    ("Integration, temperature, pressure", "rigid", "Rigid bodies (per molecule)", "Moves each molecule as a rigid body (temperature control: rigid/nvt etc.).",
     "fix {ID} all rigid molecule", ["Unit"], "fix_rigid.html"),
    # --- constraints and restraints
    ("Constraints and restraints", "shake", "SHAKE (fixed bonds and angles)",
     "After b and a come the bond / angle type numbers (bond_coeff / angle_coeff in system.settings). For making flexible SPC water rigid.",
     "fix {ID} water shake 1.0e-4 20 0 b 2 a 1", ["Tolerance", "Max iterations", "Output interval", None, "Bond type", None, "Angle type"],
     "fix_shake.html"),
    ("Constraints and restraints", "rattle", "RATTLE (SHAKE plus velocity constraints)", "", "fix {ID} water rattle 1.0e-4 20 0 b 2 a 1",
     ["Tolerance", "Max iterations", "Output interval", None, "Bond type", None, "Angle type"], "fix_shake.html"),
    ("Constraints and restraints", "setforce", "Fix the force (0 freezes the atoms)", "NULL components are left unchanged.",
     "fix {ID} all setforce 0.0 0.0 0.0", ["Fx", "Fy", "Fz"], "fix_setforce.html"),
    ("Constraints and restraints", "spring/self", "Spring to the initial position", "Tethers each atom to its starting position.",
     "fix {ID} all spring/self 10.0", ["Spring constant (kcal/mol/Å²)"], "fix_spring_self.html"),
    ("Constraints and restraints", "momentum", "Zero the total momentum", "Stops the translation (and rotation) of the system.",
     "fix {ID} all momentum 100 linear 1 1 1", ["Interval (steps)", None, "x", "y", "z"], "fix_momentum.html"),
    ("Constraints and restraints", "recenter", "Fix the center of mass", "", "fix {ID} all recenter INIT INIT INIT", ["x", "y", "z"], "fix_recenter.html"),
    ("Constraints and restraints", "restrain", "Restraints on bonds, angles, etc.", "Atom indices and a force constant.",
     "fix {ID} all restrain bond 1 2 1000.0 1000.0 1.0", [None, "Atom 1", "Atom 2", "Start strength", "End strength", "Distance (Å)"],
     "fix_restrain.html"),
    ("Constraints and restraints", "viscous", "Viscous drag", "", "fix {ID} all viscous 0.1", ["Coefficient"], "fix_viscous.html"),
    # --- external fields and walls
    ("External fields and walls", "efield", "Electric field", "V/Å in real units.", "fix {ID} all efield 0.0 0.0 0.1", ["Ex", "Ey", "Ez"], "fix_efield.html"),
    ("External fields and walls", "addforce", "Add a force", "", "fix {ID} all addforce 0.0 0.0 -0.1", ["Fx", "Fy", "Fz"], "fix_addforce.html"),
    ("External fields and walls", "aveforce", "Average the force over the group", "", "fix {ID} all aveforce NULL NULL 0.0", ["Fx", "Fy", "Fz"], "fix_aveforce.html"),
    ("External fields and walls", "wall/reflect", "Reflecting wall", "Keeps atoms from being lost through non-periodic faces.", "fix {ID} all wall/reflect zlo EDGE zhi EDGE",
     ["Face", "Position", "Face", "Position"], "fix_wall_reflect.html"),
    ("External fields and walls", "wall/lj93", "LJ 9-3 wall", "", "fix {ID} all wall/lj93 zlo EDGE 1.0 1.0 2.5",
     ["Face", "Position", "epsilon", "sigma", "Cutoff"], "fix_wall.html"),
    ("External fields and walls", "deform", "Deform the box", "E.g. tension at a constant strain rate.", "fix {ID} all deform 1 z erate 1.0e-6 remap x",
     ["Interval", "Direction", "Method", "Value", None, "remap"], "fix_deform.html"),
    ("External fields and walls", "heat", "Add or remove heat", "", "fix {ID} all heat 10 0.1", ["Interval", "Energy (kcal/mol/fs)"], "fix_heat.html"),
    # --- output / averaging
    ("Output and averaging", "ave/time", "Time-average and write to a file", "Writes a compute value (c_ID) every Nfreq steps.",
     "fix {ID} all ave/time 10 100 1000 c_thermo_temp file ave_time.dat", ["Nevery", "Nrepeat", "Nfreq", "Value", None, "File"],
     "fix_ave_time.html"),
    ("Output and averaging", "ave/chunk", "Per-chunk average (density profiles etc.)",
     "Place compute chunk/atom (ID cz) first. E.g. the mass density profile along z.",
     "fix {ID} all ave/chunk 10 100 1000 cz density/mass file density_z.dat",
     ["Nevery", "Nrepeat", "Nfreq", "Chunk ID", "Value", None, "File"], "fix_ave_chunk.html"),
    ("Output and averaging", "ave/correlate", "Time correlation (autocorrelation etc.)", "", "fix {ID} all ave/correlate 10 100 1000 c_thermo_press[4] type auto file corr.dat",
     ["Nevery", "Nrepeat", "Nfreq", "Value", None, "Correlation type", None, "File"], "fix_ave_correlate.html"),
    ("Output and averaging", "ave/atom", "Per-atom time average", "", "fix {ID} all ave/atom 10 100 1000 x y z",
     ["Nevery", "Nrepeat", "Nfreq", "Value", "Value", "Value"], "fix_ave_atom.html"),
    ("Output and averaging", "dt/reset", "Change the time step automatically", "", "fix {ID} all dt/reset 10 NULL 1.0 0.1 units box",
     ["Interval", "Min (fs)", "Max (fs)", "Max displacement (Å)"], "fix_dt_reset.html"),
]

_COMPUTE = [
    ("Diffusion and structure", "msd", "Mean squared displacement (MSD)", "Write it out with fix ave/time (c_ID[4] is the total).",
     "compute {ID} all msd", [], "compute_msd.html"),
    ("Diffusion and structure", "msd/chunk", "MSD per chunk", "Place compute chunk/atom first.", "compute {ID} all msd/chunk cmol",
     ["Chunk ID"], "compute_msd_chunk.html"),
    ("Diffusion and structure", "vacf", "Velocity autocorrelation", "", "compute {ID} all vacf", [], "compute_vacf.html"),
    ("Diffusion and structure", "rdf", "Radial distribution function (RDF)", "Write it out with fix ave/time ... mode vector.",
     "compute {ID} all rdf 100", ["Number of bins"], "compute_rdf.html"),
    ("Diffusion and structure", "displace/atom", "Per-atom displacement", "", "compute {ID} all displace/atom", [], "compute_displace_atom.html"),
    ("Diffusion and structure", "coord/atom", "Coordination number", "", "compute {ID} all coord/atom cutoff 3.0", [None, "Cutoff (Å)"],
     "compute_coord_atom.html"),
    ("Diffusion and structure", "cluster/atom", "Clusters (connectivity)", "", "compute {ID} all cluster/atom 3.5", ["Cutoff (Å)"],
     "compute_cluster_atom.html"),
    ("Diffusion and structure", "centro/atom", "Centrosymmetry parameter", "", "compute {ID} all centro/atom 12", ["Number of neighbors"],
     "compute_centro_atom.html"),
    ("Diffusion and structure", "orientorder/atom", "Orientational order (Steinhardt)", "", "compute {ID} all orientorder/atom", [],
     "compute_orientorder_atom.html"),
    ("Diffusion and structure", "gyration", "Radius of gyration", "", "compute {ID} all gyration", [], "compute_gyration.html"),
    ("Diffusion and structure", "com", "Center of mass", "", "compute {ID} all com", [], "compute_com.html"),
    ("Diffusion and structure", "chunk/atom", "Divide atoms into chunks (z bins etc.)", "Combine with fix ave/chunk or msd/chunk.",
     "compute {ID} all chunk/atom bin/1d z lower 0.5 units box", ["Method", "Direction", "Origin", "Width (Å)", None, "Units"],
     "compute_chunk_atom.html"),
    ("Temperature, energy, pressure", "temp", "Temperature", "", "compute {ID} all temp", [], "compute_temp.html"),
    ("Temperature, energy, pressure", "temp/com", "Temperature without the center-of-mass motion", "", "compute {ID} all temp/com", [], "compute_temp_com.html"),
    ("Temperature, energy, pressure", "temp/partial", "Temperature of selected directions only", "", "compute {ID} all temp/partial 1 1 0", ["x", "y", "z"],
     "compute_temp_partial.html"),
    ("Temperature, energy, pressure", "pe", "Potential energy", "", "compute {ID} all pe", [], "compute_pe.html"),
    ("Temperature, energy, pressure", "pe/atom", "Per-atom potential energy", "", "compute {ID} all pe/atom", [], "compute_pe_atom.html"),
    ("Temperature, energy, pressure", "ke", "Kinetic energy", "", "compute {ID} all ke", [], "compute_ke.html"),
    ("Temperature, energy, pressure", "ke/atom", "Per-atom kinetic energy", "", "compute {ID} all ke/atom", [], "compute_ke_atom.html"),
    ("Temperature, energy, pressure", "pressure", "Pressure tensor", "", "compute {ID} all pressure thermo_temp", ["Temperature compute"],
     "compute_pressure.html"),
    ("Temperature, energy, pressure", "stress/atom", "Per-atom stress", "", "compute {ID} all stress/atom NULL", ["Temperature compute"],
     "compute_stress_atom.html"),
    ("Temperature, energy, pressure", "heat/flux", "Heat flux", "Place ke/atom, pe/atom and stress/atom first.",
     "compute {ID} all heat/flux ke pe stress", ["ke/atom ID", "pe/atom ID", "stress/atom ID"], "compute_heat_flux.html"),
    ("Temperature, energy, pressure", "group/group", "Interaction between groups", "", "compute {ID} water group/group all", ["Other group"],
     "compute_group_group.html"),
    ("Temperature, energy, pressure", "pair", "Energy per pair-style component", "", "compute {ID} all pair lj/cut/coul/long", ["Pair style"],
     "compute_pair.html"),
    ("Atom and bond quantities", "property/atom", "Atom quantities (charge, velocity, etc.)", "", "compute {ID} all property/atom q", ["Quantity"],
     "compute_property_atom.html"),
    ("Atom and bond quantities", "reduce", "Sum or average of per-atom quantities", "", "compute {ID} all reduce sum c_pea", ["Operation", "Value"], "compute_reduce.html"),
    ("Atom and bond quantities", "count/type", "Number of atoms per type", "", "compute {ID} all count/type atom", ["Target"], "compute_count_type.html"),
    ("Atom and bond quantities", "bond", "Bond energy", "", "compute {ID} all bond", [], "compute_bond.html"),
    ("Atom and bond quantities", "angle", "Angle energy", "", "compute {ID} all angle", [], "compute_angle.html"),
    ("Atom and bond quantities", "bond/local", "Per-bond length", "", "compute {ID} all bond/local dist", ["Quantity"], "compute_bond_local.html"),
    ("Atom and bond quantities", "angle/local", "Per-angle value", "", "compute {ID} all angle/local theta", ["Quantity"], "compute_angle_local.html"),
    ("Atom and bond quantities", "dipole", "Dipole moment", "", "compute {ID} all dipole", [], "compute_dipole.html"),
]

_OTHER = [
    ("MD basics", "run", "Run MD", "{N} steps. Time = steps x time step.", "run {N}", ["Steps"], "run.html"),
    ("MD basics", "timestep", "Time step", "0.5 fs for flexible SPC water (ClayFF).", "timestep {DT}", ["Time step (fs)"], "timestep.html"),
    ("MD basics", "velocity", "Assign initial velocities", "Maxwell distribution at {T} K.", "velocity all create {T} 4928459 mom yes rot no dist gaussian",
     ["Group", "Method", "Temperature (K)", "Random seed"], "velocity.html"),
    ("MD basics", "unfix", "Remove a fix", "", "unfix {LASTFIX}", ["ID"], "unfix.html"),
    ("MD basics", "uncompute", "Remove a compute", "", "uncompute {LASTCOMPUTE}", ["ID"], "uncompute.html"),
    ("MD basics", "reset_timestep", "Reset the step to 0", "", "reset_timestep 0", ["Step"], "reset_timestep.html"),
    ("MD basics", "minimize", "Energy minimization", "", "minimize 1.0e-6 1.0e-6 5000 50000",
     ["Energy tolerance", "Force tolerance", "Max iterations", "Max evaluations"], "minimize.html"),
    ("Output", "thermo", "Screen output interval", "", "thermo 1000", ["Interval (steps)"], "thermo.html"),
    ("Output", "thermo_style", "Screen output columns", "", "thermo_style custom step temp press pe ke etotal vol density lx ly lz",
     ["Style"], "thermo_style.html"),
    ("Output", "dump", "Trajectory (with type names and charges, for OVITO)", "Writes the ClayFF type name in the element column.",
     "dump {ID} all custom 2000 traj.lammpstrj id type element x y z ix iy iz q", ["ID", "Group", "Style", "Interval (steps)", "File"],
     "dump.html"),
    ("Output", "write_data", "Write the final structure", "", "write_data final.data", ["File"], "write_data.html"),
    ("Output", "write_restart", "Restart file", "", "write_restart restart.bin", ["File"], "write_restart.html"),
    ("Groups", "group", "Create a group by type", "Type numbers are in the Masses section of system.data.", "group {ID} type 1", ["Group name", "Method", "Type"],
     "group.html"),
]

# special entries built by the application (need the type numbers of the system)
SPECIAL = [
    ("soapff workflows", "soapff:diffusion", "MSD for the diffusion coefficient (for soapff analysis)",
     "Creates groups for water O and each ion, and writes compute msd with fix ave/time to msd_<name>.dat every 0.1 ps. "
     "\"Analyze\" derives the diffusion coefficient from these files. Place it before run."),
    ("soapff workflows", "soapff:protocol", "soapff relaxation (initial velocities -> relaxation NVT)",
     "Inserts the initial velocities, the time step, a relaxation NVT (not recorded) and the step reset in one go. "
     "The production part (fix nvt / npt, dump, diffusion MSD, run, write_data) is added afterwards by choosing items from the list. "
     "The length and temperature are the MD conditions in the calculation settings."),
    ("soapff workflows", "soapff:swelling", "Swelling curve (add water to the interlayer under NPT, record d001)",
     "Method of Cygan et al. (2004): add water molecules to the interlayer one at a time (later two at a time); at each state run "
     "a minimization and NPT (average after equilibration) and write d001 (= lz) to swelling.dat. A dehydration pass that removes "
     "two at a time can also be selected at the end. If the structure has no water type, one is added to system.data "
     "(water.mol is written too). The clay is kept at the center of the box, and the interlayer is the side across the periodic boundary. "
     "The columns of swelling.dat are step, number of added waters, phase (1 = recorded), lx ly lz xy xz yz, PE, T "
     "(the example analyze_swelling.py can plot d001 against the water content). Place it right after the minimization in in.lammps."),
]


def entries(kind):
    """list of dicts: kind ('fix', 'compute', 'other', 'special'), category, style, title, desc, template, labels, doc."""
    src = {"fix": _FIX, "compute": _COMPUTE, "other": _OTHER}.get(kind)
    if kind == "special":
        return [dict(kind="special", category=c, style=s, title=t, desc=d, template="", labels=[], doc="") for c, s, t, d in SPECIAL]
    return [dict(kind=kind, category=c, style=s, title=t, desc=d, template=tpl, labels=lab, doc=DOC + doc)
            for c, s, t, d, tpl, lab, doc in src]


def labels_for(command, args):
    """field labels of a fix / compute line from the catalog (the entry with the same style and, when there are several,
    the same keyword arguments); None when the style is not in the catalog."""
    if command not in ("fix", "compute") or len(args) < 3:
        return None
    style = args[2]
    cands = [e for e in entries(command) if e["style"] == style]
    if not cands:
        return None
    best = cands[0]
    for e in cands:                     # e.g. npt iso / z / aniso / tri: same labels, first one is fine
        toks = e["template"].split()[4:]
        if all(i + 3 < len(args) and (t.startswith("{") or t == args[i + 3]) for i, t in enumerate(toks[:1])):
            best = e
            break
    return BASE + best["labels"]
