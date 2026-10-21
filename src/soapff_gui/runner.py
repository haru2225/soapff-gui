"""What the GUI does besides the window: the bundled soapff-cui (force field, analysis, OVITO), the library, LAMMPS and
mpirun on this computer, SOVA, and running LAMMPS."""
import glob
import os
import shutil
import signal
import subprocess
import sys

from soapff_cui import core

LIBRARY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "library_clayff_bulk.npz")


def import_soapff():
    """the force-field writer (soapff_cui.core, bundled)."""
    return core


def default_library():
    return LIBRARY


def which(path_or_name):
    """absolute path of an executable given by path or name, or None."""
    if not path_or_name:
        return None
    if os.path.isfile(path_or_name):
        return os.path.abspath(path_or_name)
    return shutil.which(path_or_name)


def find_executables():
    """(lmp, mpirun): first next to this Python, then ~/.soapff-gui/lammps/bin (install.sh), then on PATH; the mpirun next to lmp (same MPI) is preferred."""
    home = os.environ.get("SOAPFF_HOME") or os.path.join(os.path.expanduser("~"), ".soapff-gui")
    dirs = [os.path.dirname(sys.executable), os.path.join(home, "lammps", "bin")] + os.environ.get("PATH", "").split(os.pathsep)

    def first(names):
        for d in dirs:
            for n in names:
                p = os.path.join(d, n)
                if os.path.isfile(p) and os.access(p, os.X_OK):
                    return p
        return None
    lmp = first(["lmp", "lmp_mpi", "lmp_serial"])
    mpirun = None
    if lmp and not lmp.endswith("lmp_serial"):
        mpirun = next((os.path.join(os.path.dirname(lmp), n) for n in ("mpirun", "mpiexec")
                       if os.path.isfile(os.path.join(os.path.dirname(lmp), n))), None) or first(["mpirun", "mpiexec"])
    return lmp, mpirun


def prepare(structure, workdir, library=None, log=print, **md):
    """soapff-cui make: structure -> workdir/{system.data, system.settings, in.lammps, charges.csv, report.json}."""
    md.setdefault("stop_on_fallback", True)        # the GUI does not write a force field with fallback atoms
    return core.make(structure, workdir, library or LIBRARY, log=log, **md)


def find_debyer(path=None):
    """the debyer executable (given path, $DEBYER, PATH, ~/*/debyer/debyer/debyer) or None."""
    from soapff_cui import analysis
    return analysis.find_debyer(path or None)


def analyze(out_dir, log=print, skip_ps=None, sova_opts=None, what=None, debyer=None, xrd_opts=None,
            traj=None):
    """diffusion from LAMMPS compute msd, SOVA PDF / RDF / S(Q), the fix ave/* files of the computes, XRD by Debyer
    -> out_dir/analysis/."""
    from soapff_cui import analysis
    return analysis.analyze(out_dir, skip_ps=skip_ps, sova_opts=sova_opts, log=log,
                            what=tuple(what) if what else analysis.WHAT, debyer=debyer or None, xrd_opts=xrd_opts,
                            traj=traj or None)


def analyzable(out_dir):
    """the outputs that the analysis can use in out_dir: dumps (*.lammpstrj) and MSD files -> list of short descriptions (empty:
    nothing to analyse yet)."""
    from soapff_cui import analysis
    if not out_dir or not os.path.isdir(out_dir):
        return []
    found = [f for f in analysis.trajectories(out_dir)]
    found += [os.path.basename(f) for f in glob.glob(os.path.join(out_dir, "msd_*.dat"))]
    return found


def has_ovito_module():
    try:
        import ovito  # noqa: F401
        return True
    except ImportError:
        return False


def replicate(structure, n, out=None):
    """supercell of a structure file: OVITO's Replicate modifier (ovito Python module; ReplicateModifier with
    adjust_box, the same as Replicate in the OVITO application) -> extended XYZ (elements, positions, cell only);
    without the ovito module, ASE Atoms.repeat (the same supercell).  n = (nx, ny, nz).
    Returns (path written, number of atoms, cell lengths a b c, 'ovito' or 'ase')."""
    nx, ny, nz = (int(v) for v in n)
    if min(nx, ny, nz) < 1:
        raise ValueError("the repetitions must be 1 or more")
    stem, _ = os.path.splitext(os.path.basename(structure))
    out = out or os.path.join(os.path.dirname(os.path.abspath(structure)), f"{stem}_{nx}x{ny}x{nz}.extxyz")
    try:
        from ovito.io import import_file, export_file
        from ovito.modifiers import ReplicateModifier
        tool = "ovito"
    except ImportError:
        tool = "ase"
    at0 = core.read_structure(structure)                      # any format soapff reads (POSCAR, CIF, LAMMPS data ...)
    if tool == "ovito":
        import tempfile
        from ase.io import write
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.extxyz")
            write(src, at0, format="extxyz")
            pl = import_file(src)
            pl.modifiers.append(ReplicateModifier(num_x=nx, num_y=ny, num_z=nz, adjust_box=True, unique_ids=True))
            export_file(pl, out, "xyz", columns=["Particle Type", "Position.X", "Position.Y", "Position.Z"])
        from ase.io import read
        at = read(out)
    else:
        from ase.io import write
        at = at0.repeat((nx, ny, nz))
        write(out, at, format="extxyz")
    return out, len(at), tuple(float(x) for x in at.cell.lengths()), tool


def open_in_ovito(path):
    return core.open_in_ovito(path)


class LammpsJob:
    """lmp -in in.lammps in workdir (with mpirun when nproc > 1); lines() yields the screen output."""

    def __init__(self, workdir, lmp, mpirun=None, nproc=1):
        cmd = [mpirun, "-np", str(nproc), lmp] if nproc > 1 and mpirun else [lmp]
        self.cmd = cmd + ["-in", "in.lammps", "-log", "log.lammps"]
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.setdefault("FI_PROVIDER", "tcp")    # MPICH/libfabric on a laptop
        self.proc = subprocess.Popen(self.cmd, cwd=workdir, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, bufsize=1, start_new_session=True)

    def lines(self):
        for line in self.proc.stdout:
            yield line.rstrip("\n")
        self.proc.wait()

    def stop(self):
        if self.proc.poll() is None:
            try:
                os.killpg(self.proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    @property
    def returncode(self):
        return self.proc.poll()
