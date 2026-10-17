"""Everything the GUI does, without the GUI: find soapff-cui / LAMMPS, write the input with soapff-cui, run LAMMPS and read
its screen output (phase, thermo values, progress)."""
import os
import shutil
import signal
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLED_LIBRARY = os.path.join(HERE, "data", "library_clayff.npz")
PACKAGES = ("soapff_cui", "soapff")     # soapff was renamed soapff-cui (package soapff_cui); the old name still works


def _soapff_homes():
    """checkouts to look in when soapff-cui is not installed: $SOAPFF_HOME, then next to this folder."""
    up = os.path.join(HERE, "..", "..")
    return [d for d in (os.environ.get("SOAPFF_HOME"), os.path.join(up, "soapff-cui"), os.path.join(up, "soapff")) if d]


def import_soapff():
    """soapff_cui.core (or the old soapff.core); when not installed, look for a checkout (see _soapff_homes)."""
    import importlib
    for pkg in PACKAGES:
        try:
            return importlib.import_module(pkg + ".core")
        except ImportError:
            pass
    for d in _soapff_homes():
        for pkg in PACKAGES:
            if os.path.isfile(os.path.join(d, pkg, "core.py")):
                sys.path.insert(0, os.path.abspath(d))
                return importlib.import_module(pkg + ".core")
    raise ImportError("soapff-cui not found: pip install -e /path/to/soapff-cui, or set SOAPFF_HOME to its folder")


def default_library():
    """the library of the soapff-cui checkout (examples/minerals/library_clayff.npz, the newest), else the copy
    bundled with the GUI."""
    try:
        core = import_soapff()
        p = os.path.join(os.path.dirname(os.path.dirname(core.__file__)), "examples", "minerals", "library_clayff.npz")
        if os.path.isfile(p):
            return os.path.abspath(p)
    except ImportError:
        pass
    return BUNDLED_LIBRARY



def find_executables():
    """(lmp, mpirun) found first in the Python environment's bin, then on PATH; None when missing."""
    dirs = [os.path.dirname(sys.executable)] + os.environ.get("PATH", "").split(os.pathsep)
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
        # the mpirun next to lmp belongs to the MPI it was built with
        for n in ("mpirun", "mpiexec"):
            p = os.path.join(os.path.dirname(lmp), n)
            if os.path.isfile(p):
                mpirun = p
                break
        mpirun = mpirun or first(["mpirun", "mpiexec"])
    return lmp, mpirun


def structure_summary(path):
    """a short description of a structure file (atoms per element, cell lengths) to show before running."""
    core = import_soapff()
    at = core.read_structure(path)
    el = at.get_chemical_symbols()
    counts = {e: el.count(e) for e in sorted(set(el))}
    a, b, c = at.cell.lengths()
    return dict(atoms=len(at), elements=counts, cell=(a, b, c), angles=tuple(at.cell.angles()))


def prepare(structure, workdir, library=None, log=print, **md):
    """soapff make: structure -> workdir/{system.data, system.settings, in.lammps, charges.csv, report.json}.
    md: temp, nvt_ps, npt_ps, press, baro, threshold, cutoff (soapff.core.make keywords)."""
    core = import_soapff()
    return core.make(structure, workdir, library or default_library(), log=log, **md)


def open_in_ovito(path):
    """soapff-cui's open_in_ovito: an output folder (traj.lammpstrj, else final.data / system.data) or a file in OVITO."""
    core = import_soapff()
    if not hasattr(core, "open_in_ovito"):
        raise RuntimeError("this soapff has no open_in_ovito; use soapff-cui")
    return core.open_in_ovito(path)


def total_md_steps(nvt_ps, npt_ps):
    """number of MD steps in.lammps runs (1 fs time step; soapff rounds each run to an integer)."""
    return int(round(nvt_ps * 1000)) * (nvt_ps > 0) + int(round(npt_ps * 1000)) * (npt_ps > 0)


class LammpsJob:
    """lmp -in in.lammps in workdir, with mpirun when nproc > 1.  lines() yields the screen output line by line."""

    def __init__(self, workdir, lmp, mpirun=None, nproc=1):
        if nproc > 1 and mpirun:
            cmd = [mpirun, "-np", str(nproc), lmp]
        else:
            cmd = [lmp]
        self.cmd = cmd + ["-in", "in.lammps", "-log", "log.lammps"]
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.setdefault("FI_PROVIDER", "tcp")    # MPICH/libfabric on a laptop (as in soapff's run_example.sh)
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


class ThermoParser:
    """follows the LAMMPS screen output: phase (minimisation / MD), thermo columns and MD progress.
    soapff's in.lammps resets the step to 0 after minimisation, and the NPT run continues the NVT step count."""

    def __init__(self, md_steps):
        self.md_steps = max(md_steps, 1)
        self.phase = "start"
        self.columns = None
        self.values = {}

    def feed(self, line):
        """returns True when phase or values changed."""
        s = line.strip()
        if s.startswith("Setting up cg style minimization"):
            self.phase, self.columns, self.values = "minimize", None, {}
            return True
        if s.startswith("Setting up Verlet run"):
            self.phase, self.columns, self.values = "md", None, {}
            return True
        if s.startswith("Loop time of"):
            self.columns = None
            return False
        if s.startswith("Step") and len(s.split()) > 1:
            self.columns = s.split()
            return False
        if self.columns:
            parts = s.split()
            if len(parts) == len(self.columns):
                try:
                    vals = [float(x) for x in parts]
                except ValueError:
                    return False
                self.values = dict(zip(self.columns, vals))
                return True
        return False

    @property
    def progress(self):
        """0..1 over the MD steps (0 during minimisation)."""
        if self.phase != "md" or "Step" not in self.values:
            return 0.0
        return min(self.values["Step"] / self.md_steps, 1.0)


def which(path_or_name):
    """absolute path of an executable given by path or name, or None."""
    if not path_or_name:
        return None
    if os.path.isfile(path_or_name):
        return os.path.abspath(path_or_name)
    return shutil.which(path_or_name)
