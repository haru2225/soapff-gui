"""soapff-GUI: choose a clay structure; soapff-cui writes the LAMMPS input (ClayFF + SOAP charges), which is shown line by
line and can be edited (select a line, change its values below); press "Start calculation" to run LAMMPS on this computer or on
another one over SSH, then analyse (diffusion coefficients; SOVA PDF, RDF, S(Q); the results of the computes)
and look at the result in OVITO.

Layout (design/soapff-gui-ui-concept.png): menu row, toolbar (New calculation / Calculation settings › Run & analysis / Start calculation), project sidebar
(New calculation, Calculation history, Structure library, Recent calculations), the input file (in.lammps) as an editor, the form of the selected
command, the run panel (status, progress, log) and a status bar."""
import json
import math
import os
import queue
import subprocess
import sys
import threading
import time
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, scrolledtext, simpledialog, ttk

import webbrowser

from . import __version__, catalog, lammps_input as L, remote, runner

SETTINGS = os.path.join(os.path.expanduser("~"), ".config", "soapff-gui", "settings.json")
STRUCTURE_TYPES = [("Structure files", "*.xyz *.extxyz *.cif *.vasp *.data *.lammps *.lmp POSCAR* CONTCAR*"), ("All files", "*")]
AXIS_STYLES = ["p", "f", "m", "s", "fm", "mf", "fs", "sf", "ms", "sm"]
BARO = {"z": "z: only the direction normal to the surface (interlayer / interface)", "aniso": "aniso: the 3 edges independently (bulk crystal)",
        "tri": "tri: everything, including the box tilt (bulk crystal)", "iso": "iso: isotropic"}
DEFAULTS = dict(temp=300.0, relax_ps=10.0, nvt_ps=50.0, npt_ps=0.0, press=1.0, baro="aniso", boundary="p p p",
                nproc=max(1, min(3, (os.cpu_count() or 2) // 2)), threshold=0.95, cutoff=10.0, library="", lmp="", mpirun="",
                charges="clayff2004", water="spc",
                last_dir=os.path.expanduser("~"), dump_ps=1.0, auto_analyze=True, skip_ps="",
                analyses=["msd", "sova"], debyer="", xrd_wavelength=1.5406, xrd_from=3.0, xrd_to=80.0,
                xrd_step=0.02, xrd_frames=3, xrd_mode="cluster", xrd_size=50.0,
                where="local", computers=[], history=[],
                sova=dict(dr=0.02, rmax=8.0, qmin=0.05, qmax=12.0, dq=0.05, radiation="xray", max_frames=40, pairs="water_ion",
                          all_elements=True),
                swelling=dict(one_to=23, two_to=99, down=True, eq_ps=30.0, avg_ps=20.0))
# charge schemes and water models (soapff_cui.clayff CHARGE_SCHEMES / WATER_MODELS): key -> label
CHARGES = {"clayff2004": "ClayFF 2004 (Cygan et al., Table 1)"}
WATERS = {"spc": "Flexible SPC (ClayFF)", "spce": "Rigid SPC/E (SHAKE)"}
# the analyses (soapff_cui.analysis): key, label
ANALYSES = (("msd", "Diffusion coefficient (when LAMMPS has a compute msd)"), ("sova", "SOVA: PDF g(r), RDF N(r), coordination numbers, S(Q)"))
XRD_MODES = {"cluster": "Crystallite: a finite crystal of tiled cells (crystals, clay minerals)",
             "periodic": "Periodic boundary: cutoff = half the box (liquids, glasses; orthorhombic boxes only)"}
# a registered computer (SSH): name, host, user, port, root (folder for the runs there), lmp, mpirun, setup (shell commands
# run first, e.g. 'module load lammps'), nproc.  The password is never stored: it is asked every time.
COMPUTER = dict(name="", host="", user="", port=22, root="soapff_runs", lmp="lmp", mpirun="mpirun", setup="", nproc=4)

# the name of the local computer in the "Run on" list (stored in settings as where="local")
LOCAL_NAME = "This computer"
LEGACY_LOCAL_NAME = "このコンピュータ"        # backward compatibility: the name used by earlier (Japanese) versions
# status words stored in the calculation history (settings.json)
STATUS_INPUT, STATUS_RUNNING, STATUS_DONE, STATUS_ERROR, STATUS_STOPPED = "Input created", "Running", "Done", "Error", "Stopped"
# backward compatibility: the status words written by earlier (Japanese) versions
LEGACY_STATUS = {"入力作成": STATUS_INPUT, "計算中": STATUS_RUNNING, "完了": STATUS_DONE, "エラー": STATUS_ERROR,
                 "停止": STATUS_STOPPED}
LOG_PLACEHOLDER = "The run log appears here."


def is_local_name(name):
    """True for the name of the local computer (also the Japanese one that may be in saved settings / history)."""
    return name in (LOCAL_NAME, LEGACY_LOCAL_NAME)


def status_of(h):
    """the status of a history entry, with the Japanese words of earlier versions mapped to the English ones."""
    st = h.get("status", "")
    return LEGACY_STATUS.get(st, st)


def where_of(h):
    """the 'run on' label of a history entry (earlier versions stored the Japanese name and full-width parentheses)."""
    w = h.get("where", "")
    return LOCAL_NAME if is_local_name(w) else w.replace("（", " (").replace("）", ")")

# colours and fonts of the design
BG, PANEL, SIDE, LINE = "#e9edf0", "#ffffff", "#f4f6f7", "#d5dbe0"
TEXT, MUTED, FAINT = "#35434c", "#77848d", "#a2abb1"
ACCENT, ACCENT_DARK, ACCENT_BG, CMD = "#28788d", "#1f6070", "#e3eef1", "#7257a3"
SEL_BG, CODE_BG, OK, ERR = "#e7f0f3", "#f6f7f8", "#3b8f5a", "#b3261e"


def load_settings():
    s = dict(DEFAULTS)
    try:
        s.update(json.load(open(SETTINGS)))
    except (OSError, ValueError):
        pass
    lib = s.get("library") or ""
    if not os.path.isfile(lib):
        s["library"] = runner.default_library()
    elif (os.path.basename(lib) == os.path.basename(runner.default_library())
          and os.path.abspath(lib) != os.path.abspath(runner.default_library())):
        # an older copy of the bundled library (it has fewer structures: e.g. no micas, so a mica gets fallback atoms):
        # the bundled one is used; a library with another name stays as the user chose it
        s["library"] = runner.default_library()
    if s.get("charges") not in CHARGES:              # a charge scheme of an earlier version (Underwood & Bourg, library)
        s["charges"] = "clayff2004"
    lmp, mpirun = runner.find_executables()
    s["sova"] = dict(DEFAULTS["sova"], **(s.get("sova") or {}))
    for k in ("sova_python", "sova_home"):          # the separate SOVA Python of earlier versions: sovapy is installed here now
        s.pop(k, None)
    if not s.get("debyer"):
        s["debyer"] = runner.find_debyer() or ""
    # analyses added in a later version are on: those not known when the settings were saved
    known = s.get("analyses_known") or ["msd", "sova", "compute", "xrd"]
    s["analyses"] = list(s.get("analyses") or []) + [k for k, _ in ANALYSES if k not in known and k not in (s.get("analyses") or [])]
    s["analyses"] = [k for k in s["analyses"] if k in dict(ANALYSES)]      # e.g. "props" of an earlier version
    s["analyses_known"] = [k for k, _ in ANALYSES]
    for k in ("props_selected", "props_params"):                          # settings of the removed "+ Properties" feature
        s.pop(k, None)
    s["lmp"] = s.get("lmp") or lmp or ""
    s["mpirun"] = s.get("mpirun") or mpirun or ""
    s["history"] = [h for h in s.get("history", []) if isinstance(h, dict)][:50]
    s["computers"] = [dict(COMPUTER, **c) for c in s.get("computers", []) if isinstance(c, dict) and c.get("name")]
    if s.get("r_host") and not s["computers"]:            # settings of an earlier version
        s["computers"].append(dict(COMPUTER, name=s["r_host"], host=s["r_host"], user=s.get("r_user", ""),
                                   port=s.get("r_port", 22)))
    for k in [k for k in s if k.startswith("r_")]:
        s.pop(k)
    if is_local_name(s.get("where")) or (s.get("where") != "local" and not any(c["name"] == s.get("where") for c in s["computers"])):
        s["where"] = "local"
    return s


def save_settings(s):
    try:
        os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
        json.dump(s, open(SETTINGS, "w"), indent=1, ensure_ascii=False)
    except OSError:
        pass


def open_folder(path):
    if sys.platform == "darwin":
        subprocess.Popen(["open", path])
    elif sys.platform.startswith("win"):
        os.startfile(path)
    else:
        subprocess.Popen(["xdg-open", path])


def ago(t):
    d = time.time() - float(t or 0)
    if d < 60:
        return "just now"
    if d < 3600:
        return f"{int(d // 60)} min ago"
    if d < 86400:
        return f"{int(d // 3600)} h ago"
    if d < 2 * 86400:
        return "yesterday"
    return f"{int(d // 86400)} d ago"


def pick_font(cands, size, weight="normal"):
    fams = set(tkfont.families())
    fam = next((c for c in cands if c in fams), "TkDefaultFont")
    return (fam, size, weight)


FILES = ("in.lammps", "system.settings", "system.data")
FILE_NOTE = {"in.lammps": "LAMMPS input", "system.settings": "Force field (included by in.lammps)",
             "system.data": "Atoms, bonds, angles (read-only)"}
READONLY = {"system.data"}
# columns of the sections of a LAMMPS data file (atom_style full)
DATA_COLS = {"Masses": ["Type no.", "Mass"], "Atoms": ["Atom no.", "Molecule no.", "Type", "Charge (e)", "x (Å)", "y (Å)", "z (Å)", "ix", "iy", "iz"],
             "Bonds": ["Bond no.", "Bond type", "Atom 1", "Atom 2"], "Angles": ["Angle no.", "Angle type", "Atom 1", "Atom 2 (center)", "Atom 3"],
             "Velocities": ["Atom no.", "vx", "vy", "vz"]}


class App:
    # the editor shows one of FILES; self.lines / self.undo are those of the shown file
    lines = property(lambda self: self.docs[self.tab], lambda self, v: self.docs.__setitem__(self.tab, v))
    undo = property(lambda self: self.undos[self.tab], lambda self, v: self.undos.__setitem__(self.tab, v))

    def __init__(self, root):
        self.root = root
        self.s = load_settings()
        self.q = queue.Queue()
        self.job = None            # local LammpsJob
        self.rjob = None           # (Remote, name, pid) of a remote run
        self.worker = None
        self.stopping = False
        self.structure = ""        # current project
        self.out = ""
        self.report = None
        self.docs = {f: [] for f in FILES}      # the lines of in.lammps and system.settings
        self.undos = {f: [] for f in FILES}
        self.sels = {f: None for f in FILES}
        self.dirty = set()
        self.tab = FILES[0]
        self.sel = None            # selected line index (of the shown file)
        self.form_vars = []
        self._form_cmd = ""
        self._save_after = None
        ui = ["Hiragino Sans", "Hiragino Kaku Gothic ProN", "Yu Gothic UI", "Noto Sans CJK JP", "Helvetica"]
        # monospaced fonts that also have Japanese glyphs first: with Menlo the Japanese of a comment falls back to another
        # font and that line (and the selected line with its hint) gets another height
        mono = ["BIZ UDGothic", "Osaka-Mono", "Noto Sans Mono CJK JP", "Sarasa Mono J", "Menlo", "SF Mono", "Consolas",
                "DejaVu Sans Mono", "Courier"]
        self.f_ui, self.f_small, self.f_bold = pick_font(ui, 12), pick_font(ui, 11), pick_font(ui, 12, "bold")
        self.f_title, self.f_head = pick_font(ui, 17, "bold"), pick_font(ui, 11, "bold")
        self.f_mono, self.f_mono_s = pick_font(mono, 12), pick_font(mono, 11)
        root.title("soapff-GUI — New calculation")
        root.geometry("1360x860")
        root.minsize(1100, 780)                    # smaller than this the editor would be squeezed out
        root.minsize(1100, 720)
        root.configure(bg=BG)
        self._styles()
        self._build()
        root.protocol("WM_DELETE_WINDOW", self.on_close)
        root.bind_all("<Command-r>" if sys.platform == "darwin" else "<Control-r>", lambda e: self.start())
        root.bind_all("<Command-n>" if sys.platform == "darwin" else "<Control-n>", lambda e: self.new_calc())
        root.bind_all("<Command-z>" if sys.platform == "darwin" else "<Control-z>", lambda e: self.do_undo())
        root.after(100, self.poll)
        self.refresh_recent()
        self.show_tab("in.lammps")
        self.update_status()

    # ================================================================================================== styles
    def _styles(self):
        st = ttk.Style(self.root)
        st.theme_use("clam")
        st.configure(".", font=self.f_ui, background=PANEL, foreground=TEXT)
        st.configure("Accent.TButton", background=ACCENT, foreground="white", borderwidth=0, focusthickness=0,
                     padding=(18, 8), font=self.f_bold)
        st.map("Accent.TButton", background=[("disabled", "#9bb8c0"), ("active", ACCENT_DARK)],
               foreground=[("disabled", "#eef3f5")])
        st.configure("Small.Accent.TButton", padding=(14, 4), font=pick_font(["Hiragino Sans", "Helvetica"], 11, "bold"))
        st.configure("Ghost.TButton", background=PANEL, foreground=TEXT, bordercolor="#c9d1d7", lightcolor=PANEL,
                     darkcolor=PANEL, borderwidth=1, padding=(16, 7))
        st.map("Ghost.TButton", background=[("active", "#f1f4f6")])
        st.configure("Tiny.TButton", background=PANEL, foreground=TEXT, bordercolor="#c9d1d7", lightcolor=PANEL,
                     darkcolor=PANEL, borderwidth=1, padding=(6, 3))
        st.map("Tiny.TButton", background=[("active", "#f1f4f6")])
        st.configure("Link.TButton", background=PANEL, foreground=ACCENT, borderwidth=1, bordercolor="#c9d1d7",
                     lightcolor=PANEL, darkcolor=PANEL, padding=(12, 5), font=self.f_small)
        st.map("Link.TButton", background=[("active", "#f1f4f6")])
        st.configure("Field.TEntry", fieldbackground=PANEL, bordercolor="#c9d1d7", lightcolor="#c9d1d7", darkcolor="#c9d1d7",
                     padding=6)
        st.configure("Accent.Horizontal.TProgressbar", troughcolor="#e3e8eb", background=ACCENT, bordercolor="#e3e8eb",
                     lightcolor=ACCENT, darkcolor=ACCENT, thickness=4)
        st.configure("TNotebook", background=PANEL)
        st.configure("TNotebook.Tab", padding=(14, 6))
        st.configure("TLabelframe", background=PANEL)
        st.configure("TLabelframe.Label", background=PANEL, foreground=MUTED)
        st.configure("TCheckbutton", background=PANEL)
        st.configure("TRadiobutton", background=PANEL)
        st.configure("Section.TLabel", font=self.f_bold, foreground=ACCENT_DARK)
        st.configure("Treeview", rowheight=26)

    # ================================================================================================== layout
    def _build(self):
        r = self.root
        # --- menu row (drawn in the window, as in the design)
        mb = tk.Frame(r, bg=PANEL, height=36, highlightthickness=0)
        mb.pack(fill="x")
        tk.Label(mb, text="soapff", bg=PANEL, fg="#23313b", font=self.f_bold).pack(side="left", padx=(18, 14))
        self.menus = {}
        for name, items in self._menu_spec():
            b = tk.Menubutton(mb, text=name, bg=PANEL, fg=TEXT, activebackground="#eef2f4", relief="flat",
                              font=self.f_small, padx=10, pady=6)
            m = tk.Menu(b, tearoff=0)
            for it in items:
                if it is None:
                    m.add_separator()
                else:
                    m.add_command(label=it[0], command=it[1], accelerator=it[2] if len(it) > 2 else "")
            b.configure(menu=m)
            b.pack(side="left")
            self.menus[name] = m
        tk.Label(mb, text=f"soapff-GUI  v{__version__}", bg=PANEL, fg=MUTED, font=self.f_small).pack(side="right", padx=18)
        tk.Frame(r, bg=LINE, height=1).pack(fill="x")

        # --- toolbar
        tb = tk.Frame(r, bg=BG, height=56)
        tb.pack(fill="x")
        ttk.Button(tb, text="New calculation", style="Ghost.TButton", command=self.new_calc).pack(side="left", padx=(18, 4), pady=10)
        ttk.Button(tb, text="Supercell", style="Ghost.TButton", command=self.new_supercell).pack(side="left", padx=(0, 4), pady=10)
        ttk.Button(tb, text="Surface", style="Ghost.TButton", command=self.new_surface).pack(side="left", padx=(0, 4), pady=10)
        ttk.Button(tb, text="Add water", style="Ghost.TButton", command=self.new_water).pack(side="left", padx=(0, 14), pady=10)
        tk.Frame(tb, bg=LINE, width=1).pack(side="left", fill="y", pady=12)
        self.crumb_settings = tk.Label(tb, text="Calculation settings", bg=BG, fg=MUTED, font=self.f_small, cursor="hand2")
        self.crumb_settings.pack(side="left", padx=(16, 8))
        self.crumb_settings.bind("<Button-1>", lambda e: self.open_settings())
        tk.Label(tb, text="›", bg=BG, fg=FAINT, font=self.f_small).pack(side="left")
        self.crumb_run = tk.Label(tb, text="Run & analysis", bg=BG, fg=MUTED, font=self.f_small, cursor="hand2")
        self.crumb_run.pack(side="left", padx=8)
        self.crumb_run.bind("<Button-1>", lambda e: self.show_results())
        self.b_start_top = ttk.Button(tb, text="Start calculation", style="Accent.TButton", command=self.start_or_stop)
        self.b_start_top.pack(side="right", padx=18, pady=10)
        tk.Frame(r, bg=LINE, height=1).pack(fill="x")

        # --- status bar (packed before the body so it stays at the bottom)
        sb = tk.Frame(r, bg=PANEL, height=30)
        sb.pack(side="bottom", fill="x")
        tk.Frame(r, bg=LINE, height=1).pack(side="bottom", fill="x")
        self.l_status = tk.Label(sb, text="", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_status.pack(side="left", padx=18, pady=5)
        self.l_status_r = tk.Label(sb, text="LAMMPS · SOAP", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_status_r.pack(side="right", padx=18)

        body = tk.Frame(r, bg=BG)
        body.pack(fill="both", expand=True)
        self._sidebar(body)
        main = tk.Frame(body, bg=BG)
        main.pack(side="left", fill="both", expand=True, padx=24, pady=(12, 10))
        self._main(main)

    def _menu_spec(self):
        return [
            ("File", [("New calculation…", self.new_calc, "⌘N"), ("New calculation from a supercell (OVITO Replicate)…", self.new_supercell),
                          ("New calculation from a surface (slab in vacuum)…", self.new_surface),
                          ("New calculation with water added (packmol)…", self.new_water),
                          ("Change structure…", self.new_calc),
                          ("Open existing output folder…", self.open_project), None,
                          ("Fetch a folder from another computer…", self.fetch_dialog), None,
                          ("Open output folder", self.open_out), None, ("Quit", self.on_close)]),
            ("Edit", [("Add line (free input)", self.new_line), ("Delete everything after minimize", self.cut_after_minimize), None,
                      ("Add fix…", lambda: self.open_catalog("fix")), ("Add compute…", lambda: self.open_catalog("compute")),
                      ("Add other command…", lambda: self.open_catalog("other")), None, ("Delete line", self.delete_line), None,
                      ("Move up", lambda: self.move_line(-1)), ("Move down", lambda: self.move_line(1)), None,
                      ("Undo", self.do_undo, "⌘Z")]),
            ("View", [("in.lammps", lambda: self.show_tab("in.lammps")), ("system.settings (force field)",
                                                                         lambda: self.show_tab("system.settings")),
                      ("system.data (atoms, bonds, angles)", lambda: self.show_tab("system.data")), None,
                      ("Analysis results", self.show_results), ("Plot (choose columns)…", self.open_plotter), ("Open in OVITO", self.open_ovito), None,
                      ("Calculation history", self.show_history), ("Structure library", self.show_library)]),
            ("Run", [("Start calculation", self.start, "⌘R"), ("Stop", self.stop), None,
                      ("Rebuild input (from the MD conditions)", self.regenerate), ("Analyze", self.analyze_now), None,
                      ("Calculation settings…", self.open_settings)]),
            ("Help", [("Open README", self.open_readme), ("About soapff-GUI", self.about)]),
        ]

    def _sidebar(self, body):
        sd = tk.Frame(body, bg=SIDE, width=228)
        sd.pack(side="left", fill="y")
        sd.pack_propagate(False)
        tk.Frame(body, bg=LINE, width=1).pack(side="left", fill="y")
        tk.Label(sd, text="Project", bg=SIDE, fg=MUTED, font=self.f_head).pack(anchor="w", padx=20, pady=(22, 8))
        self.nav = {}
        for key, label, cmd in (("new", "New calculation", lambda: None), ("history", "Calculation history", self.show_history),
                                ("library", "Structure library", self.show_library)):
            row = tk.Frame(sd, bg=SIDE, cursor="hand2")
            row.pack(fill="x", padx=10, pady=2)
            bar = tk.Frame(row, bg=SIDE, width=3)
            bar.pack(side="left", fill="y")
            lab = tk.Label(row, text=f"  ▤  {label}", bg=SIDE, fg=TEXT, font=self.f_small, anchor="w", padx=6, pady=8)
            lab.pack(side="left", fill="x", expand=True)
            for w in (row, lab):
                w.bind("<Button-1>", lambda e, c=cmd: c())
            self.nav[key] = (row, bar, lab)
        self._nav_select("new")
        tk.Frame(sd, bg=LINE, height=1).pack(fill="x", padx=18, pady=14)
        tk.Label(sd, text="Recent calculations", bg=SIDE, fg=MUTED, font=self.f_head).pack(anchor="w", padx=20, pady=(0, 6))
        self.recent = tk.Frame(sd, bg=SIDE)
        self.recent.pack(fill="x", padx=14)

    def _nav_select(self, key):
        if not self.root.winfo_exists() or not self.nav["new"][0].winfo_exists():
            return
        for k, (row, bar, lab) in self.nav.items():
            on = k == key
            for w in (row, lab):
                w.configure(bg=ACCENT_BG if on else SIDE)
            bar.configure(bg=ACCENT if on else SIDE)
            lab.configure(font=self.f_bold if on else self.f_small)

    def _card(self, parent, **pack):
        outer = tk.Frame(parent, bg=PANEL, highlightthickness=1, highlightbackground=LINE)
        outer.pack(**pack)
        return outer

    def _main(self, main):
        tk.Label(main, text="Input file", bg=BG, fg="#23313b", font=self.f_title).pack(anchor="w")
        tk.Label(main, text="Select a command line to edit it.", bg=BG, fg=MUTED, font=self.f_small).pack(anchor="w", pady=(0, 8))

        # structure card
        c = self._card(main, fill="x")
        left = tk.Frame(c, bg=PANEL)
        left.pack(side="left", padx=16, pady=6, fill="x", expand=True)
        tk.Label(left, text="Structure", bg=PANEL, fg=FAINT, font=pick_font(["Hiragino Sans"], 10)).pack(anchor="w")
        self.l_struct = tk.Label(left, text="None selected (choose a structure file with New calculation)", bg=PANEL, fg=TEXT, font=self.f_bold)
        self.l_struct.pack(anchor="w")
        self.l_struct2 = tk.Label(left, text="", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_struct2.pack(anchor="w")
        tk.Frame(c, bg=LINE, width=1).pack(side="left", fill="y", pady=8)
        mid = tk.Frame(c, bg=PANEL)
        mid.pack(side="left", padx=16, pady=6, fill="x", expand=True)
        tk.Label(mid, text="Output files", bg=PANEL, fg=FAINT, font=pick_font(["Hiragino Sans"], 10)).pack(anchor="w")
        self.l_outfile = tk.Label(mid, text="in.lammps", bg=PANEL, fg=TEXT, font=self.f_mono)
        self.l_outfile.pack(anchor="w")
        self.l_outdir = tk.Label(mid, text="", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_outdir.pack(anchor="w")
        ttk.Button(c, text="Change structure", style="Link.TButton", command=self.new_calc).pack(side="right", padx=16)

        # bottom row: selected command | run
        # editor above, command form and run below, with a draggable divider (drag it to make the editor bigger)
        self.pw = tk.PanedWindow(main, orient="vertical", bg=BG, sashwidth=10, sashrelief="flat", bd=0, opaqueresize=True)
        self.pw.pack(fill="both", expand=True, pady=(10, 0))
        bottom = tk.Frame(self.pw, bg=BG, height=250)
        self._bottom = bottom
        rc = tk.Frame(bottom, bg=PANEL, highlightthickness=1, highlightbackground=LINE, width=440, height=250)
        rc.pack(side="right", fill="y", padx=(14, 0))
        rc.pack_propagate(False)
        # the command form has a fixed height (its fields scroll): a form with many fields used to grow and squeeze the editor
        # (the in.lammps tabs and the lines) down to nothing
        self.fc = tk.Frame(bottom, bg=PANEL, highlightthickness=1, highlightbackground=LINE, height=250)
        self.fc.pack(side="left", fill="both", expand=True)
        self.fc.pack_propagate(False)
        self._form_area()
        self._run_area(rc)

        # editor card
        ec = tk.Frame(self.pw, bg=PANEL, highlightthickness=1, highlightbackground=LINE)
        self.pw.add(ec, stretch="always", minsize=90)
        self.pw.add(bottom, stretch="never", minsize=250, height=252)
        head = tk.Frame(ec, bg=PANEL)
        head.pack(fill="x", padx=14, pady=(10, 0))
        self.tabs = {}
        for f in FILES:
            t = tk.Label(head, text=f" {f} ", bg=PANEL, fg="#23313b", font=pick_font(["Menlo"], 11, "bold"),
                         highlightthickness=1, highlightbackground="#c9d1d7", padx=8, pady=4, cursor="hand2")
            t.pack(side="left", padx=(0, 6))
            t.bind("<Button-1>", lambda e, f=f: self.show_tab(f))
            self.tabs[f] = t
        self.l_save = tk.Label(head, text="LAMMPS input · autosave", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_save.pack(side="right")
        ttk.Button(head, text="Open as text (large, paste)", style="Tiny.TButton", command=self.open_text_editor).pack(
            side="right", padx=(0, 12))
        self.b_expand = ttk.Button(head, text="⤢ Expand", style="Tiny.TButton", command=self.toggle_expand)
        self.b_expand.pack(side="right", padx=(0, 6))
        # boundary of the box (system.settings); changing it rebuilds the force field (k-space solver, walls, centring)
        bf = tk.Frame(head, bg=PANEL)
        bf.pack(side="right", padx=(0, 18))
        tk.Label(bf, text="Boundary", bg=PANEL, fg=MUTED, font=self.f_small).pack(side="left", padx=(0, 6))
        self.v_bnd = []
        for ax in "xyz":
            tk.Label(bf, text=ax, bg=PANEL, fg=FAINT, font=self.f_small).pack(side="left", padx=(4, 2))
            vb = tk.StringVar(value="p")
            cb = ttk.Combobox(bf, textvariable=vb, values=AXIS_STYLES, width=3, state="readonly", font=self.f_mono_s)
            cb.pack(side="left")
            cb.bind("<<ComboboxSelected>>", lambda e: self.change_boundary([v.get() for v in self.v_bnd]))
            self.v_bnd.append(vb)
        tk.Frame(ec, bg=LINE, height=1).pack(fill="x", pady=(8, 0))
        tf = tk.Frame(ec, bg=PANEL)
        tf.pack(fill="both", expand=True)
        self.ed = tk.Text(tf, bg=PANEL, fg=TEXT, font=self.f_mono, relief="flat", highlightthickness=0, wrap="none",
                          cursor="arrow", spacing1=3, spacing3=3, padx=0, insertwidth=0, height=6,
                          tabs=("560p",), selectbackground=PANEL, inactiveselectbackground=PANEL)
        ys = ttk.Scrollbar(tf, orient="vertical", command=self.ed.yview)
        self.ed.configure(yscrollcommand=ys.set)
        ys.pack(side="right", fill="y")
        self.ed.pack(side="left", fill="both", expand=True)
        self.ed.tag_configure("lnum", foreground=FAINT)
        self.ed.tag_configure("cmd", foreground=CMD)
        self.ed.tag_configure("arg", foreground=TEXT)
        self.ed.tag_configure("comment", foreground=MUTED)
        self.ed.tag_configure("stripe", background="#fafbfc")
        self.ed.tag_configure("selline", background=SEL_BG)
        self.ed.tag_configure("bar", foreground=ACCENT)
        self.ed.tag_configure("hint", foreground=ACCENT, font=(self.f_mono[0], self.f_mono[1] - 1))
        self.ed.tag_raise("selline")
        self.ed.bind("<Button-1>", self.on_click)
        self.ed.bind("<Double-Button-1>", lambda e: (self.on_click(e), self.focus_form()))
        for key in ("<Up>", "<Down>"):
            self.ed.bind(key, lambda e, d=(-1 if key == "<Up>" else 1): self.move_sel(d))
        self.ed.bind("<Key>", self.on_key)
        for ev in ("<<Paste>>", "<Command-v>", "<Control-v>"):
            self.ed.bind(ev, self.paste_lines)
        for ev in ("<<Copy>>", "<Command-c>", "<Control-c>"):
            self.ed.bind(ev, self.copy_line)
        foot = tk.Frame(ec, bg=PANEL)
        foot.pack(fill="x", padx=14, pady=6)
        for text, kind in (("+ fix", "fix"), ("+ compute", "compute"), ("+ Other (run, timestep, ...)", "other")):
            ttk.Button(foot, text=text, style="Link.TButton", command=lambda k=kind: self.open_catalog(k)).pack(side="left", padx=(0, 6))
        ttk.Button(foot, text="+ Add line (free input)", style="Link.TButton", command=self.new_line).pack(side="left", padx=(0, 6))
        tk.Label(foot, text="Select a line to edit · Up/Down: move · Enter: add line · Delete: remove", bg=PANEL, fg=MUTED,
                 font=self.f_small).pack(side="right")
        self.b_insert_anchor = foot


    def _form_area(self):
        head = tk.Frame(self.fc, bg=PANEL)
        head.pack(fill="x", padx=16, pady=(10, 0))
        self.l_formtitle = tk.Label(head, text="Selected command", bg=PANEL, fg="#23313b", font=self.f_bold)
        self.l_formtitle.pack(side="left")
        self.l_formline = tk.Label(head, text="", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_formline.pack(side="right")
        tk.Frame(self.fc, bg=LINE, height=1).pack(fill="x", pady=(8, 0))
        self.l_preview = tk.Label(self.fc, text="", bg=CODE_BG, fg=TEXT, font=self.f_mono_s, anchor="w", padx=12, pady=8)
        self.l_preview.pack(side="bottom", fill="x", padx=16, pady=(0, 12))
        tk.Label(self.fc, text="Generated line", bg=PANEL, fg=MUTED, font=self.f_small).pack(side="bottom", anchor="w", padx=16,
                                                                                         pady=(6, 2))
        area = tk.Frame(self.fc, bg=PANEL)
        area.pack(fill="both", expand=True, padx=(16, 4), pady=(8, 0))
        cv = tk.Canvas(area, bg=PANEL, highlightthickness=0, height=40)
        sb = ttk.Scrollbar(area, orient="vertical", command=cv.yview)
        cv.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        cv.pack(side="left", fill="both", expand=True)
        self.fields = tk.Frame(cv, bg=PANEL)
        win = cv.create_window(0, 0, window=self.fields, anchor="nw")
        self.fields.bind("<Configure>", lambda e: cv.configure(scrollregion=cv.bbox("all")))
        cv.bind("<Configure>", lambda e: cv.itemconfigure(win, width=e.width))
        cv.bind("<Enter>", lambda e: cv.bind_all("<MouseWheel>", lambda ev: cv.yview_scroll(-1 if ev.delta > 0 else 1, "units")))
        cv.bind("<Leave>", lambda e: cv.unbind_all("<MouseWheel>"))

    def _run_area(self, rc):
        head = tk.Frame(rc, bg=PANEL)
        head.pack(fill="x", padx=16, pady=(8, 0))
        tk.Label(head, text="Run", bg=PANEL, fg="#23313b", font=self.f_bold).pack(side="left", pady=4)
        self.b_start = ttk.Button(head, text="Start calculation", style="Small.Accent.TButton", command=self.start_or_stop)
        self.b_start.pack(side="right")
        wrow = tk.Frame(rc, bg=PANEL)                  # Run on: chosen computer, + to register one, - to remove it
        wrow.pack(fill="x", padx=16, pady=(6, 0))
        tk.Label(wrow, text="Run on", bg=PANEL, fg=MUTED, font=self.f_small).pack(side="left")
        self.v_where = tk.StringVar()
        self.cb_where = ttk.Combobox(wrow, textvariable=self.v_where, state="readonly", width=12, font=self.f_small)
        self.cb_where.pack(side="left", padx=(6, 4))
        self.cb_where.bind("<<ComboboxSelected>>", lambda e: self.choose_where())
        ttk.Button(wrow, text="+ Add", width=5, style="Tiny.TButton", command=lambda: self.quick_add_computer(self.root)).pack(
            side="left", padx=(0, 4))
        ttk.Button(wrow, text="- Remove", width=7, style="Tiny.TButton", command=self.remove_computer).pack(side="left")
        tk.Frame(rc, bg=LINE, height=1).pack(fill="x", pady=(8, 0))
        row = tk.Frame(rc, bg=PANEL)
        row.pack(fill="x", padx=16, pady=(12, 0))
        self.dot = tk.Canvas(row, width=10, height=10, bg=PANEL, highlightthickness=0)
        self.dot.pack(side="left")
        self.dot_id = self.dot.create_oval(1, 1, 9, 9, fill=FAINT, outline="")
        self.l_phase = tk.Label(row, text="Idle", bg=PANEL, fg=TEXT, font=self.f_small)
        self.l_phase.pack(side="left", padx=6)
        self.l_pct = tk.Label(row, text="0%", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_pct.pack(side="right")
        self.pb = ttk.Progressbar(rc, maximum=1.0, style="Accent.Horizontal.TProgressbar")
        self.pb.pack(fill="x", padx=16, pady=(8, 0))
        self.l_thermo = tk.Label(rc, text="", bg=PANEL, fg=MUTED, font=self.f_mono_s, anchor="w")
        self.l_thermo.pack(fill="x", padx=16, pady=(6, 0))
        acts = tk.Frame(rc, bg=PANEL)                  # packed first (at the bottom) so that "Analyze" is always visible
        acts.pack(side="bottom", fill="x", padx=16, pady=(0, 10))
        for text, cmd in (("Analyze", self.analyze_now), ("Plot", self.open_plotter), ("OVITO", self.open_ovito),
                          ("Folder", self.open_out)):
            ttk.Button(acts, text=text, style="Link.TButton", command=cmd).pack(side="left", padx=(0, 6))
        tk.Label(acts, text="↑ analyzed automatically after a run", bg=PANEL, fg=FAINT, font=self.f_small).pack(side="right")
        lf = tk.Frame(rc, bg=CODE_BG, highlightthickness=1, highlightbackground="#e1e6ea")
        lf.pack(fill="both", expand=True, padx=16, pady=(8, 8))
        self.log = tk.Text(lf, bg=CODE_BG, fg=TEXT, font=pick_font(["Menlo"], 10), relief="flat", height=4, wrap="none",
                           highlightthickness=0, padx=8, pady=6)
        self.log.pack(fill="both", expand=True)
        self.log.insert("end", LOG_PLACEHOLDER, "muted")
        self.log.tag_config("muted", foreground=FAINT)
        self.log.tag_config("warn", foreground="#a15c00")
        self.log.tag_config("err", foreground=ERR)
        self.log.config(state="disabled")

    # ================================================================================================== editor
    def data_section(self, i):
        """the section (Masses, Atoms, Bonds ...) the line i of system.data is in ('' for the header)."""
        for k in range(i, -1, -1):
            w = self.lines[k].split()
            if w and w[0] in DATA_COLS and (len(w) == 1 or w[1].startswith("#")):
                return w[0]
        return ""

    def render_data(self):
        ed = self.ed
        y = ed.yview()[0]
        ed.configure(state="normal")
        ed.delete("1.0", "end")
        text = "\n".join(f" {i + 1:>6}    {l}" for i, l in enumerate(self.lines)) + "\n"
        ed.insert("end", text, "arg")
        for i, l in enumerate(self.lines):
            w = l.split()
            if w and w[0] in DATA_COLS:
                ed.tag_add("cmd", f"{i + 1}.0", f"{i + 1}.end")
            elif "#" in l:
                ed.tag_add("comment", f"{i + 1}.{l.index('#') + 11}", f"{i + 1}.end")
        ed.configure(state="disabled")
        ed.yview_moveto(y)
        self.mark_selection()

    def mark_selection(self):
        ed = self.ed
        ed.tag_remove("selline", "1.0", "end")
        if self.sel is not None and self.lines:
            ed.tag_add("selline", f"{self.sel + 1}.0", f"{self.sel + 1}.end+1c")
            ed.see(f"{self.sel + 1}.0")

    def render(self):
        if self.tab in READONLY:
            return self.render_data()
        ed = self.ed
        y = ed.yview()[0]
        ed.configure(state="normal")
        ed.delete("1.0", "end")
        if not self.lines:
            ed.insert("end", "\n   Choose a structure file and soapff-cui writes the LAMMPS input (in.lammps, system.settings) and shows it here.\n", "comment")
            ed.insert("end", "   Start with \"New calculation\" or \"Change structure\" above.\n", "comment")
        for i, line in enumerate(self.lines):
            sel = i == self.sel
            ed.insert("end", "▎" if sel else " ", "bar")
            ed.insert("end", f"{i + 1:>4}    ", "lnum")
            c, a, cm = L.split(line)
            if c:
                ed.insert("end", c, "cmd")
                pad = " " * max(1, 13 - len(c)) if len(c) < 13 else " "
                ed.insert("end", pad + " ".join(a), "arg")
                if cm:
                    ed.insert("end", "  " + cm, "comment")
            else:
                ed.insert("end", line.strip(), "comment")
            if sel:
                ed.insert("end", "\tSelected · click to edit", "hint")
            ed.insert("end", "\n")
            tags = ("selline",) if sel else (("stripe",) if i % 2 else ())
            for t in tags:
                ed.tag_add(t, f"{i + 1}.0", f"{i + 1}.end+1c")
        ed.configure(state="disabled")
        ed.yview_moveto(y)
        if self.sel is not None and self.lines:
            ed.see(f"{self.sel + 1}.0")

    def on_click(self, e):
        if not self.lines:
            return "break"
        i = int(self.ed.index(f"@{e.x},{e.y}").split(".")[0]) - 1
        if 0 <= i < len(self.lines):
            self.select(i)
        self.ed.focus_set()
        return "break"

    def on_key(self, e):
        if e.keysym in ("Return", "KP_Enter") and self.tab not in READONLY:
            self.new_line()
            return "break"
        if e.keysym in ("Delete", "BackSpace") and self.tab not in READONLY:
            self.delete_line()
            return "break"
        if e.keysym in ("Up", "Down"):
            return None
        return "break"

    def toggle_expand(self):
        """hide / show the command form and run area under the editor (the editor takes the whole height)."""
        if not getattr(self, "_expanded", False):
            self.pw.forget(self._bottom)
            self._expanded = True
            self.b_expand.config(text="⤡ Restore")
        else:
            self.pw.add(self._bottom, stretch="never", minsize=130, height=240)
            self._expanded = False
            self.b_expand.config(text="⤢ Expand")

    def copy_line(self, e=None):
        """Copy: the selected line to the clipboard."""
        if self.sel is not None and self.lines:
            self.root.clipboard_clear()
            self.root.clipboard_append(self.lines[self.sel])
        return "break"

    def paste_lines(self, e=None):
        """Paste: the clipboard text (one or many lines) is inserted below the selected line."""
        if self.tab in READONLY:
            return "break"
        try:
            text = self.root.clipboard_get()
        except tk.TclError:
            return "break"
        new = [l.rstrip("\r") for l in text.replace("\r\n", "\n").split("\n")]
        while new and not new[-1].strip():
            new.pop()
        if new:
            self.insert_lines(new)
            self.log_add(f"Pasted {len(new)} line(s)\n", None)
        return "break"

    def open_text_editor(self):
        """the whole file of the editor tab as plain text in a large window: edit, copy and paste freely, then Apply."""
        name = self.tab
        if not self.docs.get(name):
            messagebox.showinfo("soapff-GUI", "First choose a structure file and create the input")
            return
        ro = name in READONLY
        w = tk.Toplevel(self.root)
        w.title(f"{name}" + (" (read-only)" if ro else " (text editing)"))
        sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        w.geometry(f"{min(1100, sw - 100)}x{min(780, sh - 120)}")
        f = ttk.Frame(w, padding=8)
        f.pack(fill="both", expand=True)
        bar = ttk.Frame(f)
        bar.pack(side="bottom", fill="x", pady=(8, 0))
        tf = ttk.Frame(f)
        tf.pack(fill="both", expand=True)
        t = tk.Text(tf, font=self.f_mono, wrap="none", undo=True, bg=CODE_BG, fg=TEXT, relief="flat", padx=8, pady=6)
        ys = ttk.Scrollbar(tf, orient="vertical", command=t.yview)
        xs = ttk.Scrollbar(tf, orient="horizontal", command=t.xview)
        t.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        ys.pack(side="right", fill="y")
        xs.pack(side="bottom", fill="x")
        t.pack(side="left", fill="both", expand=True)
        t.insert("1.0", "\n".join(self.docs[name]) + "\n")
        if ro:
            t.config(state="disabled")
        msg = tk.StringVar(value="" if ro else "Edit, copy and paste here, then press \"Apply\" to replace the lines of the input with this text (undo with ⌘Z)")
        ttk.Label(bar, textvariable=msg, foreground=MUTED).pack(side="left")

        def apply():
            lines = t.get("1.0", "end-1c").replace("\r\n", "\n").split("\n")
            while lines and not lines[-1].strip():
                lines.pop()
            self.push_undo()
            self.docs[name] = lines
            self.sels[name] = min(self.sel or 0, max(len(lines) - 1, 0)) if lines else None
            self.sel = self.sels[name]
            self.dirty.add(name)
            self.render()
            self.build_form()
            self.schedule_save()
            msg.set(f"Applied ({len(lines)} lines)")
        if not ro:
            ttk.Button(bar, text="Apply", style="Accent.TButton", command=apply).pack(side="right")
            w.bind("<Command-s>", lambda ev: apply())
            w.bind("<Control-s>", lambda ev: apply())
        ttk.Button(bar, text="Close", style="Ghost.TButton", command=w.destroy).pack(side="right", padx=8)
        t.focus_set()

    def move_sel(self, d):
        if not self.lines:
            return "break"
        i = 0 if self.sel is None else min(max(self.sel + d, 0), len(self.lines) - 1)
        self.select(i)
        return "break"

    def select(self, i):
        self.sel = i
        if self.tab in READONLY:
            self.mark_selection()          # thousands of lines: only move the highlight
        else:
            self.render()
        self.build_form()

    def focus_form(self):
        for w in self.fields.winfo_children():
            for c in w.winfo_children():
                if isinstance(c, ttk.Entry):
                    c.focus_set()
                    return

    def build_form(self):
        for w in self.fields.winfo_children():
            w.destroy()
        self.form_vars = []
        if self.sel is None or not self.lines:
            self.l_formtitle.config(text="Selected command")
            self.l_formline.config(text="")
            self.l_preview.config(text="")
            return
        line = self.lines[self.sel]
        if self.tab in READONLY:
            return self.data_form(line)
        c, a, cm = L.split(line)
        self.l_formline.config(text=f"Line {self.sel + 1}")
        if not c:
            self.l_formtitle.config(text="Selected line: comment / blank")
            specs = [("Content", line.strip(), None)]
        else:
            self.l_formtitle.config(text=f"Selected command: {c}")
            labs = L.labels(c, a)
            specs = [(lab, a[k] if k < len(a) else "", k) for k, lab in enumerate(labs)]
            if cm:
                specs.append(("Comment", cm, "comment"))
        self._form_cmd = c
        if c == "boundary":
            self.boundary_form(a)
            return
        if c == "include" and a and a[0] in FILES:
            tk.Button(self.fields, text=f"Open {a[0]}", command=lambda f=a[0]: self.show_tab(f), relief="flat", bg=ACCENT_BG,
                      fg=ACCENT, activebackground=SEL_BG, font=self.f_small, padx=10, pady=4, cursor="hand2").pack(
                anchor="w", pady=(0, 6))
        col = 0
        rowf = None
        self.fc.update_idletasks()
        per_row = max(2, min(6, (self.fc.winfo_width() - 48) // 140)) if self.fc.winfo_width() > 50 else 6
        for lab, val, k in specs:
            if lab is None:                       # fixed keyword (e.g. 'temp' of fix nvt)
                self.form_vars.append((k, None, val))
                continue
            if col % per_row == 0:
                rowf = tk.Frame(self.fields, bg=PANEL)
                rowf.pack(fill="x", anchor="w")
            cell = tk.Frame(rowf, bg=PANEL)
            cell.pack(side="left", padx=(0, 12), pady=(0, 6))
            tk.Label(cell, text=lab, bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w")
            v = tk.StringVar(value=val)
            ttk.Entry(cell, textvariable=v, width=14 if k != "comment" else 28, style="Field.TEntry",
                      font=self.f_mono_s).pack(anchor="w")
            v.trace_add("write", lambda *x: self.form_changed())
            if k is None:                         # free text (comment / blank line): Enter shows the command's fields
                cell.winfo_children()[-1].bind("<Return>", lambda e: (self.build_form(), self.focus_form()))
            self.form_vars.append((k, v, val))
            col += 1
        self.update_preview()

    def data_form(self, line):
        """a line of system.data, its columns named (read only); the type number is shown with its ClayFF name."""
        self._form_cmd = ""
        self.form_vars = []
        sec = self.data_section(self.sel)
        self.l_formline.config(text=f"Line {self.sel + 1} (read-only)")
        self.l_formtitle.config(text=f"system.data: {sec or 'Header'}")
        body, _, comment = line.partition("#")
        vals = body.split()
        cols = DATA_COLS.get(sec, []) if vals and vals[0] not in DATA_COLS else []
        names = {v: k for k, v in self.type_ids().items()} if sec == "Atoms" else {}
        per_row = 6
        rowf = None
        for k, v in enumerate(vals):
            if k % per_row == 0:
                rowf = tk.Frame(self.fields, bg=PANEL)
                rowf.pack(fill="x", anchor="w")
            cell = tk.Frame(rowf, bg=PANEL)
            cell.pack(side="left", padx=(0, 12), pady=(0, 6))
            lab = cols[k] if k < len(cols) else ("Value" if cols else f"{k + 1}")
            tk.Label(cell, text=lab, bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w")
            shown = v
            if sec == "Atoms" and k == 2 and v.isdigit() and int(v) in names:
                shown = f"{v} ({names[int(v)]})"
            tk.Label(cell, text=shown, bg=CODE_BG, fg=TEXT, font=self.f_mono_s, width=14, anchor="w", padx=6, pady=4).pack(anchor="w")
        if comment.strip():
            tk.Label(self.fields, text=f"Comment: {comment.strip()}", bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w")
        tk.Label(self.fields, text="system.data is the table of atoms, bonds and angles written by soapff-cui (view only here). "
                                   "To change it, change the structure file or the boundary and rebuild.", wraplength=560, justify="left",
                 bg=PANEL, fg=FAINT, font=self.f_small).pack(anchor="w", pady=(4, 0))
        self.l_preview.config(text=line.strip())

    def boundary_form(self, a):
        """the boundary line: x, y, z chosen from the LAMMPS styles; applying rebuilds the force field."""
        b = (list(a) + ["p", "p", "p"])[:3]
        row = tk.Frame(self.fields, bg=PANEL)
        row.pack(fill="x", anchor="w")
        vs = []
        for ax, x in zip("xyz", b):
            cell = tk.Frame(row, bg=PANEL)
            cell.pack(side="left", padx=(0, 12), pady=(0, 6))
            tk.Label(cell, text=ax, bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w")
            v = tk.StringVar(value=x)
            ttk.Combobox(cell, textvariable=v, values=AXIS_STYLES, width=5, state="readonly", font=self.f_mono_s).pack(anchor="w")
            vs.append(v)
        ttk.Button(row, text="Rebuild with this boundary", style="Small.Accent.TButton",
                   command=lambda: self.change_boundary([v.get() for v in vs])).pack(side="left", padx=(8, 0), pady=(14, 0))
        tk.Label(self.fields, text="p: periodic   f: fixed   s: shrinks to fit the atoms   m: expands to fit the atoms (per face, e.g. fm). "
                                   "Changing it rebuilds the force field (long-range Coulomb method etc.) and the input.",
                 bg=PANEL, fg=MUTED, font=self.f_small, wraplength=560, justify="left").pack(anchor="w")
        self.form_vars = []
        self.l_preview.config(text=f"boundary {' '.join(b)}")
        for v in vs:
            v.trace_add("write", lambda *x: self.l_preview.config(text="boundary " + " ".join(v.get() for v in vs)))

    def form_line(self):
        c = self._form_cmd
        if not c:
            return self.form_vars[0][1].get() if self.form_vars else ""
        args, comment = [], ""
        for k, v, val in self.form_vars:
            if k == "comment":
                comment = v.get().strip()
                if comment and not comment.startswith("#"):
                    comment = "# " + comment
            else:
                args.append(val if v is None else v.get().strip())
        return L.join(c, args, comment)

    def update_preview(self):
        if self._form_cmd == "boundary":
            return
        self.l_preview.config(text=self.form_line() if self.form_vars else "")

    def form_changed(self):
        if self.sel is None or not self.form_vars:
            return
        new = self.form_line()
        if new != self.lines[self.sel]:
            self.push_undo()
            self.lines[self.sel] = new
            self.render()
            self.schedule_save()
        self.update_preview()

    def push_undo(self):
        if not self.undo or self.undo[-1] != self.lines:
            self.undo.append(list(self.lines))
            self.undo = self.undo[-100:]

    def do_undo(self):
        if not self.undo:
            return
        self.lines = self.undo.pop()
        self.sel = min(self.sel or 0, len(self.lines) - 1) if self.lines else None
        self.render()
        self.build_form()
        self.schedule_save()

    def open_catalog(self, kind):
        """LAMMPS fix / compute styles (and run, timestep ...) to insert after the selected line of in.lammps."""
        if not self.docs["in.lammps"]:
            messagebox.showinfo("soapff-GUI", "First choose a structure file and create the input")
            return
        if self.tab != "in.lammps":
            self.show_tab("in.lammps")
        items = catalog.entries(kind) + (catalog.entries("special") if kind in ("other", "compute") else [])
        w = tk.Toplevel(self.root)
        w.title({"fix": "Add fix", "compute": "Add compute", "other": "Add command"}[kind])
        w.geometry("980x600")
        w.configure(bg=PANEL)
        left = tk.Frame(w, bg=PANEL)
        left.pack(side="left", fill="both", padx=(12, 6), pady=12)
        q = tk.StringVar()
        ttk.Entry(left, textvariable=q, style="Field.TEntry", width=40).pack(fill="x")
        tk.Label(left, text="Filter by name or description", bg=PANEL, fg=FAINT, font=self.f_small).pack(anchor="w", pady=(2, 6))
        tv = ttk.Treeview(left, show="tree", height=22)
        tv.column("#0", width=420)
        tv.pack(fill="both", expand=True)
        right = tk.Frame(w, bg=PANEL)
        right.pack(side="left", fill="both", expand=True, padx=(6, 12), pady=12)
        l_title = tk.Label(right, text="", bg=PANEL, fg="#23313b", font=self.f_bold, anchor="w", justify="left", wraplength=480)
        l_title.pack(fill="x")
        l_style = tk.Label(right, text="", bg=PANEL, fg=CMD, font=self.f_mono, anchor="w")
        l_style.pack(fill="x", pady=(2, 8))
        l_desc = tk.Label(right, text="", bg=PANEL, fg=TEXT, font=self.f_small, anchor="w", justify="left", wraplength=480)
        l_desc.pack(fill="x")
        tk.Label(right, text="Lines to be inserted", bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w", pady=(12, 2))
        l_tpl = tk.Label(right, text="", bg=CODE_BG, fg=TEXT, font=self.f_mono_s, anchor="w", justify="left", padx=10, pady=8,
                         wraplength=480)
        l_tpl.pack(fill="x")
        l_fields = tk.Label(right, text="", bg=PANEL, fg=MUTED, font=self.f_small, anchor="w", justify="left", wraplength=480)
        l_fields.pack(fill="x", pady=(8, 0))
        bar = tk.Frame(right, bg=PANEL)
        bar.pack(side="bottom", fill="x")
        where = tk.Label(right, text="", bg=PANEL, fg=FAINT, font=self.f_small, anchor="w")
        where.pack(side="bottom", fill="x", pady=(0, 8))
        cur = {}

        def fill(*_):
            tv.delete(*tv.get_children())
            k = q.get().strip().lower()
            cats = {}
            for n, e in enumerate(items):
                if k and k not in (e["style"] + e["title"] + e["desc"]).lower():
                    continue
                if e["category"] not in cats:
                    cats[e["category"]] = tv.insert("", "end", text=e["category"], open=True)
                tv.insert(cats[e["category"]], "end", iid=f"e{n}", text=f"{e['style']:<16}  {e['title']}")
            first = next((i for c in tv.get_children() for i in tv.get_children(c)), None)
            if first:
                tv.selection_set(first)
                tv.see(first)

        def chosen():
            sel = tv.selection()
            return items[int(sel[0][1:])] if sel and sel[0].startswith("e") else None

        def show(*_):
            e = chosen()
            cur.clear()
            if not e:
                return
            cur.update(e)
            lines = self.catalog_lines(e)
            l_title.config(text=e["title"])
            l_style.config(text=(f"{e['kind']} {e['style']}" if e["kind"] in ("fix", "compute") else e["style"]))
            l_desc.config(text=e["desc"].replace("{N}", "N").replace("{T}", "T") or " ")
            l_tpl.config(text="\n".join(lines))
            labs = [x for x in (catalog.BASE if e["kind"] in ("fix", "compute") else []) + e["labels"] if x]
            l_fields.config(text=("Fields: " + ", ".join(labs)) if labs else "")
            n = self.sel + 1 if self.sel is not None else len(self.lines)
            where.config(text=f"Inserted after line {n} (you can edit the values in the fields below afterwards)")

        def insert(*_):
            e = chosen()
            if e:
                self.insert_lines(self.catalog_lines(e))
                w.destroy()

        def manual():
            if cur.get("doc"):
                webbrowser.open(cur["doc"])
        ttk.Button(bar, text="Insert", style="Accent.TButton", command=insert).pack(side="right")
        ttk.Button(bar, text="LAMMPS manual", style="Ghost.TButton", command=manual).pack(side="right", padx=8)
        ttk.Button(bar, text="Close", style="Ghost.TButton", command=w.destroy).pack(side="left")
        q.trace_add("write", fill)
        tv.bind("<<TreeviewSelect>>", show)
        tv.bind("<Double-1>", insert)
        w.bind("<Return>", insert)
        fill()

    def cut_after_minimize(self):
        """keep in.lammps up to the last minimize (and the walls of fixed faces); the MD is added from the lists."""
        lines = self.docs["in.lammps"]
        extra = L.after_minimize(lines)
        if not extra:
            messagebox.showinfo("soapff-GUI", "There are no commands after minimize")
            return
        if not messagebox.askyesno("soapff-GUI", f"This deletes the {len(extra)} line(s) of in.lammps after minimize"
                                                 f" ({extra[0]} …). Continue? (Undo restores them)"):
            return
        if self.tab != "in.lammps":
            self.show_tab("in.lammps")
        self.push_undo()
        self.lines = L.basic_part(lines)
        self.sel = len(self.lines) - 1
        self.render()
        self.build_form()
        self.schedule_save()

    def warn_assignment(self, rep):
        """after the force field is made: a clear warning when atoms were not found in the library (fallback) or their
        SOAP type differs from the bond-graph ClayFF type (the library may lack the mineral)."""
        if not isinstance(rep, dict):
            return
        fb = int(rep.get("fallback", 0))
        tc = rep.get("type_check") or {}
        nd = int(tc.get("differ", 0))
        if not fb and not nd:
            return
        lines = [f"Library: {rep.get('library', '?')}"]
        if fb:
            per = rep.get("per_type") or {}
            kinds = ", ".join(f"{t} {v['fallback']}" for t, v in per.items() if v.get("fallback"))
            lines.append(f"- Atoms not determined by the library, typed from the bond connectivity instead: {fb} ({kinds})")
        if nd:
            kinds = ", ".join(f"{k}: {v}" for k, v in (tc.get("kinds") or {}).items())
            lines.append(f"- Atoms whose SOAP type differs from the bond-connectivity type: {nd} ({kinds}; 'bond type → SOAP type')")
        lines.append("The environment of this mineral may be missing from the library. In Calculation settings › Library and details, choose the bundled library "
                     f"({runner.default_library()}) or add reference structures of this mineral to the library. "
                     "The per-atom values are in charges.csv (columns source, bond_type).")
        self.log_add("Warning: " + " ".join(lines) + "\n", "warn")
        messagebox.showwarning("soapff-GUI", "\n".join(lines))

    def type_ids(self):
        """{ClayFF type name: type number} from the Masses of system.data."""
        out = {}
        try:
            L_ = open(os.path.join(self.out, "system.data")).read().splitlines()
            i = next(k for k, l in enumerate(L_) if l.strip() == "Masses") + 2
            while i < len(L_) and L_[i].strip():
                p = L_[i].split()
                if "#" in L_[i]:
                    out[L_[i].split("#")[1].strip()] = int(p[0])
                i += 1
        except (OSError, StopIteration):
            pass
        return out

    def catalog_lines(self, e):
        """the line(s) a catalog entry inserts, with the placeholders filled in."""
        lines = self.docs["in.lammps"]
        dt = L.timestep(lines, default=float(getattr(runner.import_soapff(), "TIMESTEP_FS", 0.5)))
        s = self.s
        if e["kind"] == "special":
            core = runner.import_soapff()
            tid = self.type_ids()
            every = max(1, int(round(100.0 / dt)))
            msd = core.analysis_lines(tid, every)
            if e["style"] == "soapff:diffusion":
                return msd or ["# There is no water or ion, so there is no MSD for the diffusion coefficient"]
            if e["style"] == "soapff:swelling":
                return self.swelling_block(dt)
            names = [t for t, _ in sorted(tid.items(), key=lambda kv: kv[1])]
            md = core.md_protocol(s["temp"], s["nvt_ps"], s["npt_ps"], s["press"], s["baro"], names, s["dump_ps"], msd,
                                  s["relax_ps"], dt).splitlines()
            T = s["temp"]
            if "# production" in md:                  # the relaxation only; the production is chosen by the user
                md = md[:md.index("# production")]
            return [f"velocity all create {T} 4928459 mom yes rot no dist gaussian", f"timestep {dt:g}"] + md
        fixes, comps = L.ids(lines, "fix"), L.ids(lines, "compute")
        want = e["style"].replace("/", "_") if e["kind"] in ("fix", "compute") else e["style"]
        kind = e["kind"] if e["kind"] in ("fix", "compute") else ("dump" if e["style"] == "dump" else
                                                                 ("group" if e["style"] == "group" else None))
        if kind == "dump":
            want = "d"
        if kind == "group":
            want = "g"
        uid = L.unique_id(lines, kind, want) if kind in ("fix", "compute", "dump") else want
        steps = int(round(float(s.get("nvt_ps", 10.0)) * 1000.0 / dt)) or 20000
        line = e["template"].format(ID=uid, T=f"{float(s['temp']):g}", P=f"{float(s['press']):g}", DT=f"{dt:g}", N=steps,
                                    LASTFIX=(fixes[-1] if fixes else "ID"), LASTCOMPUTE=(comps[-1] if comps else "ID"))
        return [line]

    def ask_swelling(self):
        """parameters of the swelling curve (dialog) -> dict, or None when cancelled."""
        p = dict(DEFAULTS["swelling"], **(self.s.get("swelling") or {}))
        w = tk.Toplevel(self.root)
        w.title("Swelling curve (method of Cygan et al. 2004)")
        w.transient(self.root)
        f = ttk.Frame(w, padding=16)
        f.pack(fill="both", expand=True)
        ttk.Label(f, text="Runs NPT while adding water molecules to the interlayer and records d001 (= lz) for each state.", wraplength=440).grid(
            row=0, column=0, columnspan=3, sticky="w", pady=(0, 10))
        v = {}
        for r, (k, lab, note) in enumerate((("one_to", "Add 1 at a time (up to, molecules)", "Cygan et al.: 23"),
                                             ("two_to", "Then add 2 at a time (up to, molecules)", "Cygan et al.: 99"),
                                             ("eq_ps", "Equilibration NPT per state (ps)", "Cygan et al.: 30"),
                                             ("avg_ps", "Averaging NPT afterwards (ps)", "Cygan et al.: 20")), start=1):
            ttk.Label(f, text=lab).grid(row=r, column=0, sticky="w", pady=2)
            v[k] = tk.StringVar(value=f"{p[k]:g}" if isinstance(p[k], float) else str(p[k]))
            ttk.Entry(f, textvariable=v[k], width=8).grid(row=r, column=1, sticky="w", padx=6)
            ttk.Label(f, text=note, foreground=MUTED).grid(row=r, column=2, sticky="w")
        v["down"] = tk.BooleanVar(value=bool(p["down"]))
        ttk.Checkbutton(f, text="Also remove 2 at a time at the end (dehydration)", variable=v["down"]).grid(row=5, column=0, columnspan=3,
                                                                                         sticky="w", pady=(6, 0))
        ttk.Label(f, text="Or specify the water counts directly").grid(row=8, column=0, sticky="w", pady=(10, 0))
        v["waters"] = tk.StringVar(value=" ".join(str(x) for x in p.get("waters", [])))
        ttk.Entry(f, textvariable=v["waters"], width=34).grid(row=8, column=1, columnspan=2, sticky="w", padx=6, pady=(10, 0))
        ttk.Label(f, text="Example: 4 8 12 16 20 25 32 40 48 59 (if empty, the settings above are used; no dehydration)", foreground=MUTED).grid(
            row=9, column=0, columnspan=3, sticky="w")
        ttk.Label(f, text=f"Temperature {self.s['temp']:g} K, pressure {self.s['press']:g} atm (the MD conditions of the Calculation settings). "
                          "If the structure has no water types, system.data and system.settings are rewritten (manual changes to system.settings are lost).",
                  foreground=MUTED, wraplength=440).grid(row=6, column=0, columnspan=3, sticky="w", pady=(8, 0))
        out = {}

        def ok():
            try:
                out.update(one_to=int(v["one_to"].get()), two_to=int(v["two_to"].get()), eq_ps=float(v["eq_ps"].get()),
                           avg_ps=float(v["avg_ps"].get()), down=bool(v["down"].get()),
                           waters=sorted({int(x) for x in v["waters"].get().replace(",", " ").split()}))
                assert 1 <= out["one_to"] <= out["two_to"] and out["eq_ps"] > 0 and out["avg_ps"] > 0
                assert all(x >= 1 for x in out["waters"])
            except (ValueError, AssertionError):
                out.clear()
                messagebox.showerror("soapff-GUI", "Check the numbers (1 ≤ upper limit for adding 1 at a time ≤ upper limit for adding 2 at a time; times > 0)", parent=w)
                return
            w.destroy()
        bar = ttk.Frame(f)
        bar.grid(row=7, column=0, columnspan=3, sticky="ew", pady=(12, 0))
        ttk.Button(bar, text="Add to in.lammps", style="Accent.TButton", command=ok).pack(side="right")
        ttk.Button(bar, text="Cancel", style="Ghost.TButton", command=w.destroy).pack(side="left")
        w.grab_set()
        self.root.wait_window(w)
        if out:
            self.s["swelling"] = dict(out)
            save_settings(self.s)
        return out or None

    def swelling_block(self, dt):
        """the LAMMPS block of the swelling curve (soapff_cui.core.swelling_lines) for this output folder; writes water.mol
        and, when the structure has no water types, rewrites system.data / system.settings with them reserved."""
        core = runner.import_soapff()
        p = self.ask_swelling()
        if not p:
            return []
        if "ow" not in self.type_ids():
            core.reserve_water_types(self.out, log=lambda m: self.log_add(m + "\n", None))
            for f in ("system.data", "system.settings"):
                self.docs[f] = open(os.path.join(self.out, f)).read().splitlines()
                self.undos[f] = []
        tid = self.type_ids()
        D = open(os.path.join(self.out, "system.data")).read().splitlines()
        nat = int(D[2].split()[0])
        k = D.index(next(l for l in D if l.startswith("Atoms")))
        atoms = [l.split() for l in D[k + 2:k + 2 + nat]]
        n_water0 = sum(1 for a in atoms if int(a[2]) == tid["ow"])
        n_nonwater = nat - 3 * n_water0
        mol_max = max(int(a[1]) for a in atoms)
        tilted = any(l.strip().endswith("xy xz yz") and any(abs(float(x)) > 1e-8 for x in l.split()[:3]) for l in D[:20])
        S = open(os.path.join(self.out, "system.settings")).read().splitlines()
        bond = int(next(l for l in S if l.startswith("bond_coeff") and "water" in l).split()[1])
        angle = int(next(l for l in S if l.startswith("angle_coeff") and "HOHw" in l).split()[1])
        try:
            water = json.load(open(os.path.join(self.out, "report.json"))).get("water_model", "spc")
        except (OSError, ValueError):
            water = "spc"
        open(os.path.join(self.out, "water.mol"), "w").write(core.water_molecule(tid, bond, angle, water=water))
        shake = (bond, angle) if water in core.RIGID_WATER else None
        if p.get("waters"):
            up, down = list(p["waters"]), []                    # the listed water counts, in increasing order; no deswelling
        else:
            up, down = core.swelling_schedule(p["one_to"], p["two_to"], p["down"])
        self.log_add(f"Swelling curve: adding water in {len(up)} states ({up[0]} to {up[-1]} molecules)" + (f", removing it in {len(down)} states" if down else "")
                     + f" (NPT {p['eq_ps']:g} + {p['avg_ps']:g} ps each). Wrote water.mol\n", None)
        return core.swelling_lines(tid, n_nonwater, n_water0, mol_max, self.s["temp"], self.s["press"], dt, tilted, up, down,
                                   eq_ps=p["eq_ps"], avg_ps=p["avg_ps"], shake=shake)

    def insert_lines(self, new):
        if self.tab in READONLY:
            return
        self.push_undo()
        i = (self.sel + 1) if self.sel is not None else len(self.lines)
        self.lines[i:i] = list(new)
        self.select(i)
        self.schedule_save()

    def new_line(self):
        """an empty line below the selected one, typed in the "Content" field of the form (e.g. 'run 10000')."""
        if self.tab in READONLY:
            return
        if not self.lines and not self.out:
            messagebox.showinfo("soapff-GUI", "First choose a structure file and create the input")
            return
        self.insert_lines([""])
        self.focus_form()

    def delete_line(self):
        if self.tab in READONLY:
            return
        if self.sel is None or not self.lines:
            return
        self.push_undo()
        del self.lines[self.sel]
        self.sel = min(self.sel, len(self.lines) - 1) if self.lines else None
        self.render()
        self.build_form()
        self.schedule_save()

    def move_line(self, d):
        if self.tab in READONLY:
            return
        if self.sel is None:
            return
        j = self.sel + d
        if 0 <= j < len(self.lines):
            self.push_undo()
            self.lines[self.sel], self.lines[j] = self.lines[j], self.lines[self.sel]
            self.select(j)
            self.schedule_save()

    def schedule_save(self):
        if self.tab in READONLY:
            return
        self.dirty.add(self.tab)
        if self._save_after:
            self.root.after_cancel(self._save_after)
        self._save_after = self.root.after(400, self.save_input)
        self.l_save.config(text=f"{FILE_NOTE[self.tab]} · save pending")

    def save_input(self):
        self._save_after = None
        if self.out and os.path.isdir(self.out):
            for f in sorted(self.dirty):
                if self.docs[f]:
                    with open(os.path.join(self.out, f), "w") as fh:
                        fh.write("\n".join(self.docs[f]) + "\n")
            self.dirty.clear()
            self.l_save.config(text=f"{FILE_NOTE[self.tab]} · autosaved ({time.strftime('%H:%M:%S')})")

    def show_tab(self, f):
        """show in.lammps or system.settings in the editor."""
        if f not in FILES:
            return
        self.sels[self.tab] = self.sel
        self.tab = f
        self.sel = self.sels[f]
        if self.sel is None and self.lines:
            self.sel = 0
        for name, t in self.tabs.items():
            on = name == f
            t.configure(fg="#23313b" if on else MUTED, highlightbackground="#c9d1d7" if on else PANEL,
                        bg=PANEL if on else "#f1f4f6")
        self.l_save.config(text=f"{FILE_NOTE[f]} · autosaved")
        self.l_outfile.config(text=f)
        self.render()
        self.build_form()

    def current_boundary(self):
        for line in self.docs["system.settings"]:
            c, a, _ = L.split(line)
            if c == "boundary" and len(a) >= 3:
                return a[:3]
        b = (self.report or {}).get("boundary") or self.s.get("boundary", "p p p")
        return (b.split() + ["p", "p", "p"])[:3]

    def change_boundary(self, new):
        """a new boundary rebuilds the force field and the input (soapff-cui chooses PPPM / PPPM+slab / MSM, centres the
        empty gap of a non-periodic direction and puts walls on fixed faces)."""
        old = self.current_boundary()
        new = [x.strip() for x in new]
        if new == old or not self.out:
            self.set_boundary_widgets(old)
            return
        if self.busy():
            messagebox.showinfo("soapff-GUI", "Cannot be changed during a run")
            self.set_boundary_widgets(old)
            return
        if new != ["p", "p", "p"] and float(self.s.get("npt_ps", 0)) > 0:
            messagebox.showerror("soapff-GUI", "NPT cannot be used when a direction is non-periodic (the pressure of a box containing vacuum is undefined).\n"
                                               "Set NPT to 0 under \"MD conditions\" in the Calculation settings, then choose again")
            self.set_boundary_widgets(old)
            return
        if not self.structure or not os.path.isfile(self.structure):
            messagebox.showerror("soapff-GUI", "The structure file is unknown, so the input cannot be rebuilt. Choose one with \"New calculation\"")
            self.set_boundary_widgets(old)
            return
        if not messagebox.askyesno("soapff-GUI", f"Change the boundary from {' '.join(old)} to {' '.join(new)} and rebuild the force field and input.\n"
                                                 "(The long-range Coulomb method and the handling of non-periodic directions also change.)\n"
                                                 "Lines of in.lammps and system.settings that you edited by hand are lost. Continue?"):
            self.set_boundary_widgets(old)
            return
        self.s["boundary"] = " ".join(new)
        save_settings(self.s)
        self.make_input(self.structure, self.out)

    def set_boundary_widgets(self, b):
        for v, x in zip(self.v_bnd, b):
            v.set(x)

    # ================================================================================================== projects
    def new_calc(self):
        if self.busy():
            messagebox.showinfo("soapff-GUI", "A run is in progress. Choose after it finishes or is stopped")
            return
        p = filedialog.askopenfilename(title="Choose a structure file", initialdir=self.s["last_dir"], filetypes=STRUCTURE_TYPES)
        if not p:
            return
        self.s["last_dir"] = os.path.dirname(p)
        stem = os.path.splitext(os.path.basename(p))[0]
        out = os.path.join(os.path.dirname(p), f"{stem}_soapff")
        if os.path.isfile(os.path.join(out, "in.lammps")):
            ans = messagebox.askyesnocancel("soapff-GUI", f"{out} contains a previous input.\n\nYes: open the previous input as it is\n"
                                                          f"No: rebuild the force field and input (overwrite)")
            if ans is None:
                return
            if ans:
                self.load_project(p, out)
                return
        self.make_input(p, out)

    def new_supercell(self):
        """a structure file -> supercell nx x ny x nz by OVITO's Replicate (ovito Python module; ASE without it), shown with
        its size, optionally opened in OVITO, then the force field and input are made from it (as New calculation)."""
        if self.busy():
            messagebox.showinfo("soapff-GUI", "A run is in progress. Choose after it finishes or is stopped")
            return
        p = filedialog.askopenfilename(title="Choose the original structure file", initialdir=self.s["last_dir"], filetypes=STRUCTURE_TYPES)
        if not p:
            return
        self.s["last_dir"] = os.path.dirname(p)
        try:
            at = runner.import_soapff().read_structure(p)
        except Exception as e:
            messagebox.showerror("soapff-GUI", f"Cannot read the structure: {e}")
            return
        L0 = at.cell.lengths()
        w = tk.Toplevel(self.root)
        w.title("Supercell (OVITO Replicate)")
        w.transient(self.root)
        f = ttk.Frame(w, padding=16)
        f.pack(fill="both", expand=True)
        ttk.Label(f, text=f"{os.path.basename(p)}: {len(at)} atoms, a b c = {L0[0]:.2f} {L0[1]:.2f} {L0[2]:.2f} Å",
                  font=self.f_bold).grid(row=0, column=0, columnspan=6, sticky="w")
        ttk.Label(f, text="Number of repeats in each direction (same as OVITO's Replicate; the box is enlarged too)", foreground=MUTED).grid(
            row=1, column=0, columnspan=6, sticky="w", pady=(4, 8))
        v = []
        for k, ax in enumerate("abc"):
            ttk.Label(f, text=f"{ax}").grid(row=2, column=2 * k, sticky="e")
            sv = tk.StringVar(value=str(max(1, math.ceil(20.0 / L0[k])) if L0[k] < 20.0 else 1))
            ttk.Spinbox(f, from_=1, to=50, width=4, textvariable=sv).grid(row=2, column=2 * k + 1, sticky="w", padx=(4, 14))
            v.append(sv)
        info = tk.StringVar()
        ttk.Label(f, textvariable=info).grid(row=3, column=0, columnspan=6, sticky="w", pady=(8, 0))
        tool = "OVITO (ovito module)" if runner.has_ovito_module() else "ASE (no ovito module, so the same tiling is done with ASE)"
        ttk.Label(f, text=f"Using: {tool}", foreground=MUTED).grid(row=4, column=0, columnspan=6, sticky="w")

        def nums():
            try:
                n = [int(x.get()) for x in v]
                return n if min(n) >= 1 else None
            except ValueError:
                return None

        def update(*_):
            n = nums()
            info.set(f"→ {len(at) * n[0] * n[1] * n[2]} atoms, a b c = {L0[0] * n[0]:.2f} {L0[1] * n[1]:.2f} {L0[2] * n[2]:.2f} Å"
                     + (" (an edge is shorter than 10 Å, twice the SOAP cutoff)" if min(L0[k] * n[k] for k in range(3)) < 10.0 else "")
                     if n else "Enter integers of 1 or more")
        for x in v:
            x.trace_add("write", update)
        update()

        def make_file():
            n = nums()
            if not n:
                return None
            try:
                out, nat, abc, used = runner.replicate(p, n)
            except Exception as e:
                messagebox.showerror("soapff-GUI", f"Could not create the supercell: {e}", parent=w)
                return None
            self.log_add(f"Supercell {n[0]}×{n[1]}×{n[2]} ({used}): {out} ({nat} atoms)\n", None)
            return out

        def preview():
            out = make_file()
            if out:
                try:
                    runner.open_in_ovito(out)
                except Exception as e:
                    messagebox.showerror("soapff-GUI", str(e), parent=w)

        def go():
            out = make_file()
            if not out:
                return
            w.destroy()
            stem = os.path.splitext(os.path.basename(out))[0]
            self.make_input(out, os.path.join(os.path.dirname(out), f"{stem}_soapff"))
        bar = ttk.Frame(f)
        bar.grid(row=5, column=0, columnspan=6, sticky="ew", pady=(14, 0))
        ttk.Button(bar, text="Build force field from this structure", style="Accent.TButton", command=go).pack(side="right")
        ttk.Button(bar, text="Check in OVITO", style="Ghost.TButton", command=preview).pack(side="right", padx=8)
        ttk.Button(bar, text="Cancel", style="Ghost.TButton", command=w.destroy).pack(side="left")

    def _derived_dialog(self, title, header):
        """the common frame of the Surface / Add water dialogs: (window, body frame, button bar)"""
        w = tk.Toplevel(self.root)
        w.title(title)
        w.transient(self.root)
        f = ttk.Frame(w, padding=16)
        f.pack(fill="both", expand=True)
        ttk.Label(f, text=header, font=self.f_bold).pack(anchor="w")
        body = ttk.Frame(f)
        body.pack(fill="both", expand=True, pady=(8, 0))
        bar = ttk.Frame(f)
        bar.pack(fill="x", pady=(14, 0))
        return w, body, bar

    def _build_from(self, path, w):
        w.destroy()
        stem = os.path.splitext(os.path.basename(path))[0]
        self.make_input(path, os.path.join(os.path.dirname(path), f"{stem}_soapff"))

    def new_surface(self):
        """the current calculation (its structure and the ClayFF types of charges.csv) -> a slab in vacuum of a face (h k l)
        (soapff_cui.surfaces: every Si and tetrahedral Al kept whole, broken O healed by bond valence, OH2 on a single
        octahedral metal removed), written next to the structure; then the force field and input are made from it."""
        if self.busy():
            messagebox.showinfo("soapff-GUI", "A run is in progress. Choose after it finishes or is stopped")
            return
        csvf = os.path.join(self.out, "charges.csv") if self.out else ""
        if not (self.structure and os.path.isfile(self.structure) and os.path.isfile(csvf)):
            messagebox.showinfo("soapff-GUI", "Build the force field of the bulk crystal first (New calculation): the surface is cut "
                                              "from the current calculation, using its ClayFF types")
            return
        core = runner.import_soapff()
        import csv as _csv
        from soapff_cui import surfaces
        at = core.read_structure(self.structure)
        types = [r["clayff_type"] for r in _csv.DictReader(open(csvf))]
        if len(types) != len(at):
            messagebox.showerror("soapff-GUI", f"charges.csv has {len(types)} atoms, the structure {len(at)}")
            return
        keep = [t not in ("ow", "hw") for t in types]
        at = at[keep]; types = [t for t, k in zip(types, keep) if k]
        w, f, bar = self._derived_dialog("Surface (slab in vacuum)", f"{os.path.basename(self.structure)}: {len(at)} atoms"
                                         + (" (its water is left out)" if not all(keep) else ""))
        ttk.Label(f, text="Face (h k l) of the unit cell. Only a face that cuts no cation-O bond can be cut: the basal plane between the "
                          "layers of a layered mineral, (0 0 1) for the bundled minerals. Edges (bonds cut, O that would need H) are not supported", foreground=MUTED, wraplength=560, justify="left").grid(row=0, column=0, columnspan=6, sticky="w")
        hkl = []
        for k, (ax, v0) in enumerate(zip("hkl", (0, 0, 1))):
            ttk.Label(f, text=ax).grid(row=1, column=2 * k, sticky="e", pady=8)
            sv = tk.StringVar(value=str(v0))
            ttk.Spinbox(f, from_=-1, to=1, width=4, textvariable=sv).grid(row=1, column=2 * k + 1, sticky="w", padx=(4, 14))
            hkl.append(sv)
        SLAB_THICKNESS, SLAB_VACUUM = 16.0, 30.0           # fixed: at least 16 A of mineral, 30 A of vacuum (room for water)
        ttk.Label(f, text="Repeat the cell (a b c) before the cut").grid(row=2, column=0, columnspan=3, sticky="w", pady=2)
        rep = []
        for k in range(3):
            sv = tk.StringVar(value="1")
            ttk.Spinbox(f, from_=1, to=8, width=4, textvariable=sv).grid(row=2, column=3 + k, sticky="w", padx=(4, 4))
            rep.append(sv)
        keepv = tk.BooleanVar(value=True)
        ttk.Checkbutton(f, text="Cut inside the cell: the cell is not changed, the slab fills it except for one unit cell along the "
                                "direction out of the plane (the empty part takes the water)", variable=keepv).grid(row=3, column=0, columnspan=6, sticky="w")
        ttk.Label(f, text="Unit cells kept (blank: all but one)").grid(row=4, column=0, columnspan=3, sticky="w", pady=2)
        perv = tk.StringVar(value="")
        ttk.Entry(f, textvariable=perv, width=6, style="Field.TEntry").grid(row=4, column=3, sticky="w")
        ttk.Label(f, text="Extra length of the cell (Å)").grid(row=5, column=0, columnspan=3, sticky="w", pady=2)
        extv = tk.StringVar(value="0")
        ttk.Entry(f, textvariable=extv, width=6, style="Field.TEntry").grid(row=5, column=3, sticky="w")
        info = tk.StringVar(value="")
        ttk.Label(f, textvariable=info, wraplength=560, justify="left").grid(row=6, column=0, columnspan=6, sticky="w", pady=(10, 0))
        made = {}

        def make_file():
            try:
                m = tuple(int(x.get()) for x in hkl)
                th, vac = SLAB_THICKNESS, SLAB_VACUUM
                rp = tuple(max(1, int(x.get())) for x in rep)
                per = (int(perv.get()) if perv.get().strip() else "auto") if keepv.get() else None
                ext = float(extv.get())
            except ValueError:
                info.set("Enter integers for h k l, the repeats and the unit cells, and a number for the extra length")
                return None
            kc = True
            if per is not None:
                vac = ext                                   # inside the cell: the vacuum is only the extra length of the cell
            if m == (0, 0, 0):
                info.set("(0 0 0) is not a face")
                return None
            key = (m, th, vac, rp, kc, per)
            if key in made:
                return made[key]
            w.config(cursor="watch"); w.update()
            try:
                slab, typ, inf = surfaces.make_slab(at, types, m, min_thickness=th, vacuum=vac, keep_cell=kc, repeat=rp, periods=per)
            except Exception as e:
                messagebox.showerror("soapff-GUI", f"Could not cut the surface: {e}", parent=w)
                return None
            finally:
                w.config(cursor="")
            stem = os.path.splitext(os.path.basename(self.structure))[0]
            path = os.path.join(os.path.dirname(self.structure), f"{stem}_{surfaces.hkl_name(m)}"
                                + ("" if rp == (1, 1, 1) else "_x" + "".join(map(str, rp))) + ".extxyz")
            import ase.io
            ase.io.write(path, slab, format="extxyz")
            txt = (f"→ {os.path.basename(path)}: {inf['atoms']} atoms, {inf['cut_bonds']} cation-O bonds cut"
                   + (" (basal cleavage: no bond cut, the interlayer ions shared by the two faces)" if inf["basal"] else
                      f", {inf['h_added']} H added, {inf['oh2_removed']} OH2 removed") + f"; net formal charge {inf['formal_charge']:+d}")
            if inf.get("cell_doubled"):
                txt += (f"\nThe cell was lengthened to {inf['cell_length_periods']} times its length along the direction out of the plane, "
                        "so that the empty part is at least 6 A wide (the two faces must not see each other)")
            if inf["formal_charge"]:
                txt += "\nThis face stays charged at this protonation: the charges are not shifted; the net charge is reported and LAMMPS adds a uniform neutralising background (a warning)"
            info.set(txt)
            self.log_add(f"Surface {m}: {path} ({inf['atoms']} atoms, formal charge {inf['formal_charge']:+d})\n", None)
            made[key] = path
            return path

        def preview():
            p = make_file()
            if p:
                try:
                    runner.open_in_ovito(p)
                except Exception as e:
                    messagebox.showerror("soapff-GUI", str(e), parent=w)

        def go():
            p = make_file()
            if p:
                self._build_from(p, w)
        ttk.Button(bar, text="Build force field from this structure", style="Accent.TButton", command=go).pack(side="right")
        ttk.Button(bar, text="Check in OVITO", style="Ghost.TButton", command=preview).pack(side="right", padx=8)
        ttk.Button(bar, text="Cancel", style="Ghost.TButton", command=w.destroy).pack(side="left")

    def new_water(self):
        """a structure (the current one, or one chosen) -> the empty slices along z (vacuum, interlayer) listed; the chosen
        ones are filled with water by packmol at a density or a number of molecules (soapff_cui.water); the result is written
        next to the structure and the force field and input are made from it (the water is typed by rule, with the water
        model of the Calculation settings)."""
        if self.busy():
            messagebox.showinfo("soapff-GUI", "A run is in progress. Choose after it finishes or is stopped")
            return
        from soapff_cui import water as W
        try:
            exe = W.find_packmol()
        except FileNotFoundError as e:
            messagebox.showerror("soapff-GUI", str(e))
            return
        p = self.structure if self.structure and os.path.isfile(self.structure) else filedialog.askopenfilename(
            title="Choose the structure to add water to", initialdir=self.s["last_dir"], filetypes=STRUCTURE_TYPES)
        if not p:
            return
        try:
            at = runner.import_soapff().read_structure(p)
        except Exception as e:
            messagebox.showerror("soapff-GUI", f"Cannot read the structure: {e}")
            return
        dens = tk.StringVar(value="1.0")
        gl = W.gaps(at, density=1.0)
        w, f, bar = self._derived_dialog("Add water (packmol)", f"{os.path.basename(p)}: {len(at)} atoms; empty regions along z (normal to the a-b plane)")
        if not gl:
            ttk.Label(f, text="No empty region thicker than 2.5 Å: make room first (Surface: vacuum; a larger c)").grid(row=0, column=0, sticky="w")
            ttk.Button(bar, text="Close", style="Ghost.TButton", command=w.destroy).pack(side="left")
            return
        top = ttk.Frame(f)
        top.grid(row=0, column=0, columnspan=6, sticky="w")
        ttk.Label(top, text="Density (g/cm³)").pack(side="left")
        ttk.Entry(top, textvariable=dens, width=6, style="Field.TEntry").pack(side="left", padx=6)
        ttk.Label(top, text="sets the numbers below (area × thickness between the atom planes); a number can be edited",
                  foreground=MUTED).pack(side="left", padx=6)
        for k, h in enumerate(("Fill", "z from (Å)", "to (Å)", "Thickness (Å)", "Area (Å²)", "Molecules")):
            ttk.Label(f, text=h, font=self.f_bold).grid(row=1, column=k, sticky="w", padx=6, pady=(10, 2))
        rows = []
        for i, g in enumerate(gl):
            use = tk.BooleanVar(value=g["thickness"] >= 6.0)
            nv = tk.StringVar(value=str(g["n_water"]))
            ttk.Checkbutton(f, variable=use).grid(row=2 + i, column=0, padx=6)
            for k, x in enumerate((g["lo"], g["hi"], g["thickness"], g["area"]), start=1):
                ttk.Label(f, text=f"{x:.2f}" if k < 4 else f"{x:.1f}").grid(row=2 + i, column=k, sticky="w", padx=6)
            ttk.Entry(f, textvariable=nv, width=7, style="Field.TEntry").grid(row=2 + i, column=5, sticky="w", padx=6)
            rows.append((g, use, nv))

        def redo(*_):
            try:
                d = float(dens.get())
            except ValueError:
                return
            for g, use, nv in rows:
                nv.set(str(W.n_water(g["area"] * g["thickness"], d)))
        dens.trace_add("write", redo)
        status = tk.StringVar(value=f"packmol: {exe}; water at least 2.0 Å from every atom")
        ttk.Label(f, textvariable=status, foreground=MUTED, wraplength=620, justify="left").grid(row=3 + len(gl), column=0, columnspan=6, sticky="w", pady=(10, 0))
        made = {}

        def make_file(then):
            try:
                chosen = [(g, int(nv.get())) for g, use, nv in rows if use.get() and int(nv.get()) > 0]
            except ValueError:
                status.set("The numbers of molecules must be integers")
                return
            if not chosen:
                status.set("Choose at least one region with molecules")
                return
            key = tuple((g["lo"], n) for g, n in chosen)
            if key in made:
                then(made[key])
                return
            status.set(f"packmol: placing {sum(n for _, n in chosen)} water molecules…")
            w.config(cursor="watch")

            def work():
                try:
                    new, rep = W.add_water(at, chosen, log=lambda m: self.q.put(("log", m + "\n", None)))
                    stem = os.path.splitext(os.path.basename(p))[0]
                    out = os.path.join(os.path.dirname(p), f"{stem}_water.extxyz")
                    import ase.io
                    ase.io.write(out, new, format="extxyz")
                    res = (out, rep, None)
                except Exception as e:
                    res = (None, None, e)

                def done():
                    if not w.winfo_exists():
                        return
                    w.config(cursor="")
                    out, rep, err = res
                    if err:
                        status.set(f"packmol failed: {err}")
                        return
                    bad = [r for r in rep if not r["converged"]]
                    status.set(f"→ {os.path.basename(out)}: " + ", ".join(f"{r['n']} molecules in z {r['lo']:.1f}-{r['hi']:.1f} Å" for r in rep)
                               + (" (packmol did not converge in a region: fewer molecules or a lower density)" if bad else ""))
                    self.log_add(f"Water added: {out}\n", None)
                    made[key] = out
                    then(out)
                self.q.put(("call", done))                       # back on the Tk thread (the queue is polled there)
            threading.Thread(target=work, daemon=True).start()

        def preview():
            def show(out):
                try:
                    runner.open_in_ovito(out)
                except Exception as e:
                    messagebox.showerror("soapff-GUI", str(e), parent=w)
            make_file(show)
        ttk.Button(bar, text="Build force field from this structure", style="Accent.TButton",
                   command=lambda: make_file(lambda out: self._build_from(out, w))).pack(side="right")
        ttk.Button(bar, text="Check in OVITO", style="Ghost.TButton", command=preview).pack(side="right", padx=8)
        ttk.Button(bar, text="Cancel", style="Ghost.TButton", command=w.destroy).pack(side="left")

    def open_project(self):
        d = filedialog.askdirectory(title="Choose the output folder (the one containing in.lammps)", initialdir=self.s["last_dir"])
        if d and os.path.isfile(os.path.join(d, "in.lammps")):
            st = next((h["structure"] for h in self.s["history"] if h.get("out") == d), "")
            self.load_project(st, d)
        elif d:
            messagebox.showerror("soapff-GUI", "This folder has no in.lammps")

    def make_input(self, structure, out):
        """soapff make in a worker thread (writes system.data, system.settings, in.lammps ...)."""
        md = {k: self.s[k] for k in ("temp", "relax_ps", "nvt_ps", "npt_ps", "press", "baro", "boundary", "threshold",
                                     "cutoff", "dump_ps", "charges", "water")}
        if not os.path.isfile(self.s["library"]):
            messagebox.showerror("soapff-GUI", f"SOAP library not found: {self.s['library']}\nChoose one in the Calculation settings")
            return
        self.log_clear()
        self.set_state("Preparing", FAINT, 0)
        self.l_phase.config(text="Building force field…")

        def go():
            put = self.q.put
            try:
                put(("log", f"soapff make {structure}\n  -> {out}\n", None))
                rep = runner.prepare(structure, out, self.s["library"],
                                     log=lambda m: put(("log", m + "\n", "warn" if m.startswith("WARNING") else None)), md=False, **md)
                put(("made", structure, out, rep))
            except Exception as e:
                put(("log", f"ERROR {type(e).__name__}: {e}\n", "err"))
                put(("make_failed", str(e)))
        threading.Thread(target=go, daemon=True).start()

    def load_project(self, structure, out):
        self.structure, self.out = structure, out
        try:
            docs = {f: open(os.path.join(out, f)).read().splitlines() for f in FILES}
        except OSError as e:
            messagebox.showerror("soapff-GUI", f"Cannot read the input: {e}")
            return
        try:
            self.report = json.load(open(os.path.join(out, "report.json")))
        except (OSError, ValueError):
            self.report = None
        self.docs = docs
        self.undos = {f: [] for f in FILES}
        self.dirty = set()
        self.sels = {f: next((i for i, l in enumerate(docs[f]) if L.split(l)[0] in ("fix", "boundary")), 0) for f in FILES}
        self.tab = "system.settings"                          # show_tab stores this selection, then switches
        self.sel = self.sels["system.settings"]
        self.show_tab("in.lammps")
        self.set_boundary_widgets(self.current_boundary())
        self.show_project()
        self.remember(status=status_of(self.history_entry()) or STATUS_INPUT)

    def show_project(self):
        name = os.path.basename(self.structure) if self.structure else os.path.basename(self.out)
        r = self.report
        atoms = r["atoms"] if r else None
        self.l_struct.config(text=f"{name}  ·  {atoms:,} atoms" if atoms else name)
        if r:
            ms = min((p["min_similarity"] for p in r["per_type"].values()), default=float("nan"))
            roles = r.get("roles", {})
            self.l_struct2.config(text=f"Clay {roles.get('mineral', 0)} · Water {roles.get('water', 0)} · Ions {roles.get('ion', 0)}"
                                       f"   SOAP {r['from_soap']} / fallback {r['fallback']} · min similarity {ms:.3f}",
                                  fg="#a15c00" if r["fallback"] else MUTED)
        else:
            self.l_struct2.config(text="")
        self.l_outdir.config(text=self.out)
        self.root.title(f"soapff-GUI — {os.path.splitext(name)[0]}")
        self.update_status()

    def history_entry(self):
        return next((h for h in self.s["history"] if h.get("out") == self.out), {})

    def remember(self, status):
        if not self.out:
            return
        h = [x for x in self.s["history"] if x.get("out") != self.out]
        name = os.path.splitext(os.path.basename(self.structure or self.out))[0]
        h.insert(0, dict(name=name, structure=self.structure, out=self.out, status=status, time=time.time(),
                         where=self.where_label()))
        self.s["history"] = h[:50]
        save_settings(self.s)
        self.refresh_recent()

    def refresh_recent(self):
        for w in self.recent.winfo_children():
            w.destroy()
        for h in self.s["history"][:6]:
            row = tk.Frame(self.recent, bg=SIDE, cursor="hand2")
            row.pack(fill="x", pady=4)
            col = OK if status_of(h) == STATUS_DONE else (ACCENT if status_of(h) == STATUS_RUNNING else
                                                        (ERR if status_of(h) == STATUS_ERROR else FAINT))
            tk.Label(row, text="●", bg=SIDE, fg=col, font=pick_font(["Hiragino Sans"], 8)).pack(side="left", anchor="n",
                                                                                                padx=(6, 6), pady=3)
            box = tk.Frame(row, bg=SIDE)
            box.pack(side="left", fill="x")
            a = tk.Label(box, text=h.get("name", "?"), bg=SIDE, fg=TEXT, font=self.f_small, anchor="w")
            a.pack(anchor="w")
            b = tk.Label(box, text=f"{status_of(h)} · {ago(h.get('time'))}", bg=SIDE, fg=FAINT,
                         font=pick_font(["Hiragino Sans"], 10), anchor="w")
            b.pack(anchor="w")
            for w in (row, box, a, b):
                w.bind("<Button-1>", lambda e, h=h: self.open_history(h))

    def open_history(self, h):
        if self.busy():
            messagebox.showinfo("soapff-GUI", "A run is in progress")
            return
        if os.path.isfile(os.path.join(h.get("out", ""), "in.lammps")):
            self.load_project(h.get("structure", ""), h["out"])
        else:
            messagebox.showerror("soapff-GUI", f"Not found: {h.get('out')}")

    def show_history(self):
        self._nav_select("history")
        w = tk.Toplevel(self.root)
        w.title("Calculation history")
        w.geometry("900x420")
        tv = ttk.Treeview(w, columns=("status", "when", "where", "out"), show="tree headings")
        for col, txt, wd in (("#0", "Name", 200), ("status", "Status", 110), ("when", "Date", 140), ("where", "Run on", 170),
                             ("out", "Output folder", 400)):
            tv.heading(col, text=txt)
            tv.column(col, width=wd, anchor="w")
        for k, h in enumerate(self.s["history"]):
            tv.insert("", "end", iid=str(k), text=h.get("name", ""),
                      values=(status_of(h), time.strftime("%Y-%m-%d %H:%M", time.localtime(h.get("time", 0))),
                              where_of(h), h.get("out", "")))
        tv.pack(fill="both", expand=True, padx=10, pady=10)
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=10, pady=(0, 10))

        def chosen():
            s = tv.selection()
            return self.s["history"][int(s[0])] if s else None

        def open_():
            h = chosen()
            if h:
                w.destroy()
                self.open_history(h)

        def folder():
            h = chosen()
            if h and os.path.isdir(h.get("out", "")):
                open_folder(h["out"])

        def forget():
            h = chosen()
            if h:
                self.s["history"].remove(h)
                save_settings(self.s)
                tv.delete(tv.selection()[0])
                self.refresh_recent()
        ttk.Button(bar, text="Open", style="Accent.TButton", command=open_).pack(side="left")
        ttk.Button(bar, text="Open folder", style="Ghost.TButton", command=folder).pack(side="left", padx=6)
        ttk.Button(bar, text="Remove from history", style="Ghost.TButton", command=forget).pack(side="left")
        tv.bind("<Double-1>", lambda e: open_())
        w.bind("<Destroy>", lambda e: self._nav_select("new") if e.widget is w else None)

    def show_library(self):
        self._nav_select("library")
        w = tk.Toplevel(self.root)
        w.title("Structure library")
        w.geometry("760x520")
        txt = scrolledtext.ScrolledText(w, font=self.f_mono_s, height=20)
        txt.pack(fill="both", expand=True, padx=10, pady=10)

        def fill():
            txt.config(state="normal")
            txt.delete("1.0", "end")
            p = self.s["library"]
            txt.insert("end", f"Library: {p}\n\n")
            try:
                import numpy as np
                import collections
                d = np.load(p, allow_pickle=False)
                el, ty, q = d["element"], d["type"], d["q"]
                org = d["origin"] if "origin" in d.files else []
                txt.insert("end", f"{len(el)} rows (1 row = 1 atom of a reference structure: SOAP vector of {d['X'].shape[1]} components, charge, type)\n")
                txt.insert("end", f"Elements: {' '.join(d['species'])}\n\nReference structures\n")
                import re
                kept = collections.Counter(org)
                groups = collections.OrderedDict()
                for src in d["source"]:
                    base = os.path.basename(str(src))
                    stem = base[:-5] if base.endswith(".data") else base
                    mineral, _, tail = stem.rpartition("_")
                    n = kept.get(base, 0)
                    if tail == "bulk":
                        groups.setdefault(mineral, dict(bulk=None, basal=None, edges=[], other=[]))["bulk"] = n
                    elif re.fullmatch(r"[0-9m]+", tail):
                        g = groups.setdefault(mineral, dict(bulk=None, basal=None, edges=[], other=[]))
                        if tail == "001":
                            g["basal"] = n
                        else:
                            g["edges"].append((tail, n))
                    elif "particle" in stem:
                        groups.setdefault(stem, dict(bulk=None, basal=None, edges=[], other=[]))["other"].append(("particle", n))
                    else:                                       # a bulk reference without "_bulk" in its name
                        groups.setdefault(stem, dict(bulk=None, basal=None, edges=[], other=[]))["bulk"] = n
                txt.insert("end", "  Atoms kept in the library (an environment already present is not stored twice)\n"
                                  "  bulk: the crystal; 001: basal plane (no bond cut; the only kind of surface supported, no edges)\n\n")
                for m, g in groups.items():
                    if g["other"]:
                        txt.insert("end", f"  {m}: edge particle (JAEA), {g['other'][0][1]} atoms\n")
                        continue
                    txt.insert("end", f"  {m}: bulk {g['bulk']}" + ((", 001 (a face with bonds cut, no layers) " if m.startswith("quartz") else ", 001 ") + str(g["basal"]) if g["basal"] is not None else "")
                               + (f", edges {len(g['edges'])} faces ({sum(n for _, n in g['edges'])} atoms)" if g["edges"] else "") + "\n")
                    if g["edges"]:
                        txt.insert("end", "      " + " ".join(f + ("" if n else "*") for f, n in g["edges"]) + "\n")
                txt.insert("end", "\n  (* = every environment of that face was already in the library)\n")
                txt.insert("end", "\nRows and charges per type\n")
                for (e, t), n in sorted(collections.Counter(zip(el, ty)).items()):
                    qs = sorted({round(float(x), 4) for x, ee, tt in zip(q, el, ty) if ee == e and tt == t})
                    txt.insert("end", f"  {e:3s} {t:6s} {n:6d}   charge {qs}\n")
            except Exception as e:
                txt.insert("end", f"Cannot read: {e}\n")
            txt.config(state="disabled")

        def choose():
            p = filedialog.askopenfilename(title="Choose a library", filetypes=[("SOAP library", "*.npz")])
            if p:
                self.s["library"] = p
                save_settings(self.s)
                fill()
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=10, pady=(0, 10))
        ttk.Button(bar, text="Choose another library…", style="Ghost.TButton", command=choose).pack(side="left")
        ttk.Button(bar, text="Close", style="Ghost.TButton", command=w.destroy).pack(side="right")
        fill()
        w.bind("<Destroy>", lambda e: self._nav_select("new") if e.widget is w else None)

    # ================================================================================================== settings
    def open_settings(self):
        w = tk.Toplevel(self.root)
        w.title("Calculation settings")
        w.withdraw()                              # shown at the size of its largest tab (below)
        w.configure(bg=PANEL)
        nb = ttk.Notebook(w)
        nb.pack(fill="both", expand=True, padx=12, pady=12)
        v = {}

        def entry(parent, row, key, label, width=12, show=None, browse=None, note=None):
            ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4, padx=(0, 10))
            v[key] = tk.StringVar(value=str(self.s.get(key, "")))
            e = ttk.Entry(parent, textvariable=v[key], width=width, style="Field.TEntry", show=show)
            e.grid(row=row, column=1, sticky="w" if width < 30 else "ew", pady=4)
            if browse:
                ttk.Button(parent, text="Browse…", style="Ghost.TButton",
                           command=lambda: v[key].set((filedialog.askdirectory() if browse == "dir" else
                                                       filedialog.askopenfilename()) or v[key].get())).grid(row=row, column=2, padx=6)
            if note:
                ttk.Label(parent, text=note, foreground=MUTED).grid(row=row, column=3 if browse else 2, sticky="w", padx=6)
            return e

        # --- MD
        md = ttk.Frame(nb, padding=16)
        nb.add(md, text="MD conditions")
        for i, (k, lab, note) in enumerate((("temp", "Temperature (K)", ""), ("relax_ps", "Relaxation NVT (ps)", "After minimization; not recorded"),
                                            ("nvt_ps", "NVT (ps)", ""), ("npt_ps", "NPT (ps)", "0 = none"),
                                            ("press", "Pressure (atm)", ""), ("dump_ps", "Output interval (ps)", "0.1 to 0.5 for diffusion coefficients"))):
            entry(md, i, k, lab, note=note)
        ttk.Label(md, text="Pressure direction").grid(row=6, column=0, sticky="w", pady=4)
        v["baro"] = tk.StringVar(value=self.s["baro"])
        ttk.Combobox(md, textvariable=v["baro"], values=list(BARO), width=8, state="readonly").grid(row=6, column=1, sticky="w")
        ttk.Label(md, text="z: interlayer / interface; aniso / tri: bulk crystal", foreground=MUTED).grid(row=6, column=2, sticky="w", padx=6)
        ttk.Label(md, text="Boundary (x y z)").grid(row=7, column=0, sticky="w", pady=4)
        bf = ttk.Frame(md)
        bf.grid(row=7, column=1, columnspan=2, sticky="w")
        b0 = (self.s["boundary"].split() + ["p", "p", "p"])[:3]
        vb = []
        for k in range(3):
            x = tk.StringVar(value=b0[k])
            ttk.Combobox(bf, textvariable=x, values=AXIS_STYLES, width=3, state="readonly").pack(side="left", padx=(0, 4))
            vb.append(x)
        ttk.Label(md, text="These conditions are used the next time an input is built (New calculation, or the button below).\n"
                           "Once built, the input can be edited line by line in the main window.", foreground=MUTED, justify="left").grid(
            row=8, column=0, columnspan=4, sticky="w", pady=(14, 6))

        # --- where to run
        ex = ttk.Frame(nb, padding=16)
        nb.add(ex, text="Run on")
        ex.columnconfigure(0, weight=1)
        lo = ttk.LabelFrame(ex, text=LOCAL_NAME, padding=10)
        lo.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        lo.columnconfigure(1, weight=1)
        entry(lo, 0, "lmp", "LAMMPS (lmp)", width=46, browse="file")
        entry(lo, 1, "mpirun", "mpirun", width=46, browse="file")
        entry(lo, 2, "nproc", "Parallel processes", width=6, note="On a laptop, at most the number of P-cores")
        rm = ttk.LabelFrame(ex, text="Registered computers (SSH)", padding=10)
        rm.grid(row=1, column=0, sticky="nsew")
        ex.rowconfigure(1, weight=1)
        rm.columnconfigure(0, weight=1)
        comps = [dict(c) for c in self.s["computers"]]
        tv = ttk.Treeview(rm, columns=("host", "user", "port", "nproc"), show="tree headings", height=5)
        for col, txt, wd in (("#0", "Name", 170), ("host", "Host", 170), ("user", "User", 130), ("port", "Port", 70),
                             ("nproc", "Processes", 90)):
            tv.heading(col, text=txt)
            tv.column(col, width=wd, anchor="w")
        tv.grid(row=0, column=0, sticky="nsew")

        def refill():
            tv.delete(*tv.get_children())
            for k, c in enumerate(comps):
                tv.insert("", "end", iid=str(k), text=c["name"], values=(c["host"], c["user"], c["port"], c["nproc"]))

        def chosen():
            sel = tv.selection()
            return int(sel[0]) if sel else None

        def edit(k=None):
            if k is None:                                     # Add: register with only the ssh destination (connect and fill in the rest)
                c = self.quick_add_computer(w)
                if c:
                    comps.append(dict(c))
                    refill()
                return
            c = self.edit_computer(w, comps[k] if k is not None else dict(COMPUTER),
                                   taken=[x["name"] for n, x in enumerate(comps) if n != k])
            if c:
                if k is None:
                    comps.append(c)
                else:
                    comps[k] = c
                refill()

        def remove():
            k = chosen()
            if k is not None and messagebox.askyesno("soapff-GUI", f"Remove the registration of \"{comps[k]['name']}\"?", parent=w):
                del comps[k]
                refill()
        res = tk.StringVar(value="Passwords are not saved. You are asked at each connection test and at the start of each run.")

        def test():
            k = chosen()
            if k is None:
                res.set("Select the computer to test")
                return
            c = comps[k]
            pw = self.ask_password(c, parent=w)
            if pw is None:
                return
            res.set(f"Connecting to {c['name']}…")

            def go():
                r = None
                try:
                    r = self.make_remote(c, pw)
                    info = r.check()
                    vers = "  ".join(f"{k[15:-1]}: {v}" for k, v in info.items() if k.startswith("mpirun_version["))
                    msg = (f"Connected: {c['name']} = {info.get('host')} ({info.get('os')}, CPU {info.get('cpus')}, "
                           f"physical cores {info.get('cores') or '?'}. Keep the processes at or below the physical cores)\n"
                           f"lmp: {info.get('lmp')} ({info.get('lmp_mpi') or 'MPI unknown'})\n"
                           f"mpirun used: {info.get('mpirun_used')}  {vers}\nFolder: {info.get('root')}")
                    if info.get("mpi_warning"):
                        msg += "\n⚠ " + info["mpi_warning"]
                    if "NOTFOUND" in (info.get("lmp", "") + info.get("mpirun", "")):
                        msg += "\n⚠ lmp or mpirun not found. Set \"Remote lmp\" / \"Remote mpirun\" or \"Setup commands\" in the registration"
                except Exception as e:
                    msg = f"Could not connect: {e}"
                finally:
                    if r:
                        r.close()
                self.q.put(("call", lambda: res.set(msg)))
            threading.Thread(target=go, daemon=True).start()
        btns = ttk.Frame(rm)
        btns.grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Button(btns, text="Add…", style="Ghost.TButton", command=edit).pack(side="left")
        ttk.Button(btns, text="Edit…", style="Ghost.TButton", command=lambda: chosen() is not None and edit(chosen())).pack(
            side="left", padx=6)
        ttk.Button(btns, text="Remove", style="Ghost.TButton", command=remove).pack(side="left")
        ttk.Button(btns, text="Test connection", style="Ghost.TButton", command=test).pack(side="left", padx=6)
        ttk.Label(rm, textvariable=res, foreground=MUTED, wraplength=660, justify="left").grid(row=2, column=0, sticky="w",
                                                                                              pady=(8, 0))
        tv.bind("<Double-1>", lambda e: chosen() is not None and edit(chosen()))
        refill()

        # --- library and details
        de = ttk.Frame(nb, padding=16)
        nb.add(de, text="Library and details")
        de.columnconfigure(1, weight=1)
        entry(de, 0, "library", "SOAP library", width=50, browse="file")
        ttk.Button(de, text="Restore bundled", style="Ghost.TButton",
                   command=lambda: v["library"].set(runner.default_library())).grid(row=0, column=3, padx=(0, 6))
        entry(de, 1, "threshold", "Similarity threshold", width=8)
        entry(de, 2, "cutoff", "Cutoff (Å)", width=8)
        for row, key, lab, names, note in ((3, "charges", "Charges", CHARGES, "the type of each atom comes from SOAP, its charge from this scheme"),
                                           (4, "water", "Water model", WATERS, "SPC/E: fix shake1 keeps the water rigid after the minimization")):
            ttk.Label(de, text=lab).grid(row=row, column=0, sticky="w", pady=4, padx=(0, 10))
            v[key] = tk.StringVar(value=names.get(self.s.get(key), next(iter(names.values()))))
            ttk.Combobox(de, textvariable=v[key], values=list(names.values()), width=30, state="readonly").grid(row=row, column=1, sticky="w")
            ttk.Label(de, text=note, foreground=MUTED).grid(row=row, column=2, columnspan=2, sticky="w", padx=6)

        # --- analysis
        an = ttk.Frame(nb, padding=16)
        nb.add(an, text="Analysis")
        self.analysis_form(an, v, self.s, auto=True)

        def collect():
            s2 = dict(self.s)
            try:
                for k in ("temp", "relax_ps", "nvt_ps", "npt_ps", "press", "dump_ps", "threshold", "cutoff"):
                    s2[k] = float(v[k].get())
                s2["nproc"] = int(float(v["nproc"].get() or 1))
            except ValueError as e:
                raise ValueError(f"Not a number: {e}")
            if s2["dump_ps"] <= 0 or s2["temp"] <= 0:
                raise ValueError("Temperature and output interval must be positive")
            for k in ("baro", "lmp", "mpirun", "library"):
                s2[k] = v[k].get().strip()
            for k, names in (("charges", CHARGES), ("water", WATERS)):
                s2[k] = next(key for key, lab in names.items() if lab == v[k].get())
            self.collect_analysis(v, s2)
            s2["boundary"] = " ".join(x.get() for x in vb)
            s2["computers"] = [dict(c) for c in comps]
            if s2["where"] != "local" and not any(c["name"] == s2["where"] for c in comps):
                s2["where"] = "local"
            if s2["npt_ps"] > 0 and s2["boundary"] != "p p p":
                raise ValueError("NPT cannot be used when a direction is non-periodic")
            return s2

        def ok(regen=False):
            try:
                s2 = collect()
            except ValueError as e:
                messagebox.showerror("soapff-GUI", str(e), parent=w)
                return
            self.s = s2
            save_settings(self.s)
            self.update_status()
            w.destroy()
            if regen:
                self.regenerate()
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=12, pady=(0, 12))
        ttk.Button(bar, text="Save", style="Accent.TButton", command=ok).pack(side="right")
        ttk.Button(bar, text="Save and rebuild input", style="Ghost.TButton", command=lambda: ok(True)).pack(side="right", padx=8)
        ttk.Button(bar, text="Close", style="Ghost.TButton", command=w.destroy).pack(side="left")
        w.update_idletasks()                      # the largest tab decides the size; not beyond the screen
        w.geometry(f"{min(w.winfo_reqwidth(), w.winfo_screenwidth() - 40)}x{min(w.winfo_reqheight(), w.winfo_screenheight() - 90)}")
        w.deiconify()

    def regenerate(self):
        if not self.structure or not os.path.isfile(self.structure):
            messagebox.showinfo("soapff-GUI", "The structure file is unknown. Choose one with \"New calculation\"")
            return
        if self.busy():
            return
        if self.lines and not messagebox.askyesno("soapff-GUI", "Rebuild the input from the MD conditions. Lines you edited by hand are lost. Continue?"):
            return
        self.make_input(self.structure, self.out or os.path.splitext(self.structure)[0] + "_soapff")

    # ================================================================================================== run
    def fetch_dialog(self):
        """copy a folder (with everything in it, like cp -r) from a registered computer to this one."""
        if not self.s["computers"]:
            messagebox.showinfo("soapff-GUI", "First register a computer under \"Run on\" in the Calculation settings")
            return
        w = tk.Toplevel(self.root)
        w.title("Fetch a folder from another computer")
        w.geometry("820x560")
        w.configure(bg=PANEL)
        f = ttk.Frame(w, padding=14)
        f.pack(fill="both", expand=True)
        f.columnconfigure(1, weight=1)
        names = [c["name"] for c in self.s["computers"]]
        v_c = tk.StringVar(value=self.s["where"] if self.s["where"] in names else names[0])
        v_path = tk.StringVar(value=(self.computer(v_c.get()) or {}).get("root", "soapff_runs"))
        v_dest = tk.StringVar(value=os.path.dirname(self.out) if self.out else self.s["last_dir"])
        ttk.Label(f, text="Computer").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Combobox(f, textvariable=v_c, values=names, state="readonly", width=30).grid(row=0, column=1, sticky="w")
        ttk.Label(f, text="Remote folder").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(f, textvariable=v_path, style="Field.TEntry").grid(row=1, column=1, sticky="ew")
        lst = tk.Listbox(f, font=self.f_mono_s, height=14, activestyle="none", selectbackground=SEL_BG, selectforeground=TEXT)
        lst.grid(row=3, column=0, columnspan=3, sticky="nsew", pady=(6, 6))
        f.rowconfigure(3, weight=1)
        ttk.Label(f, text="Destination (this computer)").grid(row=4, column=0, sticky="w", pady=4)
        ttk.Entry(f, textvariable=v_dest, style="Field.TEntry").grid(row=4, column=1, sticky="ew")
        ttk.Button(f, text="Browse…", style="Ghost.TButton",
                   command=lambda: v_dest.set(filedialog.askdirectory(parent=w) or v_dest.get())).grid(row=4, column=2, padx=6)
        res = tk.StringVar(value="\"List\" shows the contents of the remote folder. Double-click a folder to enter it.")
        ttk.Label(f, textvariable=res, foreground=MUTED, wraplength=760, justify="left").grid(row=5, column=0, columnspan=3,
                                                                                            sticky="w", pady=(4, 8))
        st = {"r": None, "items": []}

        def remote_():
            c = self.computer(v_c.get())
            if st["r"] is None or st.get("name") != c["name"]:
                if st["r"]:
                    st["r"].close()
                pw = self.ask_password(c, parent=w)
                if pw is None:
                    return None
                st["r"], st["name"] = self.make_remote(c, pw), c["name"]
            return st["r"]

        def listing():
            r = remote_()
            if not r:
                return
            res.set("Fetching the list…")
            path = v_path.get().strip() or "."

            def go():
                try:
                    items = r.list_dir(path)
                    msg = f"{path}: {sum(d for _, d in items)} folders, {sum(not d for _, d in items)} files"
                except Exception as e:
                    items, msg = [], f"Could not get the list: {e}"

                def show():
                    st["items"] = items
                    lst.delete(0, "end")
                    for n, d in items:
                        lst.insert("end", ("📁 " if d else "    ") + n)
                    res.set(msg)
                self.q.put(("call", show))
            threading.Thread(target=go, daemon=True).start()

        def enter(e=None):
            k = lst.curselection()
            if k and st["items"][k[0]][1]:
                v_path.set(v_path.get().rstrip("/") + "/" + st["items"][k[0]][0])
                listing()

        def up():
            p = v_path.get().rstrip("/")
            v_path.set(p.rsplit("/", 1)[0] if "/" in p else ".")
            listing()

        def fetch():
            r = remote_()
            if not r:
                return
            k = lst.curselection()
            src = v_path.get().rstrip("/") + ("/" + st["items"][k[0]][0] if k else "")
            dest = v_dest.get().strip()
            if not dest:
                res.set("Choose a destination")
                return
            target = os.path.join(dest, os.path.basename(src.rstrip("/")))
            if os.path.exists(target) and not messagebox.askyesno("soapff-GUI", f"{target} already exists. Overwrite?",
                                                                  parent=w):
                return
            res.set(f"Fetching {src} (scp -r)…")

            def go():
                try:
                    got = r.fetch(src, dest)
                    n = sum(len(fs) for _, _, fs in os.walk(got)) if os.path.isdir(got) else 1
                    msg = f"Fetched: {got} ({n} files)"
                    if os.path.isfile(os.path.join(got, "in.lammps")):
                        msg += ". Open this folder with \"File → Open existing output folder\""
                except Exception as e:
                    msg = f"Could not fetch: {e}"
                self.q.put(("call", lambda: res.set(msg)))
            threading.Thread(target=go, daemon=True).start()

        def close():
            if st["r"]:
                st["r"].close()
            w.destroy()
        bar = ttk.Frame(f)
        bar.grid(row=2, column=1, sticky="w")
        ttk.Button(bar, text="List", style="Ghost.TButton", command=listing).pack(side="left")
        ttk.Button(bar, text="Up one folder", style="Ghost.TButton", command=up).pack(side="left", padx=6)
        lst.bind("<Double-1>", enter)
        b2 = ttk.Frame(f)
        b2.grid(row=6, column=0, columnspan=3, sticky="ew")
        ttk.Button(b2, text="Fetch (the whole selected folder, or the folder above if none is selected)", style="Accent.TButton",
                   command=fetch).pack(side="right")
        ttk.Button(b2, text="Close", style="Ghost.TButton", command=close).pack(side="left")
        w.protocol("WM_DELETE_WINDOW", close)

    def computer(self, name=None):
        name = name or self.s.get("where")
        return next((c for c in self.s["computers"] if c["name"] == name), None)

    def where_label(self):
        c = self.computer() if self.s.get("where") != "local" else None
        if c:
            return f"{c['name']} ({c['user'] + '@' if c['user'] else ''}{c['host']})"
        return LOCAL_NAME

    def choose_where(self):
        v = self.v_where.get()
        self.s["where"] = "local" if is_local_name(v) else v
        save_settings(self.s)
        self.update_status()

    @staticmethod
    def parse_target(text):
        """'user@host', 'user@host:port', 'ssh -p 10022 user@host', 'user@host -p 10022' -> (user, host, port)."""
        w = text.replace("ssh ", " ").split()
        port, rest = 22, []
        k = 0
        while k < len(w):
            if w[k] == "-p" and k + 1 < len(w):
                port = int(w[k + 1]); k += 2
                continue
            rest.append(w[k]); k += 1
        if len(rest) != 1:
            raise ValueError("Enter it as user@host (e.g. taro@192.168.0.10; port with :10022 or -p 10022)")
        t = rest[0]
        user, _, host = t.rpartition("@")
        if ":" in host and host.count(":") == 1:
            host, p_ = host.split(":"); port = int(p_)
        if not host:
            raise ValueError("The host is missing")
        return user, host, port

    def quick_add_computer(self, parent):
        """register a computer from 'user@host' only: asks the password, connects, and fills in the rest (physical cores as
        the number of processes, the lmp and mpirun found there).  Details can be changed later (Calculation settings › Run on)."""
        if self.busy():
            messagebox.showinfo("soapff-GUI", "Cannot register during a run", parent=parent)
            return None
        w = tk.Toplevel(parent)
        w.title("Register a computer (SSH)")
        w.transient(parent)
        f = ttk.Frame(w, padding=16)
        f.pack(fill="both", expand=True)
        ttk.Label(f, text="ssh destination", font=self.f_bold).grid(row=0, column=0, sticky="w")
        v_t = tk.StringVar()
        e = ttk.Entry(f, textvariable=v_t, width=40, style="Field.TEntry")
        e.grid(row=1, column=0, sticky="w", pady=(4, 2))
        e.focus_set()
        ttk.Label(f, text="Example: taro@192.168.0.10  (if the port is not 22: taro@192.168.0.10:10022)", foreground=MUTED).grid(
            row=2, column=0, sticky="w")
        res = tk.StringVar(value="Pressing Register asks for the password, connects, and checks the remote LAMMPS and CPU (the password is not saved).")
        ttk.Label(f, textvariable=res, foreground=MUTED, wraplength=460, justify="left").grid(row=3, column=0, sticky="w",
                                                                                            pady=(10, 0))
        done = {}

        def register():
            try:
                user, host, port = self.parse_target(v_t.get().strip())
            except ValueError as ex:
                res.set(str(ex)); return
            name = f"{user + '@' if user else ''}{host}" + (f":{port}" if port != 22 else "")
            if any(c["name"] == name for c in self.s["computers"]):
                res.set(f"\"{name}\" is already registered"); return
            c = dict(COMPUTER, name=name, host=host, user=user, port=port)
            pw = self.ask_password(c, parent=w)
            if pw is None:
                return
            res.set(f"Connecting to {name}…")
            b_ok.state(["disabled"])

            def go():
                r = None
                try:
                    r = self.make_remote(c, pw)
                    info = r.check()
                    cores = str(info.get("cores") or "")
                    c["nproc"] = int(cores) if cores.isdigit() and int(cores) > 0 else 1
                    if info.get("mpirun_used") and info.get("mpirun_used") != "NOTFOUND":
                        c["mpirun"] = info["mpirun_used"]
                    msg = (f"Registered: {name} ({info.get('os')}, physical cores {info.get('cores') or '?'} → processes {c['nproc']})\n"
                           f"lmp: {info.get('lmp')}  mpirun: {c['mpirun']}")
                    if "NOTFOUND" in (info.get("lmp", "") + info.get("mpirun", "")):
                        msg += ("\n⚠ lmp or mpirun not found. In Calculation settings › Run on, set \"Remote lmp\" / \"Remote mpirun\" or "
                                "\"Setup commands\" (e.g. module load lammps)")
                    if info.get("mpi_warning"):
                        msg += "\n⚠ " + info["mpi_warning"]
                    ok = True
                except Exception as ex:
                    msg, ok = f"Could not connect: {ex}", False
                finally:
                    if r:
                        r.close()

                def finish():
                    res.set(msg)
                    b_ok.state(["!disabled"])
                    if ok:
                        self.s["computers"].append(c)
                        self.s["where"] = c["name"]
                        save_settings(self.s)
                        self.update_status()
                        done["c"] = c
                        b_ok.configure(text="Close", command=w.destroy)
                self.q.put(("call", finish))
            threading.Thread(target=go, daemon=True).start()
        bar = ttk.Frame(f)
        bar.grid(row=4, column=0, sticky="ew", pady=(14, 0))
        b_ok = ttk.Button(bar, text="Register", style="Accent.TButton", command=register)
        b_ok.pack(side="right")
        ttk.Button(bar, text="Cancel", style="Ghost.TButton", command=w.destroy).pack(side="left")
        w.bind("<Return>", lambda ev: register() if b_ok.cget("text") == "Register" else w.destroy())
        w.grab_set()
        w.wait_window()
        return done.get("c")

    def remove_computer(self):
        """remove the computer chosen as "Run on" from the registered ones."""
        name = self.v_where.get()
        if is_local_name(name) or not self.computer(name):
            messagebox.showinfo("soapff-GUI", "Choose the computer to remove under \"Run on\" first, then press the button")
            return
        if self.busy():
            messagebox.showinfo("soapff-GUI", "Cannot remove during a run")
            return
        if not messagebox.askyesno("soapff-GUI", f"Remove the registration of \"{name}\"? (Files on the remote side are not deleted)"):
            return
        self.s["computers"] = [c for c in self.s["computers"] if c["name"] != name]
        self.s["where"] = "local"
        save_settings(self.s)
        self.update_status()

    def edit_computer(self, parent, c, taken=()):
        """a dialog to register / edit a computer; returns the new dict or None."""
        w = tk.Toplevel(parent)
        w.title("Register a computer")
        w.configure(bg=PANEL)
        w.transient(parent)
        f = ttk.Frame(w, padding=16)
        f.pack(fill="both", expand=True)
        v = {}
        rows = (("name", "Name", "Name used within this app (e.g. lab server)"), ("host", "Host", "e.g. 192.168.0.10 or a host name"),
                ("user", "User", ""), ("port", "Port", "SSH port (default 22)"),
                ("root", "Remote folder", "Folder for the runs (a path relative to home is OK)"),
                ("lmp", "Remote lmp", "Name or full path"), ("mpirun", "Remote mpirun", "Name or full path"),
                ("setup", "Setup commands", "e.g. module load lammps, source ~/.bashrc (run first)"), ("nproc", "Processes", "At most the remote physical cores"))
        for i, (k, lab, note) in enumerate(rows):
            ttk.Label(f, text=lab).grid(row=i, column=0, sticky="w", pady=4, padx=(0, 10))
            v[k] = tk.StringVar(value=str(c.get(k, "")))
            ttk.Entry(f, textvariable=v[k], width=36 if k not in ("port", "nproc") else 8, style="Field.TEntry").grid(
                row=i, column=1, sticky="w", pady=4)
            ttk.Label(f, text=note, foreground=MUTED).grid(row=i, column=2, sticky="w", padx=8)
        ttk.Label(f, text="The password is not stored. It is asked at each connection.", foreground=MUTED).grid(
            row=len(rows), column=0, columnspan=3, sticky="w", pady=(10, 0))
        out = {}

        def ok():
            d = {k: v[k].get().strip() for k in v}
            try:
                d["port"] = int(d["port"] or 22)
                d["nproc"] = int(d["nproc"] or 1)
            except ValueError:
                messagebox.showerror("soapff-GUI", "Port and processes must be integers", parent=w)
                return
            if not d["name"] or not d["host"]:
                messagebox.showerror("soapff-GUI", "Enter a name and a host", parent=w)
                return
            if d["name"] in taken or is_local_name(d["name"]):
                messagebox.showerror("soapff-GUI", f"\"{d['name']}\" is already in use", parent=w)
                return
            out.update(dict(COMPUTER, **d))
            w.destroy()
        bar = ttk.Frame(f)
        bar.grid(row=len(rows) + 1, column=0, columnspan=3, sticky="e", pady=(14, 0))
        ttk.Button(bar, text="Cancel", style="Ghost.TButton", command=w.destroy).pack(side="left", padx=6)
        ttk.Button(bar, text="Register", style="Accent.TButton", command=ok).pack(side="left")
        w.grab_set()
        w.wait_window()
        return out or None

    def update_status(self):
        r = self.report
        st = f"Structure: {os.path.basename(self.structure) if self.structure else '—'}"
        if r:
            st += f"   |   {r['atoms']:,} atoms"
        self.l_status.config(text=st)
        self.l_status_r.config(text=f"Run on: {self.where_label()}   ·   LAMMPS · SOAP")
        names = [LOCAL_NAME] + [c["name"] for c in self.s["computers"]]
        self.cb_where.configure(values=names)
        self.v_where.set(self.s["where"] if self.s.get("where") in names else LOCAL_NAME)

    def busy(self):
        return bool(self.worker and self.worker.is_alive())

    def start_or_stop(self):
        if self.busy():
            self.stop()
        else:
            self.start()

    def ask_password(self, c, parent=None):
        """the SSH password of a registered computer, asked every time (never stored).  None: cancelled."""
        return simpledialog.askstring("SSH", f"Password for {c['name']} ({c['user'] + '@' if c['user'] else ''}{c['host']}:{c['port']})",
                                             show="*", parent=parent or self.root)

    def make_remote(self, c, password=None):
        return remote.Remote(c["host"], c.get("user", ""), c.get("port", 22), "", c.get("root", "soapff_runs"),
                             c.get("lmp", "lmp"), c.get("mpirun", "mpirun"), c.get("setup", ""), password=password or None)

    def start(self):
        if self.busy():
            return
        if not self.docs["in.lammps"] or not self.out:
            messagebox.showinfo("soapff-GUI", "First choose a structure file with \"New calculation\" and create the input")
            return
        if not L.run_steps(self.docs["in.lammps"]) and not messagebox.askyesno(
                "soapff-GUI", "in.lammps has no run (energy minimization only). Run it anyway?\n"
                              "To run MD, add a \"+ fix\" (such as nvt) and a run from \"+ Other\"."):
            return
        if self._save_after:
            self.root.after_cancel(self._save_after)
        self.save_input()
        s = self.s
        comp = self.computer() if s.get("where") != "local" else None
        ssh = comp is not None
        if ssh:
            pw = self.ask_password(comp)
            if pw is None:
                return
        else:
            lmp = runner.which(s["lmp"])
            if not lmp:
                messagebox.showerror("soapff-GUI", "LAMMPS (lmp) not found. Specify its location in the Calculation settings")
                return
            mpirun = runner.which(s["mpirun"]) if int(s["nproc"]) > 1 else None
            if int(s["nproc"]) > 1 and not mpirun:
                messagebox.showerror("soapff-GUI", "mpirun not found. Set the number of processes to 1 or specify its location in the Calculation settings")
                return
        self.stopping = False
        self.log_clear()
        self.set_state("Starting", ACCENT, 0)
        self.set_running(True)
        self.remember(STATUS_RUNNING)
        lines = list(self.docs["in.lammps"])
        if ssh:
            self.worker = threading.Thread(target=self.work_remote, args=(dict(s), dict(comp), pw, lines), daemon=True)
        else:
            self.worker = threading.Thread(target=self.work_local, args=(dict(s), lmp, mpirun, lines), daemon=True)
        self.worker.start()

    def _follow(self, it, tracker, put, nproc=1):
        warned = False
        for line in it:
            if not warned and nproc > 1 and "1 by 1 by 1 MPI processor grid" in line:
                warned = True
                put(("log", f"WARNING Expected {nproc} parallel processes but the grid is '1 by 1 by 1': {nproc} separate single-process LAMMPS runs are going. "
                            "Use the mpirun of the same MPI as lmp (give the full path of \"Remote mpirun\" in the Run on settings, "
                            "from the same folder as lmp)\n", "warn"))
            put(("log", line + "\n", "err" if line.startswith("ERROR") else ("warn" if line.startswith("WARNING") else None)))
            if tracker.feed(line):
                put(("thermo", tracker.phase, dict(tracker.values), tracker.progress))

    def work_local(self, s, lmp, mpirun, lines):
        put = self.q.put
        try:
            self.job = runner.LammpsJob(self.out, lmp, mpirun, int(s["nproc"]))
            put(("log", "$ " + " ".join(self.job.cmd) + f"\n  (in {self.out})\n", None))
            self._follow(self.job.lines(), L.StepTracker(lines), put, int(s.get("nproc", 1)))
            rc = self.job.returncode
            self._finish(rc, s, put)
        except Exception as e:
            put(("log", f"ERROR {type(e).__name__}: {e}\n", "err"))
            put(("done", f"Error: {e}", 1))
        finally:
            self.job = None

    def work_remote(self, s, comp, password, lines):
        put = self.q.put
        name = os.path.basename(os.path.normpath(self.out))
        r = None
        try:
            r = self.make_remote(comp, password)
            put(("phase", f"Connecting to {comp['name']}"))
            info = r.check()
            put(("log", f"Connected: {info.get('host')} ({info.get('os')}, CPU {info.get('cpus')}, physical cores {info.get('cores') or '?'})\n", None))
            put(("log", f"lmp: {info.get('lmp')} ({info.get('lmp_mpi') or 'MPI unknown'})\nmpirun: {info.get('mpirun_used')}\n", None))
            if info.get("mpi_warning") and int(comp.get("nproc", 1)) > 1:
                put(("log", "WARNING " + info["mpi_warning"] + "\n", "warn"))
            rdir = r.upload(self.out, name)
            put(("log", f"Input sent: {r.target}:{rdir}\n", None))
            nproc = int(comp.get("nproc", 1))
            cores = int(info.get("cores") or 0) if str(info.get("cores") or "").isdigit() else 0
            cpus = int(info.get("cpus") or 0) if str(info.get("cpus") or "").isdigit() else 0
            if cpus and nproc > cpus:
                put(("log", f"WARNING The number of processes {nproc} exceeds the remote CPU (thread) count {cpus}, so {cpus} processes are used\n", "warn"))
                nproc = cpus
            elif cores and nproc > cores:
                put(("log", f"The number of processes {nproc} exceeds the physical cores {cores}, so hardware threads are used too "
                            "(--use-hwthread-cpus in Open MPI; LAMMPS does not get much faster beyond the core count)\n", None))
            comp = dict(comp, nproc=nproc)
            pid = r.start(name, nproc)
            self.rjob = (r, name, pid)
            put(("log", f"LAMMPS started ({comp['name']}, pid {pid}, {comp.get('nproc')} processes). "
                        "The run continues even if the connection drops\n", None))
            put(("phase", "MD (remote)"))
            self._follow(r.follow(name, pid), L.StepTracker(lines), put, int(comp.get("nproc", 1)))
            rc = r.exit_code(name)
            put(("phase", "Fetching results"))
            got = r.download_all(name, self.out)
            put(("log", f"Fetched the whole remote folder ({len(got)} files): {', '.join(got[:12])}"
                        f"{' …' if len(got) > 12 else ''}\n", None))
            self._finish(rc if rc is not None else (0 if self.stopping else 1), s, put)
        except Exception as e:
            put(("log", f"ERROR {type(e).__name__}: {e}\n", "err"))
            put(("done", f"Error: {e}", 1))
        finally:
            self.rjob = None
            if r:
                r.close()                 # the password is asked again next time

    def _finish(self, rc, s, put):
        if self.stopping:
            put(("done", "Stopped", rc if rc else 1))
        elif rc == 0:
            if s.get("auto_analyze") and runner.analyzable(self.out):
                self._analyze(put, s)
            put(("done", "Done", 0))
        else:
            put(("done", f"LAMMPS terminated abnormally (exit code {rc}). Check the log and log.lammps", rc or 1))

    def stop(self):
        if not self.busy():
            return
        self.stopping = True
        self.l_phase.config(text="Stopping…")
        if self.job:
            self.job.stop()
        if self.rjob:
            r, name, pid = self.rjob
            threading.Thread(target=lambda: r.stop(pid), daemon=True).start()

    # ================================================================================================== analysis
    def analysis_form(self, f, v, s, auto=False):
        """the analysis choices and their parameters (settings tab and the Analyze dialog); values in v."""
        f.columnconfigure(1, weight=1)
        row = 0
        if auto:
            v["auto_analyze"] = tk.BooleanVar(value=bool(s.get("auto_analyze", True)))
            ttk.Checkbutton(f, text="Analyze when the calculation finishes", variable=v["auto_analyze"]).grid(row=row, column=0, columnspan=4,
                                                                                               sticky="w", pady=(0, 6))
            row += 1
        ttk.Label(f, text="Analysis", style="Section.TLabel").grid(row=row, column=0, columnspan=4, sticky="w")
        row += 1
        for k, lab in ANALYSES:
            v["an_" + k] = tk.BooleanVar(value=k in s.get("analyses", []))
            ttk.Checkbutton(f, text=lab, variable=v["an_" + k]).grid(row=row, column=0, columnspan=4, sticky="w", padx=(8, 0))
            row += 1

        def ent(key, label, width=10, browse=None, note=None):
            nonlocal row
            ttk.Label(f, text=label).grid(row=row, column=0, sticky="w", pady=3, padx=(0, 10))
            v[key] = tk.StringVar(value=str(s.get(key, "")))
            rf = ttk.Frame(f)                 # entry, button and note side by side: a wide entry does not push other notes
            rf.grid(row=row, column=1, columnspan=3, sticky="ew", pady=3)
            ttk.Entry(rf, textvariable=v[key], width=width, style="Field.TEntry").pack(side="left", fill="x",
                                                                                      expand=width >= 30)
            if browse:
                ttk.Button(rf, text="Browse…", style="Ghost.TButton",
                           command=lambda: v[key].set((filedialog.askdirectory() if browse == "dir" else
                                                       filedialog.askopenfilename()) or v[key].get())).pack(side="left", padx=6)
            if note:
                ttk.Label(rf, text=note, foreground=MUTED).pack(side="left", padx=6)
            row += 1
        ent("skip_ps", "Time to exclude as equilibration (ps)", 8, note="If empty, the first 20 %")
        ttk.Label(f, text="SOVA (PDF, RDF, S(Q))", style="Section.TLabel").grid(row=row, column=0, columnspan=4, sticky="w",
                                                                              pady=(10, 0))
        row += 1
        sv = s.get("sova") or DEFAULTS["sova"]

        def sent(key, label, width=8, note=None):
            nonlocal row
            ttk.Label(f, text=label).grid(row=row, column=0, sticky="w", pady=3, padx=(0, 10))
            v["sova_" + key] = tk.StringVar(value=str(sv.get(key, "")))
            rf = ttk.Frame(f)
            rf.grid(row=row, column=1, columnspan=3, sticky="ew", pady=3)
            ttk.Entry(rf, textvariable=v["sova_" + key], width=width, style="Field.TEntry").pack(side="left")
            if note:
                ttk.Label(rf, text=note, foreground=MUTED).pack(side="left", padx=6)
            row += 1
        sent("dr", "g(r) bin width dr (Å)", note="Smaller is finer but noisier")
        sent("rmax", "Maximum r (Å)", note="Limited to half the box width")
        ttk.Label(f, text="Q range (Å⁻¹)").grid(row=row, column=0, sticky="w", pady=3)
        rq = ttk.Frame(f)
        rq.grid(row=row, column=1, columnspan=3, sticky="w")
        for k, lab in (("qmin", "from"), ("qmax", "to"), ("dq", "step")):
            v["sova_" + k] = tk.StringVar(value=str(sv.get(k, "")))
            ttk.Entry(rq, textvariable=v["sova_" + k], width=6, style="Field.TEntry").pack(side="left")
            ttk.Label(rq, text=lab, foreground=MUTED).pack(side="left", padx=(3, 10))
        row += 1
        ttk.Label(f, text="S(Q) radiation").grid(row=row, column=0, sticky="w", pady=3)
        v["sova_radiation"] = tk.StringVar(value=sv.get("radiation", "xray"))
        rr = ttk.Frame(f)
        rr.grid(row=row, column=1, columnspan=3, sticky="w")
        for k, lab in (("xray", "X-ray"), ("neutron", "Neutron")):
            ttk.Radiobutton(rr, text=lab, value=k, variable=v["sova_radiation"]).pack(side="left", padx=(0, 12))
        row += 1
        sent("max_frames", "Frames to average", 6, note="Chosen at equal intervals after equilibration")
        sent("type_pairs", "Type pairs (optional)", 30, note="e.g. K-ob K-obts (g(r) split by ClayFF type; empty = none)")
        ttk.Label(f, text="Pairs to compute").grid(row=row, column=0, sticky="w", pady=3)
        v["sova_pairs"] = tk.StringVar(value=sv.get("pairs", "water_ion"))
        v["sova_all"] = tk.BooleanVar(value=bool(sv.get("all_elements", True)))
        rp = ttk.Frame(f)
        rp.grid(row=row, column=1, columnspan=3, sticky="w")
        ttk.Radiobutton(rp, text="Water O and ion pairs (Ow–Ow, ion–Ow, ion–mineral O)", value="water_ion",
                        variable=v["sova_pairs"]).pack(anchor="w")
        ttk.Radiobutton(rp, text="None", value="none", variable=v["sova_pairs"]).pack(anchor="w")
        ttk.Checkbutton(rp, text="Element pairs of all atoms (S(Q) is also computed from these)", variable=v["sova_all"]).pack(anchor="w")
        row += 1

    @staticmethod
    def collect_analysis(v, s2):
        """the values of analysis_form into the settings s2 (ValueError when a number is wrong)."""
        if "auto_analyze" in v:
            s2["auto_analyze"] = bool(v["auto_analyze"].get())
        s2["analyses"] = [k for k, _ in ANALYSES if v["an_" + k].get()]
        s2["skip_ps"] = v["skip_ps"].get().strip()
        try:
            if s2["skip_ps"]:
                float(s2["skip_ps"])
            s2["sova"] = dict(dr=float(v["sova_dr"].get()), rmax=float(v["sova_rmax"].get()), qmin=float(v["sova_qmin"].get()),
                              qmax=float(v["sova_qmax"].get()), dq=float(v["sova_dq"].get()),
                              radiation=v["sova_radiation"].get(), max_frames=max(1, int(float(v["sova_max_frames"].get()))),
                              pairs=v["sova_pairs"].get(), all_elements=bool(v["sova_all"].get()),
                              type_pairs=v["sova_type_pairs"].get().strip())
            sv_ = s2["sova"]
            if not (sv_["dr"] > 0 and sv_["rmax"] > 0 and 0 < sv_["qmin"] < sv_["qmax"] and sv_["dq"] > 0):
                raise ValueError("SOVA dr, maximum r or Q range is invalid")
            if sv_["pairs"] == "none" and not sv_["all_elements"] and not sv_["type_pairs"]:
                raise ValueError("Select at least one pair type for SOVA")
        except ValueError as e:
            raise ValueError(f"Not a number: {e}")

    def _analyze(self, put, s):
        what = [k for k, _ in ANALYSES if k in s.get("analyses", [])]
        put(("phase", "Analyzing"))
        put(("log", "\nAnalysis: " + ", ".join(lab.split(" (")[0] for k, lab in ANALYSES if k in what) + "\n", None))
        if not what:
            put(("log", "No analysis is selected (Calculation settings › Analysis)\n", "warn"))
            return
        try:
            skip = s.get("skip_ps")
            runner.analyze(self.out, log=lambda m: put(("log", m + "\n", None)), skip_ps=float(skip) if skip else None,
                           sova_opts=dict(s.get("sova") or {}), what=what,
                           debyer=s.get("debyer"), traj=s.get("_traj"),
                           xrd_opts=dict(wavelength=float(s.get("xrd_wavelength", 1.5406)),
                                         two_theta=(float(s.get("xrd_from", 3)), float(s.get("xrd_to", 80)),
                                                    float(s.get("xrd_step", 0.02))),
                                         frames=int(s.get("xrd_frames", 3)), mode=s.get("xrd_mode", "cluster"),
                                         size_A=float(s.get("xrd_size", 50))))
            put(("analysis", self.out))
        except Exception as e:
            put(("log", f"Analysis failed: {e}\n", "warn"))

    def analyze_now(self):
        if not runner.analyzable(self.out):
            messagebox.showinfo("soapff-GUI", "There is no output to analyze yet (a trajectory *.lammpstrj or an MSD file). "
                                              "Run the calculation first")
            return
        if self.busy():
            messagebox.showinfo("soapff-GUI", "A run is in progress. Analyze after it finishes")
            return
        w = tk.Toplevel(self.root)
        w.title(f"Analysis – {os.path.basename(self.out)}")
        w.configure(bg=PANEL)
        w.transient(self.root)
        from soapff_cui import analysis as AN
        trajs = AN.trajectories(self.out)
        v = {}
        if trajs:
            tf = ttk.Frame(w, padding=(16, 14, 16, 0))
            tf.pack(fill="x")
            ttk.Label(tf, text="Trajectory for SOVA and XRD").pack(side="left", padx=(0, 10))
            v["traj"] = tk.StringVar(value=trajs[0])
            ttk.Combobox(tf, textvariable=v["traj"], values=trajs, state="readonly", width=34).pack(side="left")
            ttk.Label(tf, text="(tensile / shear trajectories are structures in the middle of the deformation)", foreground=MUTED).pack(side="left", padx=8)
        f = ttk.Frame(w, padding=16)
        f.pack(fill="both", expand=True)
        self.analysis_form(f, v, self.s)

        def go():
            s2 = dict(self.s)
            try:
                self.collect_analysis(v, s2)
            except ValueError as e:
                messagebox.showerror("soapff-GUI", str(e), parent=w)
                return
            self.s = s2
            save_settings(self.s)
            s2 = dict(s2, _traj=v["traj"].get() if "traj" in v else None)
            w.destroy()
            self.worker = threading.Thread(target=lambda: (self._analyze(self.q.put, s2), self.q.put(("idle",))), daemon=True)
            self.set_running(True)
            self.worker.start()
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=16, pady=(0, 14))
        ttk.Button(bar, text="Analyze", style="Accent.TButton", command=go).pack(side="right")
        if os.path.isfile(os.path.join(self.out, "analysis", "summary.txt")):
            ttk.Button(bar, text="View previous results", style="Ghost.TButton",
                       command=lambda: (w.destroy(), self.show_analysis(self.out))).pack(side="right", padx=8)
        ttk.Button(bar, text="Close", style="Ghost.TButton", command=w.destroy).pack(side="left")

    def show_results(self):
        if self.out and os.path.isfile(os.path.join(self.out, "analysis", "summary.txt")):
            self.show_analysis(self.out)
        elif runner.analyzable(self.out):
            self.analyze_now()
        else:
            messagebox.showinfo("soapff-GUI", "There are no results yet. Run the calculation with \"Start calculation\"")

    def show_analysis(self, out):
        """the summary and every figure of out/analysis in tabs."""
        dest = os.path.join(out, "analysis")
        try:
            figs = [tuple(x) for x in json.load(open(os.path.join(dest, "summary.json"))).get("figures", [])]
        except (OSError, ValueError):
            figs = []
        if not figs and os.path.isfile(os.path.join(dest, "analysis.png")):     # analysis of an earlier version
            figs = [("MSD, PDF, RDF", "analysis.png")]
        w = tk.Toplevel(self.root)
        w.title(f"Analysis results – {os.path.basename(out)}")
        w.configure(bg=PANEL)
        nb = ttk.Notebook(w)
        nb.pack(fill="both", expand=True, padx=8, pady=(8, 4))
        tab = ttk.Frame(nb, padding=4)
        nb.add(tab, text="Summary")
        txt = scrolledtext.ScrolledText(tab, height=26, width=118, font=self.f_mono_s, wrap="none")
        txt.insert("end", open(os.path.join(dest, "summary.txt")).read())
        txt.config(state="disabled")
        txt.pack(fill="both", expand=True)
        pngs = [None]
        for title, png in figs:
            tab = ttk.Frame(nb, padding=4)
            nb.add(tab, text=title)
            path = os.path.join(dest, png)
            pngs.append(path)
            try:
                img = tk.PhotoImage(file=path)
                k = max(1, -(-img.width() // 1100), -(-img.height() // 700))
                img = img.subsample(k, k) if k > 1 else img
                lab = ttk.Label(tab, image=img)
                lab.image = img
                lab.pack()
            except tk.TclError as e:
                ttk.Label(tab, text=f"Cannot display the figure: {e}").pack()

        def current():
            return pngs[nb.index(nb.select())] or os.path.join(dest, "summary.txt")
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=8, pady=8)
        ttk.Button(bar, text="Open this figure", style="Ghost.TButton", command=lambda: open_folder(current())).pack(side="left")
        ttk.Button(bar, text="Open analysis folder (CSV)", style="Ghost.TButton", command=lambda: open_folder(dest)).pack(side="left", padx=6)
        ttk.Button(bar, text="Analyze again…", style="Ghost.TButton", command=lambda: (w.destroy(), self.analyze_now())).pack(side="left")
        ttk.Button(bar, text="Close", style="Ghost.TButton", command=w.destroy).pack(side="right")

    def open_plotter(self):
        """Plot: columns of the output files (thermo of log.lammps, fix ave/time files, csv) as curves of one figure."""
        if not self.out or not os.path.isdir(self.out):
            messagebox.showinfo("soapff-GUI", "First choose a structure file and run the calculation (or use \"Open existing output folder\")")
            return
        from . import plotter
        plotter.PlotDialog(self, self.out)

    def open_ovito(self):
        if not self.out or not os.path.isdir(self.out):
            messagebox.showinfo("soapff-GUI", "There is no output folder yet")
            return

        def go():
            try:
                p = runner.open_in_ovito(self.out)
                self.q.put(("log", f"Opened in OVITO: {p}\n", None))
            except Exception as e:
                self.q.put(("error", f"Could not open in OVITO: {e}"))
        threading.Thread(target=go, daemon=True).start()

    def open_out(self):
        if self.out and os.path.isdir(self.out):
            open_folder(self.out)

    def open_readme(self):
        p = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "README.md")
        if os.path.isfile(p):
            open_folder(p)

    def about(self):
        messagebox.showinfo("soapff-GUI", f"soapff-GUI v{__version__}\nBuilds a LAMMPS input with ClayFF (charges and types from SOAP) from a structure file, "
                                          "runs it on this computer or over SSH, and analyzes diffusion coefficients and PDF / RDF.")

    # ================================================================================================== queue
    def poll(self):
        n = 0
        try:
            while n < 500:
                msg = self.q.get_nowait()
                n += 1
                kind = msg[0]
                if kind == "log":
                    self.log_add(msg[1], msg[2])
                elif kind == "phase":
                    self.l_phase.config(text=msg[1])
                elif kind == "error":
                    messagebox.showerror("soapff-GUI", msg[1])
                elif kind == "call":
                    msg[1]()
                elif kind == "made":
                    self.set_state("Input created (up to minimization)", OK, 0)
                    p = os.path.join(msg[2], "in.lammps")
                    try:
                        basic = L.basic_part(open(p).read().splitlines())
                        with open(p, "w") as fh:
                            fh.write("\n".join(basic) + "\n")
                    except OSError:
                        pass
                    self.load_project(msg[1], msg[2])
                    self.sel = len(self.lines) - 1
                    self.render()
                    self.build_form()
                    self.log_add("The input ends at minimization. Add MD with \"+ fix\", \"+ compute\" and \"+ Other\" (run etc.).\n", None)
                    self.remember(STATUS_INPUT)
                    self.warn_assignment(msg[3])
                elif kind == "make_failed":
                    self.set_state("Could not build the force field", ERR, 0)
                    if "overlapping atoms" in msg[1]:
                        messagebox.showerror("soapff-GUI", "The structure has overlapping atoms, so the force field was not built.\n\n"
                                                           "Correct the structure (or minimise it first).\n\n" + msg[1])
                    elif "(fallback)" in msg[1]:
                        messagebox.showerror("soapff-GUI", "Atoms with environments missing from the library were found, so the force field was not built.\n\n"
                                                           "Add reference structures of that mineral / environment to the library.\n\n" + msg[1])
                    else:
                        messagebox.showerror("soapff-GUI", msg[1])
                elif kind == "thermo":
                    phase, v, prog = msg[1:]
                    self.l_phase.config(text=phase)
                    self.pb["value"] = prog
                    self.l_pct.config(text=f"{100 * prog:.0f}%")
                    keys = [("Step", "step {:.0f}"), ("Temp", "T {:.1f} K"), ("Press", "P {:.0f} atm"), ("PotEng", "PE {:.1f}")]
                    self.l_thermo.config(text="   ".join(f.format(v[k]) for k, f in keys if k in v))
                elif kind == "analysis":
                    self.show_analysis(msg[1])
                elif kind == "call":                      # a GUI update from a worker thread
                    msg[1]()
                elif kind == "idle":
                    self.set_running(False)
                elif kind == "done":
                    self.set_running(False)
                    ok = msg[2] == 0
                    self.set_state(msg[1] if len(msg[1]) < 40 else ("Done" if ok else "Finished (error)"), OK if ok else ERR,
                                   1.0 if ok else None)
                    self.log_add(f"\n== {msg[1]}\n", None if ok else "err")
                    self.remember(STATUS_DONE if ok else (STATUS_STOPPED if self.stopping else STATUS_ERROR))
                    if not ok and not self.stopping:
                        messagebox.showerror("soapff-GUI", msg[1])
        except queue.Empty:
            pass
        self.root.after(100, self.poll)

    def set_state(self, text, colour, progress):
        self.l_phase.config(text=text)
        self.dot.itemconfigure(self.dot_id, fill=colour)
        if progress is not None:
            self.pb["value"] = progress
            self.l_pct.config(text=f"{100 * progress:.0f}%")

    def set_running(self, on):
        txt = "Stop" if on else "Start calculation"
        self.b_start.config(text=txt)
        self.b_start_top.config(text=txt)
        if on:
            self.dot.itemconfigure(self.dot_id, fill=ACCENT)

    def log_add(self, text, tag):
        self.log.config(state="normal")
        if self.log.get("1.0", "1.end").startswith(LOG_PLACEHOLDER):
            self.log.delete("1.0", "end")
        self.log.insert("end", text, tag or ())
        if int(self.log.index("end-1c").split(".")[0]) > 20000:
            self.log.delete("1.0", "5000.0")
        self.log.see("end")
        self.log.config(state="disabled")

    def log_clear(self):
        self.log.config(state="normal")
        self.log.delete("1.0", "end")
        self.log.config(state="disabled")

    def on_close(self):
        if self.busy():
            if not messagebox.askyesno("soapff-GUI", "A run is in progress. Stop it and quit?\n(An SSH run continues on the remote side)"):
                return
            if self.job:
                self.job.stop()
        if self._save_after:
            self.save_input()
        save_settings(self.s)
        self.root.destroy()


def activate_macos(root):
    """bring the window to the front on macOS.  Started from soapff-GUI.app (the applet starts this Python in the
    background), the window would open behind the active application and look as if the app did not open:
    [NSApp activateIgnoringOtherApps:YES] through the Objective-C runtime (ctypes; no extra package)."""
    root.lift()
    try:
        import ctypes
        import ctypes.util
        objc = ctypes.cdll.LoadLibrary(ctypes.util.find_library("objc"))
        objc.objc_getClass.restype = objc.sel_registerName.restype = ctypes.c_void_p
        addr = ctypes.cast(objc.objc_msgSend, ctypes.c_void_p).value
        send0 = ctypes.CFUNCTYPE(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p)(addr)
        send1 = ctypes.CFUNCTYPE(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool)(addr)
        app = send0(objc.objc_getClass(b"NSApplication"), objc.sel_registerName(b"sharedApplication"))
        send1(app, objc.sel_registerName(b"activateIgnoringOtherApps:"), True)
    except Exception:
        pass


def main():
    if "--install-app" in sys.argv:                     # macOS: ~/Applications/soapff-GUI.app
        from .macapp import install
        k = sys.argv.index("--install-app")
        print("installed:", install(sys.argv[k + 1] if len(sys.argv) > k + 1 else None))
        return
    root = tk.Tk()
    root.tk.call("tk", "appname", "soapff-GUI")
    try:                                                 # window and (macOS) Dock icon
        root.iconphoto(True, tk.PhotoImage(file=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "icon.png")))
    except tk.TclError:
        pass
    App(root)
    if sys.platform == "darwin":
        root.after(300, lambda: activate_macos(root))
    root.mainloop()


if __name__ == "__main__":
    main()
