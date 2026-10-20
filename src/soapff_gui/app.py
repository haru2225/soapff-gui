"""soapff-GUI: choose a clay structure; soapff-cui writes the LAMMPS input (ClayFF + SOAP charges), which is shown line by
line and can be edited (select a line, change its values below); press 計算を開始 to run LAMMPS on this computer or on
another one over SSH, then analyse (diffusion coefficients; SOVA PDF, RDF, S(Q); the results of the computes;
XRD with Debyer) and look at the result in OVITO.

Layout (design/soapff-gui-ui-concept.png): menu row, toolbar (新規計算 / 計算設定 › 実行・解析 / 計算を開始), project sidebar
(新しい計算, 計算履歴, 構造ライブラリ, 最近使った計算), the input file (in.lammps) as an editor, the form of the selected
command, the run panel (status, progress, log) and a status bar."""
import json
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
STRUCTURE_TYPES = [("構造ファイル", "*.xyz *.extxyz *.cif *.vasp *.data *.lammps *.lmp POSCAR* CONTCAR*"), ("すべて", "*")]
AXIS_STYLES = ["p", "f", "m", "s", "fm", "mf", "fs", "sf", "ms", "sm"]
BARO = {"z": "z：面に垂直な方向だけ（層間・界面）", "aniso": "aniso：3 辺を別々に（バルク結晶）",
        "tri": "tri：箱の傾きも含めて全部（バルク結晶）", "iso": "iso：等方"}
DEFAULTS = dict(temp=300.0, relax_ps=10.0, nvt_ps=50.0, npt_ps=0.0, press=1.0, baro="aniso", boundary="p p p",
                nproc=max(1, min(3, (os.cpu_count() or 2) // 2)), threshold=0.95, cutoff=10.0, library="", lmp="", mpirun="",
                last_dir=os.path.expanduser("~"), dump_ps=1.0, auto_analyze=True, skip_ps="", sova_python="", sova_home="",
                analyses=["msd", "sova", "compute", "xrd", "props"], debyer="", props_selected=["elastic0", "tensile"],
                props_params={}, xrd_wavelength=1.5406, xrd_from=3.0, xrd_to=80.0,
                xrd_step=0.02, xrd_frames=3, xrd_mode="cluster", xrd_size=50.0,
                where="local", computers=[], history=[])
# the analyses (soapff_cui.analysis): key, label
ANALYSES = (("msd", "拡散係数（LAMMPS compute msd）"), ("sova", "SOVA：PDF g(r)・RDF N(r)・配位数・X 線 S(Q)"),
            ("compute", "compute の結果（fix ave/time・ave/chunk などのファイル）"),
            ("xrd", "XRD（debyer、粉末 X 線回折パターン）"),
            ("props", "物性（粘度・熱伝導率・弾性定数・引張・圧縮・せん断。＋ 物性 で入れたもの）"))
XRD_MODES = {"cluster": "結晶子：セルを並べた有限の結晶（結晶・粘土鉱物）",
             "periodic": "周期境界：カットオフ = 箱の半分（液体・ガラス、直方体の箱のみ）"}
# a registered computer (SSH): name, host, user, port, root (folder for the runs there), lmp, mpirun, setup (shell commands
# run first, e.g. 'module load lammps'), nproc.  The password is never stored: it is asked every time.
COMPUTER = dict(name="", host="", user="", port=22, root="soapff_runs", lmp="lmp", mpirun="mpirun", setup="", nproc=4)

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
    if not os.path.isfile(s.get("library") or ""):
        s["library"] = runner.default_library()
    lmp, mpirun = runner.find_executables()
    if not s.get("sova_python"):
        s["sova_python"], s["sova_home"] = runner.find_sova()
    if not s.get("debyer"):
        s["debyer"] = runner.find_debyer() or ""
    # analyses added in a later version are on: those not known when the settings were saved
    known = s.get("analyses_known") or ["msd", "sova", "compute", "xrd"]
    s["analyses"] = list(s.get("analyses") or []) + [k for k, _ in ANALYSES if k not in known and k not in (s.get("analyses") or [])]
    s["analyses_known"] = [k for k, _ in ANALYSES]
    s["lmp"] = s.get("lmp") or lmp or ""
    s["mpirun"] = s.get("mpirun") or mpirun or ""
    s["history"] = [h for h in s.get("history", []) if isinstance(h, dict)][:50]
    s["computers"] = [dict(COMPUTER, **c) for c in s.get("computers", []) if isinstance(c, dict) and c.get("name")]
    if s.get("r_host") and not s["computers"]:            # settings of an earlier version
        s["computers"].append(dict(COMPUTER, name=s["r_host"], host=s["r_host"], user=s.get("r_user", ""),
                                   port=s.get("r_port", 22)))
    for k in [k for k in s if k.startswith("r_")]:
        s.pop(k)
    if s.get("where") != "local" and not any(c["name"] == s.get("where") for c in s["computers"]):
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
        return "たった今"
    if d < 3600:
        return f"{int(d // 60)} 分前"
    if d < 86400:
        return f"{int(d // 3600)} 時間前"
    if d < 2 * 86400:
        return "昨日"
    return f"{int(d // 86400)} 日前"


def pick_font(cands, size, weight="normal"):
    fams = set(tkfont.families())
    fam = next((c for c in cands if c in fams), "TkDefaultFont")
    return (fam, size, weight)


FILES = ("in.lammps", "system.settings", "system.data")
FILE_NOTE = {"in.lammps": "LAMMPS input", "system.settings": "力場の設定（in.lammps から include）",
             "system.data": "原子・結合・角度（読み取り専用）"}
READONLY = {"system.data"}
# columns of the sections of a LAMMPS data file (atom_style full)
DATA_COLS = {"Masses": ["型番号", "質量"], "Atoms": ["原子番号", "分子番号", "型", "電荷 (e)", "x (Å)", "y (Å)", "z (Å)", "ix", "iy", "iz"],
             "Bonds": ["結合番号", "結合の型", "原子 1", "原子 2"], "Angles": ["角度番号", "角度の型", "原子 1", "原子 2 (中心)", "原子 3"],
             "Velocities": ["原子番号", "vx", "vy", "vz"]}


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
        root.title("soapff-GUI — 新しい計算")
        root.geometry("1360x860")
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
        ttk.Button(tb, text="新規計算", style="Ghost.TButton", command=self.new_calc).pack(side="left", padx=(18, 14), pady=10)
        tk.Frame(tb, bg=LINE, width=1).pack(side="left", fill="y", pady=12)
        self.crumb_settings = tk.Label(tb, text="計算設定", bg=BG, fg=MUTED, font=self.f_small, cursor="hand2")
        self.crumb_settings.pack(side="left", padx=(16, 8))
        self.crumb_settings.bind("<Button-1>", lambda e: self.open_settings())
        tk.Label(tb, text="›", bg=BG, fg=FAINT, font=self.f_small).pack(side="left")
        self.crumb_run = tk.Label(tb, text="実行・解析", bg=BG, fg=MUTED, font=self.f_small, cursor="hand2")
        self.crumb_run.pack(side="left", padx=8)
        self.crumb_run.bind("<Button-1>", lambda e: self.show_results())
        self.b_start_top = ttk.Button(tb, text="計算を開始", style="Accent.TButton", command=self.start_or_stop)
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
            ("ファイル", [("新規計算…", self.new_calc, "⌘N"), ("構造を変更…", self.new_calc),
                          ("既存の出力フォルダを開く…", self.open_project), None,
                          ("別のコンピュータからフォルダを受け取る…", self.fetch_dialog), None,
                          ("出力フォルダを開く", self.open_out), None, ("終了", self.on_close)]),
            ("編集", [("行を追加（自由入力）", self.new_line), ("minimize より後をすべて消す", self.cut_after_minimize), None,
                      ("物性の計算を追加…", self.open_properties),
                      ("fix を追加…", lambda: self.open_catalog("fix")), ("compute を追加…", lambda: self.open_catalog("compute")),
                      ("その他のコマンドを追加…", lambda: self.open_catalog("other")), None, ("行を削除", self.delete_line), None,
                      ("上へ移動", lambda: self.move_line(-1)), ("下へ移動", lambda: self.move_line(1)), None,
                      ("元に戻す", self.do_undo, "⌘Z")]),
            ("表示", [("in.lammps", lambda: self.show_tab("in.lammps")), ("system.settings（力場の設定）",
                                                                         lambda: self.show_tab("system.settings")),
                      ("system.data（原子・結合・角度）", lambda: self.show_tab("system.data")), None,
                      ("解析結果", self.show_results), ("OVITO で開く", self.open_ovito), None,
                      ("計算履歴", self.show_history), ("構造ライブラリ", self.show_library)]),
            ("計算", [("計算を開始", self.start, "⌘R"), ("停止", self.stop), None,
                      ("入力を作り直す（MD の条件から）", self.regenerate), ("解析", self.analyze_now), None,
                      ("計算設定…", self.open_settings)]),
            ("ヘルプ", [("README を開く", self.open_readme), ("soapff-GUI について", self.about)]),
        ]

    def _sidebar(self, body):
        sd = tk.Frame(body, bg=SIDE, width=228)
        sd.pack(side="left", fill="y")
        sd.pack_propagate(False)
        tk.Frame(body, bg=LINE, width=1).pack(side="left", fill="y")
        tk.Label(sd, text="プロジェクト", bg=SIDE, fg=MUTED, font=self.f_head).pack(anchor="w", padx=20, pady=(22, 8))
        self.nav = {}
        for key, label, cmd in (("new", "新しい計算", lambda: None), ("history", "計算履歴", self.show_history),
                                ("library", "構造ライブラリ", self.show_library)):
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
        tk.Label(sd, text="最近使った計算", bg=SIDE, fg=MUTED, font=self.f_head).pack(anchor="w", padx=20, pady=(0, 6))
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
        tk.Label(main, text="入力ファイル", bg=BG, fg="#23313b", font=self.f_title).pack(anchor="w")
        tk.Label(main, text="コマンドの行を選択して編集できます。", bg=BG, fg=MUTED, font=self.f_small).pack(anchor="w", pady=(0, 8))

        # structure card
        c = self._card(main, fill="x")
        left = tk.Frame(c, bg=PANEL)
        left.pack(side="left", padx=16, pady=6, fill="x", expand=True)
        tk.Label(left, text="構造", bg=PANEL, fg=FAINT, font=pick_font(["Hiragino Sans"], 10)).pack(anchor="w")
        self.l_struct = tk.Label(left, text="未選択（新規計算から構造ファイルを選んでください）", bg=PANEL, fg=TEXT, font=self.f_bold)
        self.l_struct.pack(anchor="w")
        self.l_struct2 = tk.Label(left, text="", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_struct2.pack(anchor="w")
        tk.Frame(c, bg=LINE, width=1).pack(side="left", fill="y", pady=8)
        mid = tk.Frame(c, bg=PANEL)
        mid.pack(side="left", padx=16, pady=6, fill="x", expand=True)
        tk.Label(mid, text="出力ファイル", bg=PANEL, fg=FAINT, font=pick_font(["Hiragino Sans"], 10)).pack(anchor="w")
        self.l_outfile = tk.Label(mid, text="in.lammps", bg=PANEL, fg=TEXT, font=self.f_mono)
        self.l_outfile.pack(anchor="w")
        self.l_outdir = tk.Label(mid, text="", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_outdir.pack(anchor="w")
        ttk.Button(c, text="構造を変更", style="Link.TButton", command=self.new_calc).pack(side="right", padx=16)

        # bottom row: selected command | run
        bottom = tk.Frame(main, bg=BG)
        bottom.pack(side="bottom", fill="x", pady=(10, 0))
        rc = tk.Frame(bottom, bg=PANEL, highlightthickness=1, highlightbackground=LINE, width=440, height=232)
        rc.pack(side="right", fill="y", padx=(14, 0))
        rc.pack_propagate(False)
        self.fc = tk.Frame(bottom, bg=PANEL, highlightthickness=1, highlightbackground=LINE)
        self.fc.pack(side="left", fill="both", expand=True)
        self._form_area()
        self._run_area(rc)

        # editor card
        ec = self._card(main, fill="both", expand=True, pady=(10, 0))
        head = tk.Frame(ec, bg=PANEL)
        head.pack(fill="x", padx=14, pady=(10, 0))
        self.tabs = {}
        for f in FILES:
            t = tk.Label(head, text=f" {f} ", bg=PANEL, fg="#23313b", font=pick_font(["Menlo"], 11, "bold"),
                         highlightthickness=1, highlightbackground="#c9d1d7", padx=8, pady=4, cursor="hand2")
            t.pack(side="left", padx=(0, 6))
            t.bind("<Button-1>", lambda e, f=f: self.show_tab(f))
            self.tabs[f] = t
        self.l_save = tk.Label(head, text="LAMMPS input · 自動保存", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_save.pack(side="right")
        # boundary of the box (system.settings); changing it rebuilds the force field (k-space solver, walls, centring)
        bf = tk.Frame(head, bg=PANEL)
        bf.pack(side="right", padx=(0, 18))
        tk.Label(bf, text="境界条件", bg=PANEL, fg=MUTED, font=self.f_small).pack(side="left", padx=(0, 6))
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
        foot = tk.Frame(ec, bg=PANEL)
        foot.pack(fill="x", padx=14, pady=6)
        for text, kind in (("＋ fix", "fix"), ("＋ compute", "compute"), ("＋ その他（run・timestep など）", "other")):
            ttk.Button(foot, text=text, style="Link.TButton", command=lambda k=kind: self.open_catalog(k)).pack(side="left", padx=(0, 6))
        ttk.Button(foot, text="＋ 行を追加（自由入力）", style="Link.TButton", command=self.new_line).pack(side="left", padx=(0, 6))
        ttk.Button(foot, text="＋ 物性", style="Small.Accent.TButton", command=self.open_properties).pack(side="left", padx=(0, 6))
        tk.Label(foot, text="行を選択して値を編集 · 上下キーで移動 · Enter で下に行を追加 · Delete で削除", bg=PANEL, fg=MUTED,
                 font=self.f_small).pack(side="right")
        self.b_insert_anchor = foot


    def _form_area(self):
        head = tk.Frame(self.fc, bg=PANEL)
        head.pack(fill="x", padx=16, pady=(10, 0))
        self.l_formtitle = tk.Label(head, text="選択したコマンド", bg=PANEL, fg="#23313b", font=self.f_bold)
        self.l_formtitle.pack(side="left")
        self.l_formline = tk.Label(head, text="", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_formline.pack(side="right")
        tk.Frame(self.fc, bg=LINE, height=1).pack(fill="x", pady=(8, 0))
        self.fields = tk.Frame(self.fc, bg=PANEL)
        self.fields.pack(fill="x", padx=16, pady=(10, 0))
        tk.Label(self.fc, text="生成される行", bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w", padx=16, pady=(10, 2))
        self.l_preview = tk.Label(self.fc, text="", bg=CODE_BG, fg=TEXT, font=self.f_mono_s, anchor="w", padx=12, pady=8)
        self.l_preview.pack(fill="x", padx=16, pady=(0, 14))

    def _run_area(self, rc):
        head = tk.Frame(rc, bg=PANEL)
        head.pack(fill="x", padx=16, pady=(8, 0))
        tk.Label(head, text="実行", bg=PANEL, fg="#23313b", font=self.f_bold).pack(side="left", pady=4)
        self.b_start = ttk.Button(head, text="計算を開始", style="Small.Accent.TButton", command=self.start_or_stop)
        self.b_start.pack(side="right")
        self.v_where = tk.StringVar()
        self.cb_where = ttk.Combobox(head, textvariable=self.v_where, state="readonly", width=18, font=self.f_small)
        self.cb_where.pack(side="right", padx=(0, 10))
        self.cb_where.bind("<<ComboboxSelected>>", lambda e: self.choose_where())
        tk.Frame(rc, bg=LINE, height=1).pack(fill="x", pady=(8, 0))
        row = tk.Frame(rc, bg=PANEL)
        row.pack(fill="x", padx=16, pady=(12, 0))
        self.dot = tk.Canvas(row, width=10, height=10, bg=PANEL, highlightthickness=0)
        self.dot.pack(side="left")
        self.dot_id = self.dot.create_oval(1, 1, 9, 9, fill=FAINT, outline="")
        self.l_phase = tk.Label(row, text="待機中", bg=PANEL, fg=TEXT, font=self.f_small)
        self.l_phase.pack(side="left", padx=6)
        self.l_pct = tk.Label(row, text="0%", bg=PANEL, fg=MUTED, font=self.f_small)
        self.l_pct.pack(side="right")
        self.pb = ttk.Progressbar(rc, maximum=1.0, style="Accent.Horizontal.TProgressbar")
        self.pb.pack(fill="x", padx=16, pady=(8, 0))
        self.l_thermo = tk.Label(rc, text="", bg=PANEL, fg=MUTED, font=self.f_mono_s, anchor="w")
        self.l_thermo.pack(fill="x", padx=16, pady=(6, 0))
        lf = tk.Frame(rc, bg=CODE_BG, highlightthickness=1, highlightbackground="#e1e6ea")
        lf.pack(fill="both", expand=True, padx=16, pady=(8, 8))
        self.log = tk.Text(lf, bg=CODE_BG, fg=TEXT, font=pick_font(["Menlo"], 10), relief="flat", height=4, wrap="none",
                           highlightthickness=0, padx=8, pady=6)
        self.log.pack(fill="both", expand=True)
        self.log.insert("end", "実行ログはここに表示されます。", "muted")
        self.log.tag_config("muted", foreground=FAINT)
        self.log.tag_config("warn", foreground="#a15c00")
        self.log.tag_config("err", foreground=ERR)
        self.log.config(state="disabled")
        acts = tk.Frame(rc, bg=PANEL)
        acts.pack(fill="x", padx=16, pady=(0, 12))
        for text, cmd in (("解析", self.analyze_now), ("OVITO で開く", self.open_ovito), ("フォルダ", self.open_out)):
            ttk.Button(acts, text=text, style="Link.TButton", command=cmd).pack(side="left", padx=(0, 6))
        tk.Label(acts, text="↑ 実行先", bg=PANEL, fg=FAINT, font=self.f_small).pack(side="right")

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
            ed.insert("end", "\n   構造ファイルを選ぶと、soapff-cui が LAMMPS の入力（in.lammps、system.settings）を作り、ここに表示します。\n", "comment")
            ed.insert("end", "   上の「新規計算」か「構造を変更」から始めてください。\n", "comment")
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
                ed.insert("end", "\t選択中 · クリックして編集", "hint")
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
            self.l_formtitle.config(text="選択したコマンド")
            self.l_formline.config(text="")
            self.l_preview.config(text="")
            return
        line = self.lines[self.sel]
        if self.tab in READONLY:
            return self.data_form(line)
        c, a, cm = L.split(line)
        self.l_formline.config(text=f"{self.sel + 1} 行目")
        if not c:
            self.l_formtitle.config(text="選択した行: コメント・空行")
            specs = [("内容", line.strip(), None)]
        else:
            self.l_formtitle.config(text=f"選択したコマンド: {c}")
            labs = L.labels(c, a)
            specs = [(lab, a[k] if k < len(a) else "", k) for k, lab in enumerate(labs)]
            if cm:
                specs.append(("コメント", cm, "comment"))
        self._form_cmd = c
        if c == "boundary":
            self.boundary_form(a)
            return
        if c == "include" and a and a[0] in FILES:
            tk.Button(self.fields, text=f"{a[0]} を開く", command=lambda f=a[0]: self.show_tab(f), relief="flat", bg=ACCENT_BG,
                      fg=ACCENT, activebackground=SEL_BG, font=self.f_small, padx=10, pady=4, cursor="hand2").pack(
                anchor="w", pady=(0, 6))
        col = 0
        rowf = None
        self.fc.update_idletasks()
        per_row = max(2, min(6, (self.fc.winfo_width() - 32) // 140)) if self.fc.winfo_width() > 50 else 6
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
        self.l_formline.config(text=f"{self.sel + 1} 行目（読み取り専用）")
        self.l_formtitle.config(text=f"system.data: {sec or 'ヘッダー'}")
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
            lab = cols[k] if k < len(cols) else ("値" if cols else f"{k + 1}")
            tk.Label(cell, text=lab, bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w")
            shown = v
            if sec == "Atoms" and k == 2 and v.isdigit() and int(v) in names:
                shown = f"{v}（{names[int(v)]}）"
            tk.Label(cell, text=shown, bg=CODE_BG, fg=TEXT, font=self.f_mono_s, width=14, anchor="w", padx=6, pady=4).pack(anchor="w")
        if comment.strip():
            tk.Label(self.fields, text=f"コメント：{comment.strip()}", bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w")
        tk.Label(self.fields, text="system.data は soapff-cui が書く原子・結合・角度の表です（ここでは見るだけ）。"
                                   "変えるときは構造ファイルか境界条件を変えて作り直します。",
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
        ttk.Button(row, text="この境界条件で作り直す", style="Small.Accent.TButton",
                   command=lambda: self.change_boundary([v.get() for v in vs])).pack(side="left", padx=(8, 0), pady=(14, 0))
        tk.Label(self.fields, text="p：周期　f：固定　s：原子に合わせて縮む　m：原子に合わせて広がる（面ごとに fm など）。"
                                   "変えると力場（長距離クーロン力の計算法など）と入力を作り直します。",
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
            messagebox.showinfo("soapff-GUI", "先に構造ファイルを選んで、入力を作ってください")
            return
        if self.tab != "in.lammps":
            self.show_tab("in.lammps")
        items = catalog.entries(kind) + (catalog.entries("special") if kind in ("other", "compute") else [])
        w = tk.Toplevel(self.root)
        w.title({"fix": "fix を追加", "compute": "compute を追加", "other": "コマンドを追加"}[kind])
        w.geometry("980x600")
        w.configure(bg=PANEL)
        left = tk.Frame(w, bg=PANEL)
        left.pack(side="left", fill="both", padx=(12, 6), pady=12)
        q = tk.StringVar()
        ttk.Entry(left, textvariable=q, style="Field.TEntry", width=40).pack(fill="x")
        tk.Label(left, text="名前・説明で絞り込み", bg=PANEL, fg=FAINT, font=self.f_small).pack(anchor="w", pady=(2, 6))
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
        tk.Label(right, text="入る行", bg=PANEL, fg=MUTED, font=self.f_small).pack(anchor="w", pady=(12, 2))
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
            l_fields.config(text=("欄：" + "、".join(labs)) if labs else "")
            n = self.sel + 1 if self.sel is not None else len(self.lines)
            where.config(text=f"{n} 行目のあとに入ります（入れたあと、下の欄で値を直せます）")

        def insert(*_):
            e = chosen()
            if e:
                self.insert_lines(self.catalog_lines(e))
                w.destroy()

        def manual():
            if cur.get("doc"):
                webbrowser.open(cur["doc"])
        ttk.Button(bar, text="挿入", style="Accent.TButton", command=insert).pack(side="right")
        ttk.Button(bar, text="LAMMPS のマニュアル", style="Ghost.TButton", command=manual).pack(side="right", padx=8)
        ttk.Button(bar, text="閉じる", style="Ghost.TButton", command=w.destroy).pack(side="left")
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
            messagebox.showinfo("soapff-GUI", "minimize より後のコマンドはありません")
            return
        if not messagebox.askyesno("soapff-GUI", f"in.lammps の minimize より後のコマンド {len(extra)} 行"
                                                 f"（{extra[0]} …）を消します。よろしいですか？（元に戻す で戻せます）"):
            return
        if self.tab != "in.lammps":
            self.show_tab("in.lammps")
        self.push_undo()
        self.lines = L.basic_part(lines)
        self.sel = len(self.lines) - 1
        self.render()
        self.build_form()
        self.schedule_save()

    # ================================================================================================== properties
    PROP_FIELDS = (("共通", (("temp", "温度 (K)", ""), ("press", "圧力 (atm)", ""), ("relax_ps", "緩和 NVT (ps)", ""),
                             ("npt_ps", "NPT 平衡化 (ps)", "0 でなし。終わると箱を NPT の平均に"),
                             ("prod_ps", "有限温度の弾性定数の MD (ps)", "数百 ps〜ns 要る"))),
                   ("弾性定数", (("e0_strain", "0 K：ひずみ", "± この大きさで 6 方向"),
                                 ("eT_stress_every", "有限温度：応力の記録間隔（ステップ）", ""),
                                 ("eT_born_every", "有限温度：Born 項の間隔（ステップ）", "1 回で力の計算 12 回分"))),
                   ("引張・圧縮", (("tensile_rate", "ひずみ速度 (1/s)", ""), ("tensile_max", "最大ひずみ", ""))),
                   ("せん断", (("shear_rate", "ひずみ速度 (1/s)", ""), ("shear_max", "最大ひずみ", ""))))

    def open_properties(self):
        """choose the properties to compute: their LAMMPS commands (soapff_cui.properties) are added to in.lammps."""
        from soapff_cui import properties as PR
        if not self.docs["in.lammps"]:
            messagebox.showinfo("soapff-GUI", "先に構造ファイルを選んで、入力を作ってください")
            return
        lines = self.docs["in.lammps"]
        dt = L.timestep(lines, default=float(getattr(runner.import_soapff(), "TIMESTEP_FS", 0.5)))
        p = dict(PR.DEFAULTS, temp=self.s["temp"], press=self.s["press"], relax_ps=self.s["relax_ps"],
                 baro=self.s.get("baro", "aniso"), dump_ps=self.s["dump_ps"], timestep=dt)
        p.update(self.s.get("props_params") or {})
        if float(p.get("e0_strain", PR.DEFAULTS["e0_strain"])) < PR.DEFAULTS["e0_strain"]:
            p["e0_strain"] = PR.DEFAULTS["e0_strain"]        # 0.001 of an earlier version: too small with PPPM (noise)
        p["timestep"] = dt
        w = tk.Toplevel(self.root)
        w.title("物性の計算を追加")
        w.configure(bg=PANEL)
        w.transient(self.root)
        body = ttk.Frame(w, padding=14)
        body.pack(fill="both", expand=True)
        left, right = ttk.Frame(body), ttk.Frame(body)
        left.pack(side="left", fill="y", padx=(0, 18))
        right.pack(side="left", fill="both", expand=True)
        ttk.Label(left, text="知りたい物性", style="Section.TLabel").pack(anchor="w")
        sel = {}
        for k, lab in PR.PROPS:
            sel[k] = tk.BooleanVar(value=k in self.s.get("props_selected", []))
            ttk.Checkbutton(left, text=lab, variable=sel[k]).pack(anchor="w", pady=2)
        note = ("in.lammps の最後に、0 K の弾性定数 → 平衡化（NVT・NPT）→ 有限温度の弾性定数の MD → 引張・圧縮、せん断\n"
                "（それぞれ平衡化した状態から）の順に入ります。周期境界（p p p）の構造ならどの鉱物でも使えます。\n"
                "NPT の方向は 計算設定 の「圧力の方向」（バルクは aniso か tri、層状で層間が広がるなら z）です。")
        ttk.Label(left, text=note.replace("\n", ""), foreground=MUTED, justify="left", font=self.f_small,
                  wraplength=380).pack(anchor="w", pady=(10, 0))
        v = {}
        nbp = ttk.Notebook(right)
        nbp.pack(fill="both", expand=True)
        for sec, fields in self.PROP_FIELDS:
            tab = ttk.Frame(nbp, padding=12)
            nbp.add(tab, text=sec.split("（")[0])
            row = 0
            for k, lab, nt in fields:
                ttk.Label(tab, text=lab).grid(row=row, column=0, sticky="w", padx=(8, 10), pady=2)
                v[k] = tk.StringVar(value=f"{p[k]:g}" if isinstance(p[k], float) else str(p[k]))
                ttk.Entry(tab, textvariable=v[k], width=10, style="Field.TEntry").grid(row=row, column=1, sticky="w", pady=2)
                if nt:
                    ttk.Label(tab, text=nt, foreground=MUTED, wraplength=300, justify="left").grid(row=row, column=2,
                                                                                                    sticky="w", padx=6)
                row += 1
            if sec == "共通":
                ttk.Label(tab, text="NPT の方向").grid(row=row, column=0, sticky="w", padx=(8, 10), pady=2)
                v["baro"] = tk.StringVar(value=p["baro"])
                ttk.Combobox(tab, textvariable=v["baro"], values=list(BARO), width=8, state="readonly").grid(
                    row=row, column=1, sticky="w")
                ttk.Label(tab, text=BARO.get(p["baro"], ""), foreground=MUTED).grid(row=row, column=2, sticky="w", padx=6)
                row += 1
            if sec in ("引張・圧縮", "せん断"):
                f = ttk.Frame(tab)
                f.grid(row=row, column=0, columnspan=3, sticky="w", padx=(8, 0), pady=2)
                key, axes = ("tensile_dirs", ("x", "y", "z")) if sec == "引張・圧縮" else ("shear_planes", ("xy", "xz", "yz"))
                ttk.Label(f, text="方向" if sec == "引張・圧縮" else "面").pack(side="left", padx=(0, 8))
                chosen = str(p[key]).replace(",", " ").split()
                for a in axes:
                    v[f"{key}:{a}"] = tk.BooleanVar(value=a in chosen)
                    ttk.Checkbutton(f, text=a, variable=v[f"{key}:{a}"]).pack(side="left")
                mkey, modes = (("tensile_mode", (("tension", "引張"), ("compression", "圧縮"), ("both", "両方")))
                               if sec == "引張・圧縮" else ("shear_mode", (("solid", "固体（fix nvt）"),
                                                                       ("flow", "流れ（fix nvt/sllod → 粘度）"))))
                v[mkey] = tk.StringVar(value=p[mkey])
                ttk.Label(f, text="　").pack(side="left")
                for val, lab in modes:
                    ttk.Radiobutton(f, text=lab, value=val, variable=v[mkey]).pack(side="left")
                row += 1
        est = tk.StringVar()
        ttk.Label(w, textvariable=est, foreground=MUTED, padding=(14, 0)).pack(anchor="w")

        def collect():
            q = dict(p)
            for _, fields in self.PROP_FIELDS:
                for k, lab, _ in fields:
                    try:
                        q[k] = float(v[k].get())
                    except ValueError:
                        raise ValueError(f"{lab} が数値ではありません")
            for k in ("eT_stress_every", "eT_born_every"):
                q[k] = max(1, int(q[k]))
            q["baro"] = v["baro"].get()
            q["tensile_dirs"] = " ".join(a for a in "xyz" if v[f"tensile_dirs:{a}"].get()) or "x"
            q["shear_planes"] = " ".join(a for a in ("xy", "xz", "yz") if v[f"shear_planes:{a}"].get()) or "xy"
            q["tensile_mode"], q["shear_mode"] = v["tensile_mode"].get(), v["shear_mode"].get()
            props = [k for k, _ in PR.PROPS if sel[k].get()]
            return props, q

        def make(props, q):
            tid = self.type_ids()
            existing = set(L.ids(self.docs["in.lammps"], "fix")) | set(L.ids(self.docs["in.lammps"], "compute"))
            return PR.lines(props, q, tid, runner.import_soapff().analysis_lines, existing)

        def update(*_):
            try:
                props, q = collect()
            except ValueError as e:
                est.set(str(e))
                return
            if not props:
                est.set("物性を選んでください")
                return
            new = make(props, q)
            n = sum(L.run_steps(new))
            nmin = sum(1 for l in new if l.startswith("minimize"))
            est.set(f"追加される行 {len(new)}、MD 合計 {n:,} ステップ（{n * q['timestep'] / 1000:g} ps）"
                    + (f"、最小化 {nmin} 回" if nmin else ""))
        for var in list(v.values()) + list(sel.values()):
            var.trace_add("write", update)
        update()

        def add():
            try:
                props, q = collect()
            except ValueError as e:
                messagebox.showerror("soapff-GUI", str(e), parent=w)
                return
            if not props:
                messagebox.showinfo("soapff-GUI", "物性を選んでください", parent=w)
                return
            bnd = self.s.get("boundary", "p p p")
            if bnd.split() != ["p", "p", "p"] and any(k in props for k in ("elastic0", "elasticT", "tensile", "shear")):
                messagebox.showerror("soapff-GUI", f"境界条件が {bnd} です。弾性定数・引張・せん断は周期境界（p p p）でだけ計算できます",
                                     parent=w)
                return
            if any(l.startswith("# soapff-property:") for l in self.docs["in.lammps"]) and not messagebox.askyesno(
                    "soapff-GUI", "in.lammps にはもう物性の計算があります。後ろにもう 1 つ足しますか？"
                                  "（出力ファイルは後のもので上書きされます）", parent=w):
                return
            new = make(props, q)
            if self.tab != "in.lammps":
                self.show_tab("in.lammps")
            self.push_undo()
            cur = list(self.lines)
            if PR.needs_triclinic(props) and not any(l.split()[:3] == ["change_box", "all", "triclinic"] for l in cur):
                k = max((i for i, l in enumerate(cur) if L.split(l)[0] == "minimize"), default=len(cur) - 1)
                cur[k + 1:k + 1] = PR.triclinic_lines(self.out)
            self.lines = cur + new
            self.sel = len(self.lines) - 1
            self.render()
            self.build_form()
            self.schedule_save()
            self.s["props_selected"], self.s["props_params"] = props, {k: q[k] for k in q if k != "timestep"}
            save_settings(self.s)
            self.log_add(f"物性の計算を in.lammps に追加しました：{', '.join(props)}（{len(new)} 行）。"
                         "計算が終わったら「解析」で結果の図が出ます\n", None)
            w.destroy()
        bar = ttk.Frame(w, padding=(14, 8, 14, 14))
        bar.pack(fill="x")
        ttk.Button(bar, text="in.lammps に追加", style="Accent.TButton", command=add).pack(side="right")
        ttk.Button(bar, text="閉じる", style="Ghost.TButton", command=w.destroy).pack(side="left")

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
        lines = [f"ライブラリ：{rep.get('library', '?')}"]
        if fb:
            per = rep.get("per_type") or {}
            kinds = "、".join(f"{t} {v['fallback']}" for t, v in per.items() if v.get("fallback"))
            lines.append(f"・ライブラリで決まらず、結合のつながりから型を決めた原子：{fb} 個（{kinds}）")
        if nd:
            kinds = "、".join(f"{k} が {v} 個" for k, v in (tc.get("kinds") or {}).items())
            lines.append(f"・SOAP の型が結合のつながりの型と違う原子：{nd} 個（{kinds}。「結合の型→SOAP の型」）")
        lines.append("この鉱物の環境がライブラリにない可能性があります。計算設定 › ライブラリと詳細 で同梱のライブラリ"
                     f"（{runner.default_library()}）を選ぶか、この鉱物の参照をライブラリに加えてください。"
                     "原子ごとの値は charges.csv（列 source、bond_type）にあります。")
        self.log_add("警告：" + " ".join(lines) + "\n", "warn")
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
                return msd or ["# 水もイオンもないので、拡散係数の MSD はありません"]
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

    def insert_lines(self, new):
        if self.tab in READONLY:
            return
        self.push_undo()
        i = (self.sel + 1) if self.sel is not None else len(self.lines)
        self.lines[i:i] = list(new)
        self.select(i)
        self.schedule_save()

    def new_line(self):
        """an empty line below the selected one, typed in the 内容 field of the form (e.g. 'run 10000')."""
        if self.tab in READONLY:
            return
        if not self.lines and not self.out:
            messagebox.showinfo("soapff-GUI", "先に構造ファイルを選んで、入力を作ってください")
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
        self.l_save.config(text=f"{FILE_NOTE[self.tab]} · 保存待ち")

    def save_input(self):
        self._save_after = None
        if self.out and os.path.isdir(self.out):
            for f in sorted(self.dirty):
                if self.docs[f]:
                    with open(os.path.join(self.out, f), "w") as fh:
                        fh.write("\n".join(self.docs[f]) + "\n")
            self.dirty.clear()
            self.l_save.config(text=f"{FILE_NOTE[self.tab]} · 自動保存（{time.strftime('%H:%M:%S')}）")

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
        self.l_save.config(text=f"{FILE_NOTE[f]} · 自動保存")
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
            messagebox.showinfo("soapff-GUI", "計算中は変えられません")
            self.set_boundary_widgets(old)
            return
        if new != ["p", "p", "p"] and float(self.s.get("npt_ps", 0)) > 0:
            messagebox.showerror("soapff-GUI", "周期でない方向があるときは NPT を使えません（真空を含む箱の圧力は決まらないため）。\n"
                                               "計算設定の「MD の条件」で NPT を 0 にしてから選んでください")
            self.set_boundary_widgets(old)
            return
        if not self.structure or not os.path.isfile(self.structure):
            messagebox.showerror("soapff-GUI", "構造ファイルがわからないので作り直せません。「新規計算」から選んでください")
            self.set_boundary_widgets(old)
            return
        if not messagebox.askyesno("soapff-GUI", f"境界条件を {' '.join(old)} から {' '.join(new)} に変えて、力場と入力を作り直します。\n"
                                                 "（長距離クーロン力の計算法や、周期でない方向の扱いも変わります）\n"
                                                 "in.lammps と system.settings を手で直した行は消えます。よろしいですか？"):
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
            messagebox.showinfo("soapff-GUI", "計算中です。終わるか停止してから選んでください")
            return
        p = filedialog.askopenfilename(title="構造ファイルを選ぶ", initialdir=self.s["last_dir"], filetypes=STRUCTURE_TYPES)
        if not p:
            return
        self.s["last_dir"] = os.path.dirname(p)
        stem = os.path.splitext(os.path.basename(p))[0]
        out = os.path.join(os.path.dirname(p), f"{stem}_soapff")
        if os.path.isfile(os.path.join(out, "in.lammps")):
            ans = messagebox.askyesnocancel("soapff-GUI", f"{out} に前の入力があります。\n\nはい：前の入力をそのまま開く\n"
                                                          f"いいえ：力場と入力を作り直す（上書き）")
            if ans is None:
                return
            if ans:
                self.load_project(p, out)
                return
        self.make_input(p, out)

    def open_project(self):
        d = filedialog.askdirectory(title="出力フォルダ（in.lammps のあるフォルダ）を選ぶ", initialdir=self.s["last_dir"])
        if d and os.path.isfile(os.path.join(d, "in.lammps")):
            st = next((h["structure"] for h in self.s["history"] if h.get("out") == d), "")
            self.load_project(st, d)
        elif d:
            messagebox.showerror("soapff-GUI", "このフォルダには in.lammps がありません")

    def make_input(self, structure, out):
        """soapff make in a worker thread (writes system.data, system.settings, in.lammps ...)."""
        md = {k: self.s[k] for k in ("temp", "relax_ps", "nvt_ps", "npt_ps", "press", "baro", "boundary", "threshold",
                                     "cutoff", "dump_ps")}
        if not os.path.isfile(self.s["library"]):
            messagebox.showerror("soapff-GUI", f"SOAP ライブラリがありません：{self.s['library']}\n計算設定で選んでください")
            return
        self.log_clear()
        self.set_state("準備中", FAINT, 0)
        self.l_phase.config(text="力場を作成中…")

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
            messagebox.showerror("soapff-GUI", f"入力を読めません：{e}")
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
        self.remember(status=self.history_entry().get("status", "入力作成"))

    def show_project(self):
        name = os.path.basename(self.structure) if self.structure else os.path.basename(self.out)
        r = self.report
        atoms = r["atoms"] if r else None
        self.l_struct.config(text=f"{name}  ·  {atoms:,} atoms" if atoms else name)
        if r:
            ms = min((p["min_similarity"] for p in r["per_type"].values()), default=float("nan"))
            roles = r.get("roles", {})
            self.l_struct2.config(text=f"粘土 {roles.get('mineral', 0)} · 水 {roles.get('water', 0)} · イオン {roles.get('ion', 0)}"
                                       f"   SOAP {r['from_soap']} / fallback {r['fallback']} · 最小類似度 {ms:.3f}",
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
            col = OK if h.get("status") == "完了" else (ACCENT if h.get("status") == "計算中" else
                                                        (ERR if h.get("status") in ("エラー",) else FAINT))
            tk.Label(row, text="●", bg=SIDE, fg=col, font=pick_font(["Hiragino Sans"], 8)).pack(side="left", anchor="n",
                                                                                                padx=(6, 6), pady=3)
            box = tk.Frame(row, bg=SIDE)
            box.pack(side="left", fill="x")
            a = tk.Label(box, text=h.get("name", "?"), bg=SIDE, fg=TEXT, font=self.f_small, anchor="w")
            a.pack(anchor="w")
            b = tk.Label(box, text=f"{h.get('status', '')} · {ago(h.get('time'))}", bg=SIDE, fg=FAINT,
                         font=pick_font(["Hiragino Sans"], 10), anchor="w")
            b.pack(anchor="w")
            for w in (row, box, a, b):
                w.bind("<Button-1>", lambda e, h=h: self.open_history(h))

    def open_history(self, h):
        if self.busy():
            messagebox.showinfo("soapff-GUI", "計算中です")
            return
        if os.path.isfile(os.path.join(h.get("out", ""), "in.lammps")):
            self.load_project(h.get("structure", ""), h["out"])
        else:
            messagebox.showerror("soapff-GUI", f"見つかりません：{h.get('out')}")

    def show_history(self):
        self._nav_select("history")
        w = tk.Toplevel(self.root)
        w.title("計算履歴")
        w.geometry("900x420")
        tv = ttk.Treeview(w, columns=("status", "when", "where", "out"), show="tree headings")
        for col, txt, wd in (("#0", "名前", 200), ("status", "状態", 90), ("when", "日時", 150), ("where", "実行先", 180),
                             ("out", "出力フォルダ", 400)):
            tv.heading(col, text=txt)
            tv.column(col, width=wd, anchor="w")
        for k, h in enumerate(self.s["history"]):
            tv.insert("", "end", iid=str(k), text=h.get("name", ""),
                      values=(h.get("status", ""), time.strftime("%Y-%m-%d %H:%M", time.localtime(h.get("time", 0))),
                              h.get("where", ""), h.get("out", "")))
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
        ttk.Button(bar, text="開く", style="Accent.TButton", command=open_).pack(side="left")
        ttk.Button(bar, text="フォルダを開く", style="Ghost.TButton", command=folder).pack(side="left", padx=6)
        ttk.Button(bar, text="履歴から外す", style="Ghost.TButton", command=forget).pack(side="left")
        tv.bind("<Double-1>", lambda e: open_())
        w.bind("<Destroy>", lambda e: self._nav_select("new") if e.widget is w else None)

    def show_library(self):
        self._nav_select("library")
        w = tk.Toplevel(self.root)
        w.title("構造ライブラリ")
        w.geometry("760x520")
        txt = scrolledtext.ScrolledText(w, font=self.f_mono_s, height=20)
        txt.pack(fill="both", expand=True, padx=10, pady=10)

        def fill():
            txt.config(state="normal")
            txt.delete("1.0", "end")
            p = self.s["library"]
            txt.insert("end", f"ライブラリ：{p}\n\n")
            try:
                import numpy as np
                import collections
                d = np.load(p, allow_pickle=False)
                el, ty, q = d["element"], d["type"], d["q"]
                org = d["origin"] if "origin" in d.files else []
                txt.insert("end", f"{len(el)} 行（1 行 = 参照の構造の原子 1 個：SOAP ベクトル {d['X'].shape[1]} 成分、電荷、型）\n")
                txt.insert("end", f"元素：{' '.join(d['species'])}\n\n参照の構造\n")
                for o, n in collections.Counter(org).items():
                    txt.insert("end", f"  {o:30s} {n:6d} 原子\n")
                txt.insert("end", "\n型ごとの行数と電荷\n")
                for (e, t), n in sorted(collections.Counter(zip(el, ty)).items()):
                    qs = sorted({round(float(x), 4) for x, ee, tt in zip(q, el, ty) if ee == e and tt == t})
                    txt.insert("end", f"  {e:3s} {t:6s} {n:6d}   電荷 {qs}\n")
            except Exception as e:
                txt.insert("end", f"読めません：{e}\n")
            txt.config(state="disabled")

        def choose():
            p = filedialog.askopenfilename(title="ライブラリを選ぶ", filetypes=[("SOAP ライブラリ", "*.npz")])
            if p:
                self.s["library"] = p
                save_settings(self.s)
                fill()
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=10, pady=(0, 10))
        ttk.Button(bar, text="別のライブラリを選ぶ…", style="Ghost.TButton", command=choose).pack(side="left")
        ttk.Button(bar, text="閉じる", style="Ghost.TButton", command=w.destroy).pack(side="right")
        fill()
        w.bind("<Destroy>", lambda e: self._nav_select("new") if e.widget is w else None)

    # ================================================================================================== settings
    def open_settings(self):
        w = tk.Toplevel(self.root)
        w.title("計算設定")
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
                ttk.Button(parent, text="選ぶ…", style="Ghost.TButton",
                           command=lambda: v[key].set((filedialog.askdirectory() if browse == "dir" else
                                                       filedialog.askopenfilename()) or v[key].get())).grid(row=row, column=2, padx=6)
            if note:
                ttk.Label(parent, text=note, foreground=MUTED).grid(row=row, column=3 if browse else 2, sticky="w", padx=6)
            return e

        # --- MD
        md = ttk.Frame(nb, padding=16)
        nb.add(md, text="MD の条件")
        for i, (k, lab, note) in enumerate((("temp", "温度 (K)", ""), ("relax_ps", "緩和 NVT (ps)", "最小化のあと。記録しない"),
                                            ("nvt_ps", "NVT (ps)", ""), ("npt_ps", "NPT (ps)", "0 でなし"),
                                            ("press", "圧力 (atm)", ""), ("dump_ps", "出力間隔 (ps)", "拡散係数なら 0.1〜0.5"))):
            entry(md, i, k, lab, note=note)
        ttk.Label(md, text="圧力の方向").grid(row=6, column=0, sticky="w", pady=4)
        v["baro"] = tk.StringVar(value=self.s["baro"])
        ttk.Combobox(md, textvariable=v["baro"], values=list(BARO), width=8, state="readonly").grid(row=6, column=1, sticky="w")
        ttk.Label(md, text="z：層間・界面、aniso / tri：バルク結晶", foreground=MUTED).grid(row=6, column=2, sticky="w", padx=6)
        ttk.Label(md, text="境界条件 (x y z)").grid(row=7, column=0, sticky="w", pady=4)
        bf = ttk.Frame(md)
        bf.grid(row=7, column=1, columnspan=2, sticky="w")
        b0 = (self.s["boundary"].split() + ["p", "p", "p"])[:3]
        vb = []
        for k in range(3):
            x = tk.StringVar(value=b0[k])
            ttk.Combobox(bf, textvariable=x, values=AXIS_STYLES, width=3, state="readonly").pack(side="left", padx=(0, 4))
            vb.append(x)
        ttk.Label(md, text="この条件は、次に入力を作るとき（新規計算、または下のボタン）に使います。\n"
                           "作ったあとの入力は、メイン画面で行ごとに直せます。", foreground=MUTED, justify="left").grid(
            row=8, column=0, columnspan=4, sticky="w", pady=(14, 6))

        # --- where to run
        ex = ttk.Frame(nb, padding=16)
        nb.add(ex, text="実行先")
        ex.columnconfigure(0, weight=1)
        lo = ttk.LabelFrame(ex, text="このコンピュータ", padding=10)
        lo.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        lo.columnconfigure(1, weight=1)
        entry(lo, 0, "lmp", "LAMMPS (lmp)", width=46, browse="file")
        entry(lo, 1, "mpirun", "mpirun", width=46, browse="file")
        entry(lo, 2, "nproc", "並列数", width=6, note="ノート PC では P コアの数以下に")
        rm = ttk.LabelFrame(ex, text="登録したコンピュータ（SSH）", padding=10)
        rm.grid(row=1, column=0, sticky="nsew")
        ex.rowconfigure(1, weight=1)
        rm.columnconfigure(0, weight=1)
        comps = [dict(c) for c in self.s["computers"]]
        tv = ttk.Treeview(rm, columns=("host", "user", "port", "nproc"), show="tree headings", height=5)
        for col, txt, wd in (("#0", "名前", 170), ("host", "ホスト", 170), ("user", "ユーザー", 150), ("port", "ポート", 70),
                             ("nproc", "並列数", 70)):
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
            if k is not None and messagebox.askyesno("soapff-GUI", f"「{comps[k]['name']}」の登録を消しますか？", parent=w):
                del comps[k]
                refill()
        res = tk.StringVar(value="パスワードは保存しません。接続テストと計算を始めるたびに聞きます。")

        def test():
            k = chosen()
            if k is None:
                res.set("テストするコンピュータを選んでください")
                return
            c = comps[k]
            pw = self.ask_password(c, parent=w)
            if pw is None:
                return
            res.set(f"{c['name']} に接続しています…")

            def go():
                r = None
                try:
                    r = self.make_remote(c, pw)
                    info = r.check()
                    vers = "　".join(f"{k[15:-1]}：{v}" for k, v in info.items() if k.startswith("mpirun_version["))
                    msg = (f"接続できました：{c['name']} = {info.get('host')}（{info.get('os')}、CPU {info.get('cpus')}、"
                           f"物理コア {info.get('cores') or '?'}。並列数は物理コア数以下に）\n"
                           f"lmp：{info.get('lmp')}（{info.get('lmp_mpi') or 'MPI 不明'}）\n"
                           f"使う mpirun：{info.get('mpirun_used')}　{vers}\nフォルダ：{info.get('root')}")
                    if info.get("mpi_warning"):
                        msg += "\n⚠ " + info["mpi_warning"]
                    if "NOTFOUND" in (info.get("lmp", "") + info.get("mpirun", "")):
                        msg += "\n⚠ lmp か mpirun が見つかりません。登録の「相手側の lmp / mpirun」か「最初に実行するコマンド」を設定してください"
                except Exception as e:
                    msg = f"接続できませんでした：{e}"
                finally:
                    if r:
                        r.close()
                self.q.put(("call", lambda: res.set(msg)))
            threading.Thread(target=go, daemon=True).start()
        btns = ttk.Frame(rm)
        btns.grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Button(btns, text="追加…", style="Ghost.TButton", command=edit).pack(side="left")
        ttk.Button(btns, text="編集…", style="Ghost.TButton", command=lambda: chosen() is not None and edit(chosen())).pack(
            side="left", padx=6)
        ttk.Button(btns, text="削除", style="Ghost.TButton", command=remove).pack(side="left")
        ttk.Button(btns, text="接続テスト", style="Ghost.TButton", command=test).pack(side="left", padx=6)
        ttk.Label(rm, textvariable=res, foreground=MUTED, wraplength=660, justify="left").grid(row=2, column=0, sticky="w",
                                                                                              pady=(8, 0))
        tv.bind("<Double-1>", lambda e: chosen() is not None and edit(chosen()))
        refill()

        # --- library and details
        de = ttk.Frame(nb, padding=16)
        nb.add(de, text="ライブラリと詳細")
        de.columnconfigure(1, weight=1)
        entry(de, 0, "library", "SOAP ライブラリ", width=50, browse="file")
        ttk.Button(de, text="同梱のものに戻す", style="Ghost.TButton",
                   command=lambda: v["library"].set(runner.default_library())).grid(row=0, column=3, padx=(0, 6))
        entry(de, 1, "threshold", "類似度のしきい値", width=8)
        entry(de, 2, "cutoff", "カットオフ (Å)", width=8)

        # --- analysis
        an = ttk.Frame(nb, padding=16)
        nb.add(an, text="解析")
        self.analysis_form(an, v, self.s, auto=True)

        def collect():
            s2 = dict(self.s)
            try:
                for k in ("temp", "relax_ps", "nvt_ps", "npt_ps", "press", "dump_ps", "threshold", "cutoff"):
                    s2[k] = float(v[k].get())
                s2["nproc"] = int(float(v["nproc"].get() or 1))
            except ValueError as e:
                raise ValueError(f"数値ではありません：{e}")
            if s2["dump_ps"] <= 0 or s2["temp"] <= 0:
                raise ValueError("温度と出力間隔は正の値にしてください")
            for k in ("baro", "lmp", "mpirun", "library"):
                s2[k] = v[k].get().strip()
            self.collect_analysis(v, s2)
            s2["boundary"] = " ".join(x.get() for x in vb)
            s2["computers"] = [dict(c) for c in comps]
            if s2["where"] != "local" and not any(c["name"] == s2["where"] for c in comps):
                s2["where"] = "local"
            if s2["npt_ps"] > 0 and s2["boundary"] != "p p p":
                raise ValueError("周期でない方向があるときは NPT を使えません")
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
        ttk.Button(bar, text="保存", style="Accent.TButton", command=ok).pack(side="right")
        ttk.Button(bar, text="保存して入力を作り直す", style="Ghost.TButton", command=lambda: ok(True)).pack(side="right", padx=8)
        ttk.Button(bar, text="閉じる", style="Ghost.TButton", command=w.destroy).pack(side="left")
        w.update_idletasks()                      # the largest tab decides the size; not beyond the screen
        w.geometry(f"{min(w.winfo_reqwidth(), w.winfo_screenwidth() - 40)}x{min(w.winfo_reqheight(), w.winfo_screenheight() - 90)}")
        w.deiconify()

    def regenerate(self):
        if not self.structure or not os.path.isfile(self.structure):
            messagebox.showinfo("soapff-GUI", "構造ファイルがわかりません。「新規計算」から選んでください")
            return
        if self.busy():
            return
        if self.lines and not messagebox.askyesno("soapff-GUI", "MD の条件から入力を作り直します。手で直した行は消えます。よろしいですか？"):
            return
        self.make_input(self.structure, self.out or os.path.splitext(self.structure)[0] + "_soapff")

    # ================================================================================================== run
    def fetch_dialog(self):
        """copy a folder (with everything in it, like cp -r) from a registered computer to this one."""
        if not self.s["computers"]:
            messagebox.showinfo("soapff-GUI", "計算設定の「実行先」で、先にコンピュータを登録してください")
            return
        w = tk.Toplevel(self.root)
        w.title("別のコンピュータからフォルダを受け取る")
        w.geometry("820x560")
        w.configure(bg=PANEL)
        f = ttk.Frame(w, padding=14)
        f.pack(fill="both", expand=True)
        f.columnconfigure(1, weight=1)
        names = [c["name"] for c in self.s["computers"]]
        v_c = tk.StringVar(value=self.s["where"] if self.s["where"] in names else names[0])
        v_path = tk.StringVar(value=(self.computer(v_c.get()) or {}).get("root", "soapff_runs"))
        v_dest = tk.StringVar(value=os.path.dirname(self.out) if self.out else self.s["last_dir"])
        ttk.Label(f, text="コンピュータ").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Combobox(f, textvariable=v_c, values=names, state="readonly", width=30).grid(row=0, column=1, sticky="w")
        ttk.Label(f, text="相手側のフォルダ").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(f, textvariable=v_path, style="Field.TEntry").grid(row=1, column=1, sticky="ew")
        lst = tk.Listbox(f, font=self.f_mono_s, height=14, activestyle="none", selectbackground=SEL_BG, selectforeground=TEXT)
        lst.grid(row=3, column=0, columnspan=3, sticky="nsew", pady=(6, 6))
        f.rowconfigure(3, weight=1)
        ttk.Label(f, text="受け取り先（このコンピュータ）").grid(row=4, column=0, sticky="w", pady=4)
        ttk.Entry(f, textvariable=v_dest, style="Field.TEntry").grid(row=4, column=1, sticky="ew")
        ttk.Button(f, text="選ぶ…", style="Ghost.TButton",
                   command=lambda: v_dest.set(filedialog.askdirectory(parent=w) or v_dest.get())).grid(row=4, column=2, padx=6)
        res = tk.StringVar(value="「一覧」で相手側のフォルダの中身を表示します。フォルダはダブルクリックで中に入ります。")
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
            res.set("一覧を取得しています…")
            path = v_path.get().strip() or "."

            def go():
                try:
                    items = r.list_dir(path)
                    msg = f"{path}：フォルダ {sum(d for _, d in items)}、ファイル {sum(not d for _, d in items)}"
                except Exception as e:
                    items, msg = [], f"一覧を取得できませんでした：{e}"

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
                res.set("受け取り先を選んでください")
                return
            target = os.path.join(dest, os.path.basename(src.rstrip("/")))
            if os.path.exists(target) and not messagebox.askyesno("soapff-GUI", f"{target} はもうあります。上書きしますか？",
                                                                  parent=w):
                return
            res.set(f"{src} を受け取っています（scp -r）…")

            def go():
                try:
                    got = r.fetch(src, dest)
                    n = sum(len(fs) for _, _, fs in os.walk(got)) if os.path.isdir(got) else 1
                    msg = f"受け取りました：{got}（{n} ファイル）"
                    if os.path.isfile(os.path.join(got, "in.lammps")):
                        msg += "。このフォルダは「ファイル → 既存の出力フォルダを開く」で開けます"
                except Exception as e:
                    msg = f"受け取れませんでした：{e}"
                self.q.put(("call", lambda: res.set(msg)))
            threading.Thread(target=go, daemon=True).start()

        def close():
            if st["r"]:
                st["r"].close()
            w.destroy()
        bar = ttk.Frame(f)
        bar.grid(row=2, column=1, sticky="w")
        ttk.Button(bar, text="一覧", style="Ghost.TButton", command=listing).pack(side="left")
        ttk.Button(bar, text="上のフォルダへ", style="Ghost.TButton", command=up).pack(side="left", padx=6)
        lst.bind("<Double-1>", enter)
        b2 = ttk.Frame(f)
        b2.grid(row=6, column=0, columnspan=3, sticky="ew")
        ttk.Button(b2, text="受け取る（選んだフォルダを丸ごと。選ばなければ上のフォルダ）", style="Accent.TButton",
                   command=fetch).pack(side="right")
        ttk.Button(b2, text="閉じる", style="Ghost.TButton", command=close).pack(side="left")
        w.protocol("WM_DELETE_WINDOW", close)

    def computer(self, name=None):
        name = name or self.s.get("where")
        return next((c for c in self.s["computers"] if c["name"] == name), None)

    def where_label(self):
        c = self.computer() if self.s.get("where") != "local" else None
        if c:
            return f"{c['name']}（{c['user'] + '@' if c['user'] else ''}{c['host']}）"
        return "このコンピュータ"

    def choose_where(self):
        v = self.v_where.get()
        self.s["where"] = "local" if v == "このコンピュータ" else v
        save_settings(self.s)
        self.update_status()

    def edit_computer(self, parent, c, taken=()):
        """a dialog to register / edit a computer; returns the new dict or None."""
        w = tk.Toplevel(parent)
        w.title("コンピュータの登録")
        w.configure(bg=PANEL)
        w.transient(parent)
        f = ttk.Frame(w, padding=16)
        f.pack(fill="both", expand=True)
        v = {}
        rows = (("name", "名前", "このアプリの中での呼び名（例：研究室のサーバー）"), ("host", "ホスト", "例：192.168.0.10 やホスト名"),
                ("user", "ユーザー", ""), ("port", "ポート", "SSH のポート（既定 22）"),
                ("root", "相手側のフォルダ", "計算を置くフォルダ（ホームからの相対パス可）"),
                ("lmp", "相手側の lmp", "名前かフルパス"), ("mpirun", "相手側の mpirun", "名前かフルパス"),
                ("setup", "最初に実行するコマンド", "例：module load lammps、source ~/.bashrc"), ("nproc", "並列数", "相手の物理コア数以下"))
        for i, (k, lab, note) in enumerate(rows):
            ttk.Label(f, text=lab).grid(row=i, column=0, sticky="w", pady=4, padx=(0, 10))
            v[k] = tk.StringVar(value=str(c.get(k, "")))
            ttk.Entry(f, textvariable=v[k], width=36 if k not in ("port", "nproc") else 8, style="Field.TEntry").grid(
                row=i, column=1, sticky="w", pady=4)
            ttk.Label(f, text=note, foreground=MUTED).grid(row=i, column=2, sticky="w", padx=8)
        ttk.Label(f, text="パスワードは登録しません。接続するたびに聞きます。", foreground=MUTED).grid(
            row=len(rows), column=0, columnspan=3, sticky="w", pady=(10, 0))
        out = {}

        def ok():
            d = {k: v[k].get().strip() for k in v}
            try:
                d["port"] = int(d["port"] or 22)
                d["nproc"] = int(d["nproc"] or 1)
            except ValueError:
                messagebox.showerror("soapff-GUI", "ポートと並列数は整数にしてください", parent=w)
                return
            if not d["name"] or not d["host"]:
                messagebox.showerror("soapff-GUI", "名前とホストを入れてください", parent=w)
                return
            if d["name"] in taken or d["name"] == "このコンピュータ":
                messagebox.showerror("soapff-GUI", f"「{d['name']}」はもう使われています", parent=w)
                return
            out.update(dict(COMPUTER, **d))
            w.destroy()
        bar = ttk.Frame(f)
        bar.grid(row=len(rows) + 1, column=0, columnspan=3, sticky="e", pady=(14, 0))
        ttk.Button(bar, text="取消", style="Ghost.TButton", command=w.destroy).pack(side="left", padx=6)
        ttk.Button(bar, text="登録", style="Accent.TButton", command=ok).pack(side="left")
        w.grab_set()
        w.wait_window()
        return out or None

    def update_status(self):
        r = self.report
        st = f"構造: {os.path.basename(self.structure) if self.structure else '—'}"
        if r:
            st += f"   |   {r['atoms']:,} atoms"
        self.l_status.config(text=st)
        self.l_status_r.config(text=f"実行先: {self.where_label()}   ·   LAMMPS · SOAP")
        names = ["このコンピュータ"] + [c["name"] for c in self.s["computers"]]
        self.cb_where.configure(values=names)
        self.v_where.set(self.s["where"] if self.s.get("where") in names else "このコンピュータ")

    def busy(self):
        return bool(self.worker and self.worker.is_alive())

    def start_or_stop(self):
        if self.busy():
            self.stop()
        else:
            self.start()

    def ask_password(self, c, parent=None):
        """the SSH password of a registered computer, asked every time (never stored).  None: cancelled."""
        return simpledialog.askstring("SSH", f"{c['name']}（{c['user'] + '@' if c['user'] else ''}{c['host']}:{c['port']}）の"
                                             "パスワード", show="*", parent=parent or self.root)

    def make_remote(self, c, password=None):
        return remote.Remote(c["host"], c.get("user", ""), c.get("port", 22), "", c.get("root", "soapff_runs"),
                             c.get("lmp", "lmp"), c.get("mpirun", "mpirun"), c.get("setup", ""), password=password or None)

    def start(self):
        if self.busy():
            return
        if not self.docs["in.lammps"] or not self.out:
            messagebox.showinfo("soapff-GUI", "先に「新規計算」で構造ファイルを選んで、入力を作ってください")
            return
        if not L.run_steps(self.docs["in.lammps"]) and not messagebox.askyesno(
                "soapff-GUI", "in.lammps に run がありません（エネルギー最小化だけです）。このまま回しますか？\n"
                              "MD を回すには、＋ fix（nvt など）と ＋ その他 の run を追加してください。"):
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
                messagebox.showerror("soapff-GUI", "LAMMPS (lmp) が見つかりません。計算設定で場所を指定してください")
                return
            mpirun = runner.which(s["mpirun"]) if int(s["nproc"]) > 1 else None
            if int(s["nproc"]) > 1 and not mpirun:
                messagebox.showerror("soapff-GUI", "mpirun が見つかりません。並列数を 1 にするか、計算設定で場所を指定してください")
                return
        self.stopping = False
        self.log_clear()
        self.set_state("開始", ACCENT, 0)
        self.set_running(True)
        self.remember("計算中")
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
                put(("log", f"WARNING {nproc} 並列のはずが「1 by 1 by 1」です：1 並列の LAMMPS が {nproc} 個別々に動いています。"
                            "lmp と同じ MPI の mpirun を使ってください（計算設定の実行先の「相手側の mpirun」をフルパスで、"
                            "lmp と同じフォルダのものに）\n", "warn"))
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
            put(("done", f"エラー：{e}", 1))
        finally:
            self.job = None

    def work_remote(self, s, comp, password, lines):
        put = self.q.put
        name = os.path.basename(os.path.normpath(self.out))
        r = None
        try:
            r = self.make_remote(comp, password)
            put(("phase", f"{comp['name']} に接続中"))
            info = r.check()
            put(("log", f"接続：{info.get('host')}（{info.get('os')}、CPU {info.get('cpus')}、物理コア {info.get('cores') or '?'}）\n", None))
            put(("log", f"lmp：{info.get('lmp')}（{info.get('lmp_mpi') or 'MPI 不明'}）\nmpirun：{info.get('mpirun_used')}\n", None))
            if info.get("mpi_warning") and int(comp.get("nproc", 1)) > 1:
                put(("log", "WARNING " + info["mpi_warning"] + "\n", "warn"))
            rdir = r.upload(self.out, name)
            put(("log", f"入力を送りました：{r.target}:{rdir}\n", None))
            nproc = int(comp.get("nproc", 1))
            cores = int(info.get("cores") or 0) if str(info.get("cores") or "").isdigit() else 0
            cpus = int(info.get("cpus") or 0) if str(info.get("cpus") or "").isdigit() else 0
            if cpus and nproc > cpus:
                put(("log", f"WARNING 並列数 {nproc} が相手の CPU（スレッド）数 {cpus} より多いので、{cpus} 並列にします\n", "warn"))
                nproc = cpus
            elif cores and nproc > cores:
                put(("log", f"並列数 {nproc} は物理コア数 {cores} より多いので、ハードウェアスレッドも使います"
                            "（Open MPI では --use-hwthread-cpus。LAMMPS はコア数を超えるとあまり速くなりません）\n", None))
            comp = dict(comp, nproc=nproc)
            pid = r.start(name, nproc)
            self.rjob = (r, name, pid)
            put(("log", f"LAMMPS を起動しました（{comp['name']}、pid {pid}、{comp.get('nproc')} 並列）。"
                        "接続が切れても計算は続きます\n", None))
            put(("phase", "MD（リモート）"))
            self._follow(r.follow(name, pid), L.StepTracker(lines), put, int(comp.get("nproc", 1)))
            rc = r.exit_code(name)
            put(("phase", "結果を受け取り中"))
            got = r.download_all(name, self.out)
            put(("log", f"相手側のフォルダを丸ごと受け取りました（{len(got)} ファイル）：{', '.join(got[:12])}"
                        f"{' …' if len(got) > 12 else ''}\n", None))
            self._finish(rc if rc is not None else (0 if self.stopping else 1), s, put)
        except Exception as e:
            put(("log", f"ERROR {type(e).__name__}: {e}\n", "err"))
            put(("done", f"エラー：{e}", 1))
        finally:
            self.rjob = None
            if r:
                r.close()                 # the password is asked again next time

    def _finish(self, rc, s, put):
        if self.stopping:
            put(("done", "停止しました", rc if rc else 1))
        elif rc == 0:
            if s.get("auto_analyze") and runner.analyzable(self.out):
                self._analyze(put, s)
            put(("done", "完了", 0))
        else:
            put(("done", f"LAMMPS が異常終了しました（終了コード {rc}）。ログと log.lammps を確認してください", rc or 1))

    def stop(self):
        if not self.busy():
            return
        self.stopping = True
        self.l_phase.config(text="停止中…")
        if self.job:
            self.job.stop()
        if self.rjob:
            r, name, pid = self.rjob
            threading.Thread(target=lambda: r.stop(pid), daemon=True).start()

    # ================================================================================================== analysis
    def analysis_form(self, f, v, s, auto=False):
        """the analysis choices and their parameters (settings tab and the 解析 dialog); values in v."""
        f.columnconfigure(1, weight=1)
        row = 0
        if auto:
            v["auto_analyze"] = tk.BooleanVar(value=bool(s.get("auto_analyze", True)))
            ttk.Checkbutton(f, text="計算が終わったら解析する", variable=v["auto_analyze"]).grid(row=row, column=0, columnspan=4,
                                                                                               sticky="w", pady=(0, 6))
            row += 1
        ttk.Label(f, text="解析", style="Section.TLabel").grid(row=row, column=0, columnspan=4, sticky="w")
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
                ttk.Button(rf, text="選ぶ…", style="Ghost.TButton",
                           command=lambda: v[key].set((filedialog.askdirectory() if browse == "dir" else
                                                       filedialog.askopenfilename()) or v[key].get())).pack(side="left", padx=6)
            if note:
                ttk.Label(rf, text=note, foreground=MUTED).pack(side="left", padx=6)
            row += 1
        ent("skip_ps", "平衡化として除く時間 (ps)", 8, note="空欄で最初の 20 %")
        ttk.Label(f, text="XRD（debyer）", style="Section.TLabel").grid(row=row, column=0, columnspan=4, sticky="w", pady=(10, 0))
        row += 1
        ent("xrd_wavelength", "波長 λ (Å)", 8, note="Cu Kα 1.5406、Mo Kα 0.7107")
        ttk.Label(f, text="2θ の範囲 (°)").grid(row=row, column=0, sticky="w", pady=3)
        rg = ttk.Frame(f)
        rg.grid(row=row, column=1, columnspan=3, sticky="w")
        for i, (k, lab) in enumerate((("xrd_from", "から"), ("xrd_to", "まで"), ("xrd_step", "刻み"))):
            v[k] = tk.StringVar(value=str(s.get(k, "")))
            ttk.Entry(rg, textvariable=v[k], width=6, style="Field.TEntry").pack(side="left")
            ttk.Label(rg, text=lab, foreground=MUTED).pack(side="left", padx=(3, 10))
        row += 1
        ent("xrd_frames", "平均するフレーム数", 6, note="最後のフレームから")
        ttk.Label(f, text="方式").grid(row=row, column=0, sticky="nw", pady=3)
        v["xrd_mode"] = tk.StringVar(value=s.get("xrd_mode", "cluster"))
        mf = ttk.Frame(f)
        mf.grid(row=row, column=1, columnspan=3, sticky="w")
        for k, lab in XRD_MODES.items():
            ttk.Radiobutton(mf, text=lab, value=k, variable=v["xrd_mode"]).pack(anchor="w")
        row += 1
        ent("xrd_size", "結晶子の大きさ (Å)", 8, note="各辺がこれ以上になるようセルを並べる（時間は原子数の 2 乗）")
        ttk.Label(f, text="プログラム", style="Section.TLabel").grid(row=row, column=0, columnspan=4, sticky="w", pady=(10, 0))
        row += 1
        ent("sova_python", "SOVA の Python", 40, browse="file")
        ent("sova_home", "SOVA のフォルダ", 40, browse="dir")
        ent("debyer", "debyer", 40, browse="file", note="" if runner.find_debyer(s.get("debyer")) else "見つかりません")

    @staticmethod
    def collect_analysis(v, s2):
        """the values of analysis_form into the settings s2 (ValueError when a number is wrong)."""
        if "auto_analyze" in v:
            s2["auto_analyze"] = bool(v["auto_analyze"].get())
        s2["analyses"] = [k for k, _ in ANALYSES if v["an_" + k].get()]
        for k in ("skip_ps", "sova_python", "sova_home", "debyer", "xrd_mode"):
            s2[k] = v[k].get().strip()
        try:
            if s2["skip_ps"]:
                float(s2["skip_ps"])
            for k in ("xrd_wavelength", "xrd_from", "xrd_to", "xrd_step", "xrd_size"):
                s2[k] = float(v[k].get())
            s2["xrd_frames"] = max(1, int(float(v["xrd_frames"].get())))
        except ValueError as e:
            raise ValueError(f"数値ではありません：{e}")
        if not (0 < s2["xrd_from"] < s2["xrd_to"] < 180 and s2["xrd_step"] > 0 and s2["xrd_wavelength"] > 0):
            raise ValueError("XRD の 2θ の範囲か波長が正しくありません")

    def _analyze(self, put, s):
        what = [k for k, _ in ANALYSES if k in s.get("analyses", [])]
        put(("phase", "解析中"))
        put(("log", "\n解析：" + "、".join(lab.split("（")[0] for k, lab in ANALYSES if k in what) + "\n", None))
        if not what:
            put(("log", "解析が 1 つも選ばれていません（計算設定 › 解析）\n", "warn"))
            return
        try:
            skip = s.get("skip_ps")
            runner.analyze(self.out, log=lambda m: put(("log", m + "\n", None)), skip_ps=float(skip) if skip else None,
                           sova_python=s.get("sova_python"), sova_home=s.get("sova_home"), what=what,
                           debyer=s.get("debyer"), traj=s.get("_traj"),
                           xrd_opts=dict(wavelength=float(s.get("xrd_wavelength", 1.5406)),
                                         two_theta=(float(s.get("xrd_from", 3)), float(s.get("xrd_to", 80)),
                                                    float(s.get("xrd_step", 0.02))),
                                         frames=int(s.get("xrd_frames", 3)), mode=s.get("xrd_mode", "cluster"),
                                         size_A=float(s.get("xrd_size", 50))))
            put(("analysis", self.out))
        except Exception as e:
            put(("log", f"解析できませんでした：{e}\n", "warn"))

    def analyze_now(self):
        if not runner.analyzable(self.out):
            messagebox.showinfo("soapff-GUI", "解析できる出力（軌道 *.lammpstrj、MSD、物性の計算や fix ave/* のファイル）が"
                                              "まだありません。先に計算してください")
            return
        if self.busy():
            messagebox.showinfo("soapff-GUI", "計算中です。終わってから解析してください")
            return
        w = tk.Toplevel(self.root)
        w.title(f"解析 – {os.path.basename(self.out)}")
        w.configure(bg=PANEL)
        w.transient(self.root)
        from soapff_cui import analysis as AN
        trajs = AN.trajectories(self.out)
        v = {}
        if trajs:
            tf = ttk.Frame(w, padding=(16, 14, 16, 0))
            tf.pack(fill="x")
            ttk.Label(tf, text="SOVA と XRD に使う軌道").pack(side="left", padx=(0, 10))
            v["traj"] = tk.StringVar(value=trajs[0])
            ttk.Combobox(tf, textvariable=v["traj"], values=trajs, state="readonly", width=34).pack(side="left")
            ttk.Label(tf, text="（引張・せん断の軌道は変形の途中の構造）", foreground=MUTED).pack(side="left", padx=8)
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
        ttk.Button(bar, text="解析する", style="Accent.TButton", command=go).pack(side="right")
        if os.path.isfile(os.path.join(self.out, "analysis", "summary.txt")):
            ttk.Button(bar, text="前回の結果を見る", style="Ghost.TButton",
                       command=lambda: (w.destroy(), self.show_analysis(self.out))).pack(side="right", padx=8)
        ttk.Button(bar, text="閉じる", style="Ghost.TButton", command=w.destroy).pack(side="left")

    def show_results(self):
        if self.out and os.path.isfile(os.path.join(self.out, "analysis", "summary.txt")):
            self.show_analysis(self.out)
        elif runner.analyzable(self.out):
            self.analyze_now()
        else:
            messagebox.showinfo("soapff-GUI", "まだ結果がありません。「計算を開始」で計算してください")

    def show_analysis(self, out):
        """the summary and every figure of out/analysis in tabs."""
        dest = os.path.join(out, "analysis")
        try:
            figs = [tuple(x) for x in json.load(open(os.path.join(dest, "summary.json"))).get("figures", [])]
        except (OSError, ValueError):
            figs = []
        if not figs and os.path.isfile(os.path.join(dest, "analysis.png")):     # analysis of an earlier version
            figs = [("MSD・PDF・RDF", "analysis.png")]
        w = tk.Toplevel(self.root)
        w.title(f"解析結果 – {os.path.basename(out)}")
        w.configure(bg=PANEL)
        nb = ttk.Notebook(w)
        nb.pack(fill="both", expand=True, padx=8, pady=(8, 4))
        tab = ttk.Frame(nb, padding=4)
        nb.add(tab, text="まとめ")
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
                ttk.Label(tab, text=f"図を表示できません：{e}").pack()

        def current():
            return pngs[nb.index(nb.select())] or os.path.join(dest, "summary.txt")
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=8, pady=8)
        ttk.Button(bar, text="この図を開く", style="Ghost.TButton", command=lambda: open_folder(current())).pack(side="left")
        ttk.Button(bar, text="解析フォルダを開く（CSV）", style="Ghost.TButton", command=lambda: open_folder(dest)).pack(side="left", padx=6)
        ttk.Button(bar, text="もう一度解析…", style="Ghost.TButton", command=lambda: (w.destroy(), self.analyze_now())).pack(side="left")
        ttk.Button(bar, text="閉じる", style="Ghost.TButton", command=w.destroy).pack(side="right")

    def open_ovito(self):
        if not self.out or not os.path.isdir(self.out):
            messagebox.showinfo("soapff-GUI", "出力フォルダがまだありません")
            return

        def go():
            try:
                p = runner.open_in_ovito(self.out)
                self.q.put(("log", f"OVITO で開きました：{p}\n", None))
            except Exception as e:
                self.q.put(("error", f"OVITO で開けませんでした：{e}"))
        threading.Thread(target=go, daemon=True).start()

    def open_out(self):
        if self.out and os.path.isdir(self.out):
            open_folder(self.out)

    def open_readme(self):
        p = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "README.md")
        if os.path.isfile(p):
            open_folder(p)

    def about(self):
        messagebox.showinfo("soapff-GUI", f"soapff-GUI v{__version__}\n構造ファイルから ClayFF（SOAP で電荷と型）の LAMMPS 入力を作り、"
                                          "このコンピュータか SSH の先で計算して、拡散係数と PDF・RDF を解析します。")

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
                elif kind == "made":
                    self.set_state("入力を作成しました（最小化まで）", OK, 0)
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
                    self.log_add("入力は最小化までです。＋ fix、＋ compute、＋ その他（run など）から MD を追加してください。\n", None)
                    self.remember("入力作成")
                    self.warn_assignment(msg[3])
                elif kind == "make_failed":
                    self.set_state("力場を作れませんでした", ERR, 0)
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
                    self.set_state(msg[1] if len(msg[1]) < 40 else ("完了" if ok else "終了（エラー）"), OK if ok else ERR,
                                   1.0 if ok else None)
                    self.log_add(f"\n== {msg[1]}\n", None if ok else "err")
                    self.remember("完了" if ok else ("停止" if self.stopping else "エラー"))
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
        txt = "停止" if on else "計算を開始"
        self.b_start.config(text=txt)
        self.b_start_top.config(text=txt)
        if on:
            self.dot.itemconfigure(self.dot_id, fill=ACCENT)

    def log_add(self, text, tag):
        self.log.config(state="normal")
        if self.log.get("1.0", "1.end").startswith("実行ログはここに"):
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
            if not messagebox.askyesno("soapff-GUI", "計算中です。止めて終了しますか？\n（SSH の計算は相手側で続きます）"):
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
