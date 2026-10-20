"""soapff-GUI: choose a clay structure, press Run.  soapff writes the LAMMPS input in the output folder and LAMMPS runs
it; the user never edits an input file."""
import json
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from . import runner

SETTINGS = os.path.join(os.path.expanduser("~"), ".config", "soapff-gui", "settings.json")
STRUCTURE_TYPES = [("構造ファイル", "*.xyz *.extxyz *.cif *.vasp *.data *.lammps *.lmp POSCAR* CONTCAR*"), ("すべて", "*")]
BARO_HELP = {"z": "z：面に垂直な方向だけ（界面・スラブ・層間）", "aniso": "aniso：3 辺を別々に（バルク結晶）",
             "tri": "tri：箱の傾きも含めて全部（バルク結晶）", "iso": "iso：等方"}
# per axis: p, or the lo/hi faces from f / s / m (one letter = both faces), as in the LAMMPS boundary command
AXIS_STYLES = ["p", "f", "m", "s", "fm", "mf", "fs", "sf", "ms", "sm"]
FACE_JP = {"f": "固定", "s": "原子に合わせて縮む", "m": "原子に合わせて広がる（初めの大きさより小さくならない）"}


def boundary_help(b):
    """short Japanese description of a boundary such as ['p', 'p', 'f']."""
    if b == ["p", "p", "p"]:
        return "3 方向とも周期（層間・水で埋めた界面・バルク）。PPPM"
    open_ = [ax for ax, x in zip("xyz", b) if x != "p"]
    faces = sorted({c for x in b if x != "p" for c in x})
    txt = f"{'・'.join(open_)} は周期なし（{'、'.join(f'{c}：{FACE_JP[c]}' for c in faces)}）。NVT のみ"
    if b[:2] == ["p", "p"] and b[2] == "f":
        txt += "。直方体なら PPPM＋slab 補正、傾いた箱は MSM"
    else:
        txt += "。MSM（PPPM の 10 倍以上遅い）"
    return txt


PHASE_JP = {"start": "準備中", "minimize": "構造最小化", "relax": "緩和 NVT", "md": "MD", "analysis": "解析中"}
DEFAULTS = dict(temp=300.0, relax_ps=10.0, nvt_ps=50.0, npt_ps=0.0, press=1.0, baro="aniso", boundary="p p p", nproc=max(1, min(4, (os.cpu_count() or 2) // 2)),
                threshold=0.95, cutoff=10.0, library="", lmp="", mpirun="", last_dir=os.path.expanduser("~"),
                dump_ps=1.0, auto_analyze=True, skip_ps="", sova_python="", sova_home="")


def load_settings():
    s = dict(DEFAULTS)
    try:
        s.update(json.load(open(SETTINGS)))
    except (OSError, ValueError):
        pass
    if not os.path.isfile(s["library"]):
        s["library"] = runner.default_library()
    lmp, mpirun = runner.find_executables()
    if not s.get("sova_python"):
        s["sova_python"], s["sova_home"] = runner.find_sova()
    s["lmp"] = s["lmp"] or lmp or ""
    s["mpirun"] = s["mpirun"] or mpirun or ""
    return s


def save_settings(s):
    try:
        os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
        json.dump(s, open(SETTINGS, "w"), indent=1)
    except OSError:
        pass


def open_folder(path):
    if sys.platform == "darwin":
        subprocess.Popen(["open", path])
    elif sys.platform.startswith("win"):
        os.startfile(path)
    else:
        subprocess.Popen(["xdg-open", path])


class App:
    def __init__(self, root):
        self.root = root
        self.s = load_settings()
        self.q = queue.Queue()
        self.job = None
        self.worker = None
        self.stopping = False
        root.title("soapff-GUI")
        root.minsize(760, 640)
        self._build()
        root.protocol("WM_DELETE_WINDOW", self.on_close)
        root.after(100, self.poll)
        try:
            runner.import_soapff()
        except ImportError as e:
            self.l_struct.config(text=f"soapff-cui が見つかりません：{e}", foreground="#c00000")
            root.after(200, lambda: messagebox.showerror("soapff-GUI", f"soapff-cui が見つかりません。\n{e}"))

    # ------------------------------------------------------------------------------------------------------- layout
    def _build(self):
        f = ttk.Frame(self.root, padding=10)
        f.pack(fill="both", expand=True)
        f.columnconfigure(1, weight=1)

        self.v_struct = tk.StringVar()
        self.v_out = tk.StringVar()
        ttk.Label(f, text="初期構造").grid(row=0, column=0, sticky="w")
        ttk.Entry(f, textvariable=self.v_struct).grid(row=0, column=1, sticky="ew", padx=4)
        ttk.Button(f, text="選ぶ…", command=self.pick_structure).grid(row=0, column=2)
        self.l_struct = ttk.Label(f, text="extxyz（Lattice 付き）、cif、POSCAR、LAMMPS データ。元素・座標・箱だけで足ります",
                                  foreground="gray")
        self.l_struct.grid(row=1, column=1, columnspan=2, sticky="w", padx=4)
        ttk.Label(f, text="出力フォルダ").grid(row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(f, textvariable=self.v_out).grid(row=2, column=1, sticky="ew", padx=4, pady=(6, 0))
        ttk.Button(f, text="選ぶ…", command=self.pick_out).grid(row=2, column=2, pady=(6, 0))

        md = ttk.LabelFrame(f, text="MD の条件", padding=8)
        md.grid(row=3, column=0, columnspan=3, sticky="ew", pady=8)
        self.v = {}
        fields = [("temp", "温度 (K)"), ("relax_ps", "緩和 NVT (ps)"), ("nvt_ps", "NVT (ps)"), ("npt_ps", "NPT (ps, 0 でなし)"), ("press", "圧力 (atm)"),
                  ("dump_ps", "出力間隔 (ps)")]
        for i, (k, lab) in enumerate(fields):
            ttk.Label(md, text=lab).grid(row=0, column=2 * i, sticky="e", padx=(8 if i else 0, 2))
            self.v[k] = tk.StringVar(value=str(self.s[k]))
            ttk.Entry(md, textvariable=self.v[k], width=8).grid(row=0, column=2 * i + 1, sticky="w")
        ttk.Label(md, text="圧力の方向").grid(row=1, column=0, sticky="e", pady=(6, 0), padx=(0, 2))
        self.v["baro"] = tk.StringVar(value=self.s["baro"])
        cb = ttk.Combobox(md, textvariable=self.v["baro"], values=list(BARO_HELP), width=6, state="readonly")
        cb.grid(row=1, column=1, sticky="w", pady=(6, 0))
        self.l_baro = ttk.Label(md, text=BARO_HELP[self.s["baro"]], foreground="gray")
        self.l_baro.grid(row=1, column=2, columnspan=4, sticky="w", pady=(6, 0), padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self.l_baro.config(text=BARO_HELP[self.v["baro"].get()]))
        ttk.Label(md, text="並列数").grid(row=1, column=6, sticky="e", pady=(6, 0), padx=(8, 2))
        self.v["nproc"] = tk.StringVar(value=str(self.s["nproc"]))
        ttk.Spinbox(md, from_=1, to=os.cpu_count() or 64, textvariable=self.v["nproc"], width=5).grid(
            row=1, column=7, sticky="w", pady=(6, 0))
        ttk.Label(md, text="境界条件").grid(row=2, column=0, sticky="e", pady=(6, 0), padx=(0, 2))
        bf = ttk.Frame(md)
        bf.grid(row=2, column=1, columnspan=7, sticky="w", pady=(6, 0))
        b0 = self.s["boundary"].split()
        b0 = b0 if len(b0) == 3 and all(x in AXIS_STYLES for x in b0) else ["p", "p", "p"]
        self.v_bnd = []
        for k, ax in enumerate("xyz"):
            ttk.Label(bf, text=ax).pack(side="left", padx=(0 if k == 0 else 6, 2))
            v = tk.StringVar(value=b0[k])
            ttk.Combobox(bf, textvariable=v, values=AXIS_STYLES, width=3, state="readonly").pack(side="left")
            v.trace_add("write", lambda *a: self.l_bnd.config(text=boundary_help([v.get() for v in self.v_bnd])))
            self.v_bnd.append(v)
        self.l_bnd = ttk.Label(md, text=boundary_help(b0), foreground="gray", wraplength=640, justify="left")
        self.l_bnd.grid(row=3, column=0, columnspan=8, sticky="w", pady=(2, 0))

        self.adv_shown = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, text="詳細設定", variable=self.adv_shown, command=self.toggle_adv).grid(
            row=4, column=0, columnspan=3, sticky="w")
        self.adv = ttk.Frame(f, padding=(16, 2, 0, 4))
        self.adv.columnconfigure(1, weight=1)
        for i, (k, lab, pick) in enumerate([("library", "SOAP ライブラリ", "file"), ("lmp", "LAMMPS (lmp)", "file"),
                                            ("mpirun", "mpirun", "file"), ("threshold", "類似度のしきい値", None),
                                            ("cutoff", "カットオフ (Å)", None),
                                            ("skip_ps", "解析で除く平衡化 (ps、空欄で 20%)", None),
                                            ("sova_python", "SOVA の Python（空欄でこの Python）", "file"),
                                            ("sova_home", "SOVA のフォルダ（ソースのとき）", "dir")]):
            ttk.Label(self.adv, text=lab).grid(row=i, column=0, sticky="w")
            self.v[k] = tk.StringVar(value=str(self.s[k]))
            ttk.Entry(self.adv, textvariable=self.v[k], width=10 if pick is None else 50).grid(
                row=i, column=1, sticky="ew" if pick else "w", padx=4)
            if pick == "dir":
                ttk.Button(self.adv, text="選ぶ…", command=lambda k=k: self.pick_dir(k)).grid(row=i, column=2)
            elif pick:
                ttk.Button(self.adv, text="選ぶ…", command=lambda k=k: self.pick_file(k)).grid(row=i, column=2)

        bar = ttk.Frame(f)
        bar.grid(row=6, column=0, columnspan=3, sticky="ew", pady=6)
        self.b_run = ttk.Button(bar, text="▶ 実行", command=self.run)
        self.b_run.pack(side="left")
        self.b_stop = ttk.Button(bar, text="■ 停止", command=self.stop, state="disabled")
        self.b_stop.pack(side="left", padx=4)
        self.b_inp = ttk.Button(bar, text="入力だけ作る", command=lambda: self.run(prepare_only=True))
        self.b_inp.pack(side="left", padx=4)
        self.v_auto = tk.BooleanVar(value=bool(self.s.get("auto_analyze", True)))
        ttk.Checkbutton(bar, text="終わったら解析", variable=self.v_auto).pack(side="left", padx=(8, 0))
        self.b_open = ttk.Button(bar, text="出力フォルダを開く", command=self.open_out)
        self.b_open.pack(side="right")
        self.b_ovito = ttk.Button(bar, text="OVITO で開く", command=self.open_ovito)
        self.b_ovito.pack(side="right", padx=4)
        self.b_an = ttk.Button(bar, text="解析", command=self.analyze_now)
        self.b_an.pack(side="right", padx=4)

        st = ttk.LabelFrame(f, text="状態", padding=8)
        st.grid(row=7, column=0, columnspan=3, sticky="ew")
        st.columnconfigure(1, weight=1)
        self.l_phase = ttk.Label(st, text="待機中", width=14)
        self.l_phase.grid(row=0, column=0, sticky="w")
        self.pb = ttk.Progressbar(st, maximum=1.0)
        self.pb.grid(row=0, column=1, sticky="ew", padx=6)
        self.l_pct = ttk.Label(st, text="", width=7)
        self.l_pct.grid(row=0, column=2)
        self.l_thermo = ttk.Label(st, text="", font=("Menlo", 11))
        self.l_thermo.grid(row=1, column=0, columnspan=3, sticky="w", pady=(4, 0))
        self.l_report = ttk.Label(st, text="", justify="left")
        self.l_report.grid(row=2, column=0, columnspan=3, sticky="w", pady=(4, 0))

        self.log = scrolledtext.ScrolledText(f, height=14, font=("Menlo", 10), state="disabled")
        self.log.grid(row=8, column=0, columnspan=3, sticky="nsew", pady=(8, 0))
        self.log.tag_config("warn", foreground="#b35900")
        self.log.tag_config("err", foreground="#c00000")
        f.rowconfigure(8, weight=1)

    def toggle_adv(self):
        if self.adv_shown.get():
            self.adv.grid(row=5, column=0, columnspan=3, sticky="ew")
        else:
            self.adv.grid_forget()

    # ------------------------------------------------------------------------------------------------------ inputs
    def pick_structure(self):
        p = filedialog.askopenfilename(title="初期構造を選ぶ", initialdir=self.s["last_dir"], filetypes=STRUCTURE_TYPES)
        if not p:
            return
        self.s["last_dir"] = os.path.dirname(p)
        self.v_struct.set(p)
        stem = os.path.splitext(os.path.basename(p))[0]
        self.v_out.set(os.path.join(os.path.dirname(p), f"{stem}_soapff"))
        try:
            d = runner.structure_summary(p)
            el = "  ".join(f"{e} {n}" for e, n in d["elements"].items())
            self.l_struct.config(text=f"{d['atoms']} 原子（{el}）  箱 {d['cell'][0]:.2f} × {d['cell'][1]:.2f} × {d['cell'][2]:.2f} Å",
                                 foreground="")
        except Exception as e:
            self.l_struct.config(text=f"読めません：{e}", foreground="#c00000")

    def pick_out(self):
        p = filedialog.askdirectory(title="出力フォルダを選ぶ", initialdir=self.s["last_dir"])
        if p:
            self.v_out.set(p)

    def pick_file(self, key):
        p = filedialog.askopenfilename(title=key)
        if p:
            self.v[key].set(p)

    def pick_dir(self, key):
        p = filedialog.askdirectory(title=key)
        if p:
            self.v[key].set(p)

    def open_out(self):
        d = self.v_out.get()
        if d and os.path.isdir(d):
            open_folder(d)

    def open_ovito(self):
        """the result (traj.lammpstrj, else final.data / system.data) in OVITO; works during a run too (frames so far)."""
        d = self.v_out.get().strip()
        if not d or not os.path.isdir(d):
            messagebox.showerror("soapff-GUI", "出力フォルダがまだありません。先に実行してください")
            return
        def go():
            try:
                path = runner.open_in_ovito(d)
                self.q.put(("log", f"OVITO で開きました：{path}\n", None))
            except Exception as e:
                self.q.put(("error", f"OVITO で開けませんでした：{e}"))
        threading.Thread(target=go, daemon=True).start()   # the `ovito` command (Linux) blocks until OVITO closes

    def collect(self):
        """read and check the form; returns a settings dict or raises ValueError with a Japanese message."""
        s = dict(self.s)
        for k in ("temp", "relax_ps", "nvt_ps", "npt_ps", "press", "threshold", "cutoff", "dump_ps"):
            try:
                s[k] = float(self.v[k].get())
            except ValueError:
                raise ValueError(f"{k} が数値ではありません：{self.v[k].get()!r}")
        if s["dump_ps"] <= 0:
            raise ValueError("出力間隔は正の値にしてください（拡散係数には 0.1〜0.5 ps）")
        s["skip_ps"] = self.v["skip_ps"].get().strip()
        if s["skip_ps"]:
            try:
                float(s["skip_ps"])
            except ValueError:
                raise ValueError("解析で除く平衡化の時間は数値か空欄にしてください")
        s["auto_analyze"] = bool(self.v_auto.get())
        try:
            s["nproc"] = int(self.v["nproc"].get())
        except ValueError:
            raise ValueError("並列数は整数で指定してください")
        for k in ("baro", "library", "lmp", "mpirun", "sova_python", "sova_home"):
            s[k] = self.v[k].get().strip()
        b = [v.get() for v in self.v_bnd]
        s["boundary"] = " ".join(b)
        if s["npt_ps"] > 0 and b != ["p", "p", "p"]:
            raise ValueError("周期でない方向があるときは NPT を使えません（真空を含む箱の圧力は決まらないため）。NPT を 0 にしてください")
        if s["temp"] <= 0 or s["relax_ps"] < 0 or s["nvt_ps"] < 0 or s["npt_ps"] < 0 or s["nvt_ps"] + s["npt_ps"] <= 0:
            raise ValueError("温度は正、緩和・NVT・NPT の長さは 0 以上（NVT か NPT のどちらかは正）にしてください")
        if not os.path.isfile(self.v_struct.get()):
            raise ValueError("初期構造のファイルを選んでください")
        if not self.v_out.get().strip():
            raise ValueError("出力フォルダを指定してください")
        if not os.path.isfile(s["library"]):
            raise ValueError(f"SOAP ライブラリがありません：{s['library']}")
        return s

    # --------------------------------------------------------------------------------------------------------- run
    def run(self, prepare_only=False):
        try:
            s = self.collect()
        except ValueError as e:
            messagebox.showerror("soapff-GUI", str(e))
            return
        lmp = runner.which(s["lmp"])
        if not prepare_only and not lmp:
            messagebox.showerror("soapff-GUI", "LAMMPS (lmp) が見つかりません。詳細設定で場所を指定してください")
            return
        mpirun = runner.which(s["mpirun"]) if s["nproc"] > 1 else None
        if not prepare_only and s["nproc"] > 1 and not mpirun:
            messagebox.showerror("soapff-GUI", "mpirun が見つかりません。並列数を 1 にするか、詳細設定で場所を指定してください")
            return
        out = os.path.abspath(self.v_out.get().strip())
        if os.path.isfile(os.path.join(out, "in.lammps")) and not messagebox.askyesno(
                "soapff-GUI", f"{out} には前の結果があります。上書きしますか？"):
            return
        self.s = s
        save_settings(s)
        self.log_clear()
        self.set_running(True)
        self.l_report.config(text="")
        self.l_thermo.config(text="")
        self.pb["value"] = 0
        self.l_pct.config(text="")
        self.stopping = False
        md = {k: s[k] for k in ("temp", "relax_ps", "nvt_ps", "npt_ps", "press", "baro", "boundary", "threshold", "cutoff", "dump_ps")}
        self.worker = threading.Thread(target=self.work, args=(self.v_struct.get(), out, s, md, lmp, mpirun, prepare_only),
                                       daemon=True)
        self.worker.start()

    def work(self, structure, out, s, md, lmp, mpirun, prepare_only):
        """worker thread: soapff make, then LAMMPS.  Talks to the GUI only through the queue."""
        put = self.q.put
        try:
            put(("phase", "力場を作成中"))
            put(("log", f"soapff make {structure} -> {out}\n", None))
            rep = runner.prepare(structure, out, s["library"],
                                 log=lambda m: put(("log", m + "\n", "warn" if m.startswith("WARNING") else None)), **md)
            put(("report", rep))
            if prepare_only:
                put(("log", f"\n入力を書き出しました：{out}\n実行するには  cd {out} && lmp -in in.lammps\n", None))
                put(("done", "入力を作成しました", 0))
                return
            if self.stopping:
                put(("done", "停止しました", 1))
                return
            self.job = runner.LammpsJob(out, lmp, mpirun, s["nproc"])
            put(("log", "\n$ " + " ".join(self.job.cmd) + f"   (in {out})\n", None))
            tp = runner.ThermoParser(runner.total_md_steps(md["nvt_ps"], md["npt_ps"]), runner.md_steps(md["relax_ps"]))
            for line in self.job.lines():
                put(("log", line + "\n", "err" if line.startswith("ERROR") else ("warn" if line.startswith("WARNING") else None)))
                if tp.feed(line):
                    put(("thermo", tp.phase, dict(tp.values), tp.progress))
            rc = self.job.returncode
            if self.stopping:
                put(("done", "停止しました", rc))
            elif rc == 0:
                if s.get("auto_analyze"):
                    self._analyze(out, s.get("skip_ps"), put, s.get("sova_python"), s.get("sova_home"))
                put(("done", "完了", 0))
            else:
                put(("done", f"LAMMPS が異常終了しました（終了コード {rc}）。ログと {out}/log.lammps を確認してください", rc))
        except Exception as e:
            put(("log", f"\nERROR {type(e).__name__}: {e}\n", "err"))
            put(("done", f"エラー：{e}", 1))
        finally:
            self.job = None

    def _analyze(self, out, skip, put, sova_python=None, sova_home=None):
        """run the analyses (worker thread) and ask the GUI to show them."""
        put(("phase", "解析中"))
        put(("log", "\n解析（LAMMPS compute msd による拡散係数、SOVA による PDF と RDF）\n", None))
        try:
            runner.analyze(out, log=lambda m: put(("log", m + "\n", None)), skip_ps=float(skip) if skip else None,
                           sova_python=sova_python, sova_home=sova_home)
            put(("analysis", out))
        except Exception as e:
            put(("log", f"解析できませんでした：{e}\n", "warn"))

    def analyze_now(self):
        """analyse the output folder now (a finished run, or the frames written so far)."""
        d = self.v_out.get().strip()
        if not d or not os.path.isfile(os.path.join(d, "traj.lammpstrj")):
            messagebox.showerror("soapff-GUI", "出力フォルダにトラジェクトリ（traj.lammpstrj）がありません。先に MD を実行してください")
            return
        skip = self.v["skip_ps"].get().strip()
        sp, sh = self.v["sova_python"].get().strip(), self.v["sova_home"].get().strip()
        self.b_an.config(state="disabled")
        def go():
            self._analyze(os.path.abspath(d), skip, self.q.put, sp, sh)
            self.q.put(("analysis_end",))
        threading.Thread(target=go, daemon=True).start()

    def show_analysis(self, out):
        """a window with the summary and the figure of out/analysis."""
        dest = os.path.join(out, "analysis")
        w = tk.Toplevel(self.root)
        w.title(f"解析結果 – {os.path.basename(out)}")
        txt = scrolledtext.ScrolledText(w, height=12, width=110, font=("Menlo", 11))
        txt.insert("end", open(os.path.join(dest, "summary.txt")).read())
        txt.config(state="disabled")
        txt.pack(fill="x", padx=8, pady=(8, 4))
        png = os.path.join(dest, "analysis.png")
        try:
            img = tk.PhotoImage(file=png)
            k = max(1, -(-img.width() // 900))                  # shrink to <= 900 px wide
            img = img.subsample(k, k) if k > 1 else img
            lab = ttk.Label(w, image=img)
            lab.image = img
            lab.pack(padx=8)
        except tk.TclError as e:
            ttk.Label(w, text=f"図を表示できません：{e}").pack()
        bar = ttk.Frame(w)
        bar.pack(fill="x", padx=8, pady=8)
        ttk.Button(bar, text="図を開く", command=lambda: open_folder(png)).pack(side="left")
        ttk.Button(bar, text="解析フォルダを開く（CSV）", command=lambda: open_folder(dest)).pack(side="left", padx=4)
        ttk.Button(bar, text="閉じる", command=w.destroy).pack(side="right")

    def stop(self):
        self.stopping = True
        if self.job:
            self.job.stop()
        self.l_phase.config(text="停止中…")

    # ----------------------------------------------------------------------------------------------- GUI updates
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
                elif kind == "report":
                    self.show_report(msg[1])
                elif kind == "thermo":
                    self.show_thermo(*msg[1:])
                elif kind == "analysis":
                    self.show_analysis(msg[1])
                elif kind == "analysis_end":
                    self.b_an.config(state="normal")
                elif kind == "done":
                    self.set_running(False)
                    self.l_phase.config(text="完了" if msg[2] == 0 else "終了")
                    if msg[2] == 0 and msg[1] == "完了":
                        self.pb["value"] = 1.0
                        self.l_pct.config(text="100%")
                    self.log_add(f"\n== {msg[1]}\n", None if msg[2] == 0 else "err")
                    if msg[2] != 0 and not self.stopping:
                        messagebox.showerror("soapff-GUI", msg[1])
        except queue.Empty:
            pass
        self.root.after(100, self.poll)

    def show_thermo(self, phase, v, progress):
        self.l_phase.config(text=PHASE_JP.get(phase, phase))
        if phase in ("relax", "md"):
            self.pb["value"] = progress
            self.l_pct.config(text=f"{100 * progress:.0f}%")
        keys = [("Step", "step {:.0f}"), ("Temp", "T {:.1f} K"), ("Press", "P {:.0f} atm"), ("PotEng", "PE {:.1f}"),
                ("Density", "ρ {:.3f} g/cm³")]
        self.l_thermo.config(text="   ".join(fmt.format(v[k]) for k, fmt in keys if k in v))

    def show_report(self, r):
        roles = r["roles"]
        minsim = min((p["min_similarity"] for p in r["per_type"].values()), default=float("nan"))
        txt = (f"{r['atoms']} 原子（粘土 {roles['mineral']}、水 {roles['water']}、イオン {roles['ion']}）   "
               f"SOAP で電荷を決めた原子 {r['from_soap']}、fallback {r['fallback']}   最小類似度 {minsim:.3f}")
        notes = []
        if r["fallback"]:
            notes.append(f"⚠ {r['fallback']} 原子はライブラリに近い環境がなく、ClayFF の型の電荷を使いました（charges.csv）")
        if r["warnings"]:
            notes.append(f"⚠ 警告 {len(r['warnings'])} 件（ログを参照）")
        if r.get("boundary", "p p p") != "p p p":
            txt += f"\n境界条件 {r['boundary']}：" + "、".join(r["kspace"])
            if any("msm" in k for k in r["kspace"]):
                notes.append("⚠ 傾いた箱か周期でない方向があるので、長距離クーロン力は MSM で計算します（PPPM の 10 倍以上遅い）")
        self.l_report.config(text="\n".join([txt] + notes), foreground="#b35900" if notes else "")

    def log_add(self, text, tag):
        self.log.config(state="normal")
        self.log.insert("end", text, tag or ())
        if int(self.log.index("end-1c").split(".")[0]) > 20000:
            self.log.delete("1.0", "5000.0")
        self.log.see("end")
        self.log.config(state="disabled")

    def log_clear(self):
        self.log.config(state="normal")
        self.log.delete("1.0", "end")
        self.log.config(state="disabled")

    def set_running(self, on):
        self.b_run.config(state="disabled" if on else "normal")
        self.b_inp.config(state="disabled" if on else "normal")
        self.b_stop.config(state="normal" if on else "disabled")

    def on_close(self):
        if self.worker and self.worker.is_alive():
            if not messagebox.askyesno("soapff-GUI", "計算中です。LAMMPS を止めて終了しますか？"):
                return
            self.stop()
        self.root.destroy()


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
