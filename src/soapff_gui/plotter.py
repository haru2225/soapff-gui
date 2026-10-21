"""グラフ: choose columns of the output files of a run (the thermo output of log.lammps, the files written by fix ave/time
or fix print, csv files) and plot them, e.g. stress against strain (with the unit changed by a factor), temperature or
pressure against time, one curve from each of several runs (temperatures) in one figure.  Saves a PNG."""
import glob
import os
import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np

NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


def _is_num(w):
    return bool(NUM.match(w)) or w.lower() in ("nan", "inf", "-inf")


def read_blocks(path):
    """the tables of a text file -> [(name, columns, data (rows x cols))].
    log.lammps: every thermo block ('Step Temp ...' up to 'Loop time'); fix ave/time and fix print files: the columns named in
    the last '#' line; csv: the first line."""
    L = open(path, errors="replace").read().splitlines()
    out = []
    if any(l.startswith("Loop time") for l in L) or any(l.split()[:2] == ["Step", "Temp"] for l in L):
        cols, rows = None, []
        for l in L:
            w = l.split()
            if w[:1] == ["Step"] and not all(_is_num(x) for x in w):
                cols, rows = w, []
            elif cols is not None and l.startswith("Loop time"):
                if rows:
                    out.append((f"run {len(out) + 1}（step {int(rows[0][0])}–{int(rows[-1][0])}）", cols, np.array(rows)))
                cols = None
            elif cols is not None and w and len(w) == len(cols) and all(_is_num(x) for x in w):
                rows.append([float(x) for x in w])
        if out:
            return out
    if path.endswith(".csv"):
        head = L[0].split(",") if L else []
        data = []
        for l in L[1:]:
            w = l.split(",")
            if len(w) == len(head) and all(_is_num(x.strip()) for x in w):
                data.append([float(x) for x in w])
        if data:
            return [("csv", [h.strip() for h in head], np.array(data))]
    names, rows = None, []
    for l in L:
        if l.startswith("#"):
            w = l.lstrip("#").split()
            if len(w) >= 2 and not all(_is_num(x) for x in w):
                names = w
            continue
        w = l.split()
        if w and all(_is_num(x) for x in w):
            rows.append([float(x) for x in w])
    if rows:
        n = len(rows[0])
        rows = [r for r in rows if len(r) == n]
        data = np.array(rows)
        if names is None or len(names) != n:
            names = [f"列 {k + 1}" for k in range(n)]
        out.append(("data", names, data))
    return out


def candidate_files(folder):
    """the files of a run folder that hold tables."""
    fs = []
    for pat in ("log.lammps", "*.dat", "analysis/*.csv", "*.csv"):
        fs += sorted(glob.glob(os.path.join(folder, pat)))
    return fs


class PlotDialog:
    def __init__(self, app, folder):
        self.app, self.folder = app, folder
        self.cache = {}
        self.series = []                      # dicts: label, path, block, x, y, xs, ys
        w = self.w = tk.Toplevel(app.root)
        w.title("グラフ（出力ファイルの列を選んで描く）")
        sw, sh = app.root.winfo_screenwidth(), app.root.winfo_screenheight()
        w.geometry(f"{min(1250, sw - 80)}x{min(820, sh - 100)}")
        main = ttk.Frame(w, padding=8)
        main.pack(fill="both", expand=True)
        left = ttk.Frame(main)
        left.pack(side="left", fill="y")
        ttk.Label(left, text="1. ファイル", font=app.f_bold).pack(anchor="w")
        self.lb = tk.Listbox(left, height=9, width=38, exportselection=False, font=app.f_mono_s)
        self.lb.pack(anchor="w", pady=(2, 4))
        self.files = []
        self.refill()
        self.lb.bind("<<ListboxSelect>>", lambda e: self.on_file())
        ttk.Button(left, text="ファイルを足す…（別の計算の出力も）", style="Ghost.TButton", command=self.add_file).pack(anchor="w")
        ttk.Label(left, text="2. 表（log.lammps は run ごと）").pack(anchor="w", pady=(10, 0))
        self.v_block = tk.StringVar()
        self.cb_block = ttk.Combobox(left, textvariable=self.v_block, state="readonly", width=34)
        self.cb_block.pack(anchor="w")
        self.cb_block.bind("<<ComboboxSelected>>", lambda e: self.on_block())
        ttk.Label(left, text="3. 横軸 x と縦軸 y の列").pack(anchor="w", pady=(10, 0))
        self.v_x, self.v_y = tk.StringVar(), tk.StringVar()
        fx = ttk.Frame(left); fx.pack(anchor="w")
        ttk.Label(fx, text="x").grid(row=0, column=0)
        self.cb_x = ttk.Combobox(fx, textvariable=self.v_x, state="readonly", width=22); self.cb_x.grid(row=0, column=1, padx=4, pady=2)
        ttk.Label(fx, text="×").grid(row=0, column=2)
        self.v_xs = tk.StringVar(value="1")
        ttk.Entry(fx, textvariable=self.v_xs, width=9).grid(row=0, column=3, padx=4)
        ttk.Label(fx, text="y").grid(row=1, column=0)
        self.cb_y = ttk.Combobox(fx, textvariable=self.v_y, state="readonly", width=22); self.cb_y.grid(row=1, column=1, padx=4, pady=2)
        ttk.Label(fx, text="×").grid(row=1, column=2)
        self.v_ys = tk.StringVar(value="1")
        ttk.Entry(fx, textvariable=self.v_ys, width=9).grid(row=1, column=3, padx=4)
        ttk.Label(left, text="× は掛ける数（符号も変えられます）。圧力 atm → GPa は 1.01325e-4、応力 = −圧力なら −1.01325e-4",
                  foreground="#77848d", wraplength=330, justify="left").pack(anchor="w", pady=(2, 2))
        fp = ttk.Frame(left); fp.pack(anchor="w", pady=2)
        for txt, f in (("y: atm→GPa", "1.01325e-4"), ("y: −atm→GPa（応力）", "-1.01325e-4"), ("x: step→ps（0.5 fs）", "0.0005")):
            ttk.Button(fp, text=txt, style="Tiny.TButton", command=lambda t=txt, f=f: (self.v_xs if t.startswith("x") else self.v_ys).set(f)).pack(
                side="left", padx=(0, 4))
        self.v_label = tk.StringVar()
        ttk.Label(left, text="4. 凡例の名前").pack(anchor="w", pady=(8, 0))
        ttk.Entry(left, textvariable=self.v_label, width=36).pack(anchor="w")
        ttk.Button(left, text="＋ この曲線を図に足す", style="Accent.TButton", command=self.add_series).pack(anchor="w", pady=(8, 4))
        ttk.Label(left, text="図に入れた曲線").pack(anchor="w", pady=(6, 0))
        self.tv = ttk.Treeview(left, columns=("what",), show="tree headings", height=6)
        self.tv.heading("#0", text="凡例"); self.tv.heading("what", text="x ← y")
        self.tv.column("#0", width=130); self.tv.column("what", width=200)
        self.tv.pack(anchor="w")
        fb = ttk.Frame(left); fb.pack(anchor="w", pady=4)
        ttk.Button(fb, text="選んだ曲線を消す", style="Ghost.TButton", command=self.remove_series).pack(side="left")
        ttk.Button(fb, text="全部消す", style="Ghost.TButton", command=self.clear_series).pack(side="left", padx=6)
        right = ttk.Frame(main)
        right.pack(side="left", fill="both", expand=True, padx=(10, 0))
        opt = ttk.Frame(right); opt.pack(fill="x")
        self.v_title, self.v_xl, self.v_yl, self.v_smooth = tk.StringVar(), tk.StringVar(), tk.StringVar(), tk.StringVar(value="1")
        self.v_xlim, self.v_ylim = tk.StringVar(), tk.StringVar()
        for r, (lab, var, wd) in enumerate((("タイトル", self.v_title, 30), ("x 軸の名前", self.v_xl, 22), ("y 軸の名前", self.v_yl, 22))):
            ttk.Label(opt, text=lab).grid(row=0, column=2 * r, sticky="e", padx=(8, 2))
            ttk.Entry(opt, textvariable=var, width=wd).grid(row=0, column=2 * r + 1, sticky="w")
        ttk.Label(opt, text="平均する点数").grid(row=1, column=0, sticky="e", padx=(8, 2), pady=3)
        ttk.Entry(opt, textvariable=self.v_smooth, width=6).grid(row=1, column=1, sticky="w")
        ttk.Label(opt, text="x の範囲 から,まで").grid(row=1, column=2, sticky="e", padx=(8, 2))
        ttk.Entry(opt, textvariable=self.v_xlim, width=14).grid(row=1, column=3, sticky="w")
        ttk.Label(opt, text="y の範囲").grid(row=1, column=4, sticky="e", padx=(8, 2))
        ttk.Entry(opt, textvariable=self.v_ylim, width=14).grid(row=1, column=5, sticky="w")
        bar = ttk.Frame(right); bar.pack(side="bottom", fill="x", pady=(6, 0))
        ttk.Button(bar, text="描く", style="Accent.TButton", command=self.draw).pack(side="left")
        ttk.Button(bar, text="PNG に保存…", style="Ghost.TButton", command=self.save).pack(side="left", padx=8)
        self.msg = tk.StringVar(value="曲線を足して「描く」を押してください")
        ttk.Label(bar, textvariable=self.msg, foreground="#77848d").pack(side="left", padx=8)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        self.fig = Figure(figsize=(7, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, pady=(8, 0))
        if self.files:
            self.lb.selection_set(0)
            self.on_file()

    # ------------------------------------------------------------------------------------------ files and columns
    def refill(self):
        self.files = candidate_files(self.folder) + [f for f in getattr(self, "extra", [])]
        self.lb.delete(0, "end")
        for f in self.files:
            self.lb.insert("end", os.path.relpath(f, self.folder) if f.startswith(self.folder) else f)

    def add_file(self):
        fs = filedialog.askopenfilenames(parent=self.w, title="表のファイル（log.lammps、*.dat、*.csv）", initialdir=self.folder)
        if not fs:
            return
        self.extra = list(getattr(self, "extra", [])) + list(fs)
        self.refill()
        self.lb.selection_clear(0, "end"); self.lb.selection_set("end"); self.on_file()

    def blocks(self, path):
        if path not in self.cache:
            try:
                self.cache[path] = read_blocks(path)
            except Exception as e:
                messagebox.showerror("soapff-GUI", f"読めません：{e}", parent=self.w)
                self.cache[path] = []
        return self.cache[path]

    def on_file(self):
        sel = self.lb.curselection()
        if not sel:
            return
        bl = self.blocks(self.files[sel[0]])
        self.cb_block.configure(values=[b[0] for b in bl])
        if bl:
            self.v_block.set(bl[-1][0])          # the last run of log.lammps
            self.on_block()
        else:
            self.v_block.set("")
            self.cb_x.configure(values=[]); self.cb_y.configure(values=[])
            self.msg.set("このファイルから表を読めませんでした")

    def current_block(self):
        sel = self.lb.curselection()
        if not sel:
            return None
        for b in self.blocks(self.files[sel[0]]):
            if b[0] == self.v_block.get():
                return self.files[sel[0]], b
        return None

    def on_block(self):
        cb = self.current_block()
        if not cb:
            return
        cols = cb[1][1]
        self.cb_x.configure(values=cols); self.cb_y.configure(values=cols)
        self.v_x.set(cols[0]); self.v_y.set(cols[1] if len(cols) > 1 else cols[0])
        self.v_label.set(f"{os.path.basename(os.path.dirname(cb[0])) or ''} {self.v_y.get()}".strip())
        self.msg.set(f"{len(cb[1][2])} 行、{len(cols)} 列")

    # ------------------------------------------------------------------------------------------ series
    def add_series(self):
        cb = self.current_block()
        if not cb:
            self.msg.set("ファイルと表を選んでください")
            return
        try:
            xs, ys = float(self.v_xs.get()), float(self.v_ys.get())
        except ValueError:
            self.msg.set("× は数にしてください")
            return
        s = dict(label=self.v_label.get() or self.v_y.get(), path=cb[0], block=cb[1][0], x=self.v_x.get(), y=self.v_y.get(), xs=xs, ys=ys)
        self.series.append(s)
        self.tv.insert("", "end", iid=str(len(self.series) - 1), text=s["label"], values=(f"{s['x']} ← {s['y']}",))
        if not self.v_xl.get():
            self.v_xl.set(f"{s['x']}" + (f" × {xs:g}" if xs != 1 else ""))
        if not self.v_yl.get():
            self.v_yl.set(f"{s['y']}" + (f" × {ys:g}" if ys != 1 else ""))
        self.draw()

    def remove_series(self):
        sel = self.tv.selection()
        if sel:
            del self.series[int(sel[0])]
            self.reload_tree()
            self.draw()

    def clear_series(self):
        self.series = []
        self.reload_tree()
        self.draw()

    def reload_tree(self):
        self.tv.delete(*self.tv.get_children())
        for k, s in enumerate(self.series):
            self.tv.insert("", "end", iid=str(k), text=s["label"], values=(f"{s['x']} ← {s['y']}",))

    # ------------------------------------------------------------------------------------------ drawing
    @staticmethod
    def _lim(text):
        try:
            a, b = [float(v) for v in text.replace("、", ",").split(",")]
            return a, b
        except ValueError:
            return None

    def draw(self):
        self.ax.clear()
        try:
            n = max(1, int(float(self.v_smooth.get() or 1)))
        except ValueError:
            n = 1
        for s in self.series:
            blk = next((b for b in self.blocks(s["path"]) if b[0] == s["block"]), None)
            if blk is None or s["x"] not in blk[1] or s["y"] not in blk[1]:
                continue
            x = blk[2][:, blk[1].index(s["x"])] * s["xs"]
            y = blk[2][:, blk[1].index(s["y"])] * s["ys"]
            if n > 1 and len(y) > n:
                k = np.ones(n) / n
                y = np.convolve(y, k, mode="valid"); x = x[n // 2:n // 2 + len(y)]
            self.ax.plot(x, y, lw=1.4, label=s["label"])
        self.ax.set_xlabel(self.v_xl.get()); self.ax.set_ylabel(self.v_yl.get())
        if self.v_title.get():
            self.ax.set_title(self.v_title.get())
        for lim, setter in ((self._lim(self.v_xlim.get()), self.ax.set_xlim), (self._lim(self.v_ylim.get()), self.ax.set_ylim)):
            if lim:
                setter(*lim)
        if self.series:
            self.ax.legend(fontsize=9)
        self.ax.grid(alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw()
        self.msg.set(f"{len(self.series)} 本の曲線")

    def save(self):
        p = filedialog.asksaveasfilename(parent=self.w, defaultextension=".png", filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")],
                                         initialdir=self.folder, initialfile="plot.png")
        if p:
            self.fig.savefig(p, dpi=200)
            self.msg.set(f"保存しました：{p}")
