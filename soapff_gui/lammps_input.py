"""in.lammps as a list of lines for the editor of soapff-GUI: split a line into command / arguments / comment, name the
arguments of the common commands (the form under the editor), count the MD steps of the runs (progress bar)."""
import re

# command -> labels of its first arguments; the remaining arguments get generic labels
SCHEMA = {
    "units": ["単位系"],
    "atom_style": ["原子のスタイル"],
    "boundary": ["x", "y", "z"],
    "read_data": ["データファイル"],
    "include": ["ファイル"],
    "timestep": ["時間刻み (fs)"],
    "thermo": ["出力間隔 (ステップ)"],
    "thermo_style": ["スタイル"],
    "neighbor": ["スキン (Å)", "方式"],
    "neigh_modify": ["キーワード"],
    "minimize": ["エネルギーの許容値", "力の許容値", "最大反復", "最大評価回数"],
    "min_style": ["方式"],
    "min_modify": ["キーワード", "値"],
    "run": ["ステップ数"],
    "reset_timestep": ["ステップ"],
    "velocity": ["対象グループ", "方法", "温度 (K)", "乱数の種"],
    "group": ["グループ名", "方法"],
    "unfix": ["ID"],
    "undump": ["ID"],
    "dump": ["ID", "対象グループ", "スタイル", "間隔 (ステップ)", "ファイル"],
    "dump_modify": ["ID"],
    "compute": ["ID", "対象グループ", "スタイル"],
    "write_data": ["ファイル"],
    "write_restart": ["ファイル"],
    "change_box": ["対象グループ"],
}
FIX_STYLE = {
    # None: a fixed keyword (kept, not shown as a field)
    "nvt": ["ID", "対象グループ", "スタイル", None, "開始温度 (K)", "終了温度 (K)", "減衰 (fs)"],
    "npt": ["ID", "対象グループ", "スタイル", None, "開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "方向", "開始圧力 (atm)",
            "終了圧力 (atm)", "減衰 (fs)"],
    "nve": ["ID", "対象グループ", "スタイル"],
    "langevin": ["ID", "対象グループ", "スタイル", "開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "乱数の種"],
    "ave/time": ["ID", "対象グループ", "スタイル", "Nevery", "Nrepeat", "Nfreq", "値"],
    "shake": ["ID", "対象グループ", "スタイル", "許容誤差", "反復", "出力間隔"],
    "wall/reflect": ["ID", "対象グループ", "スタイル"],
}

# commands offered by "＋ コマンドを挿入" (label, line)
TEMPLATES = [
    ("NVT（温度一定）", "fix nvt all nvt temp 300 300 100"),
    ("NPT（温度・圧力一定、面に垂直）", "fix npt all npt temp 300 300 100 z 1.0 1.0 1000"),
    ("fix を外す", "unfix nvt"),
    ("MD を回す", "run 20000"),
    ("時間刻み", "timestep 0.5"),
    ("出力間隔（thermo）", "thermo 1000"),
    ("速度を与える", "velocity all create 300 4928459 mom yes rot no dist gaussian"),
    ("エネルギー最小化", "minimize 1.0e-6 1.0e-6 5000 50000"),
    ("トラジェクトリの出力", "dump d2 all custom 2000 traj2.lammpstrj id type x y z"),
    ("平均二乗変位", "compute msd_all all msd"),
    ("構造の書き出し", "write_data final.data"),
    ("ステップを 0 に戻す", "reset_timestep 0"),
    ("コメント", "# コメント"),
]


def split(line):
    """(command, args, comment) of one input line; command '' for a blank or comment-only line."""
    body, sep, comment = line.partition("#")
    toks = body.split()
    return (toks[0] if toks else ""), toks[1:], (sep + comment).rstrip() if sep else ""


def join(command, args, comment=""):
    s = " ".join([command] + [a for a in args if a != ""])
    return (s + ("  " + comment if comment and s else comment)).rstrip()


def labels(command, args):
    """labels of the fields of a command (one per argument; at least the schema's); None for a fixed keyword."""
    from .catalog import labels_for
    base = labels_for(command, args)
    if base is not None:
        n = max(len(args), len(base))
        return [base[i] if i < len(base) else f"引数 {i + 1}" for i in range(n)]
    if command == "fix":
        style = args[2] if len(args) > 2 else ""
        base = FIX_STYLE.get(style, ["ID", "対象グループ", "スタイル"])
    else:
        base = SCHEMA.get(command, [])
    n = max(len(args), len(base))
    return [base[i] if i < len(base) else f"引数 {i + 1}" for i in range(n)]


def run_steps(lines):
    """the step counts of the `run` commands in order (variables or invalid values count as 0)."""
    out = []
    for line in lines:
        c, a, _ = split(line)
        if c == "run" and a:
            try:
                out.append(int(float(a[0])))
            except ValueError:
                out.append(0)
    return out


def timestep(lines, default=1.0):
    """the (last) timestep in fs."""
    t = default
    for line in lines:
        c, a, _ = split(line)
        if c == "timestep" and a:
            try:
                t = float(a[0])
            except ValueError:
                pass
    return t


def describe_run(lines, index):
    """a short Japanese description of what the run at line `index` is (from the fix that is active there)."""
    fixes = {}
    for line in lines[:index]:
        c, a, _ = split(line)
        if c == "fix" and len(a) >= 3 and a[2] in ("nvt", "npt", "nve", "langevin"):
            fixes[a[0]] = a[2]
        elif c == "unfix" and a:
            fixes.pop(a[0], None)
    kinds = sorted(set(fixes.values()))
    return "・".join(k.upper() for k in kinds) if kinds else "MD"


KEYWORD = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_/]*)")


class StepTracker:
    """progress over all `run` commands of an input from the LAMMPS screen output.  Each 'Setting up Verlet run'
    starts the next run; its first thermo step is the origin (the step counter may continue or have been reset)."""

    def __init__(self, lines):
        self.lines = list(lines)
        self.runs = run_steps(self.lines)
        self.run_lines = [i for i, l in enumerate(self.lines) if split(l)[0] == "run"]
        self.total = max(sum(self.runs), 1)
        self.k = -1                    # index of the current run (-1: before the first)
        self.start = None
        self.phase = "準備中"
        self.columns = None
        self.values = {}

    def feed(self, line):
        s = line.strip()
        if s.startswith("Setting up cg style minimization") or s.startswith("Setting up") and "minimization" in s:
            self.phase, self.columns, self.values = "最小化", None, {}
            return True
        if s.startswith("Setting up Verlet run"):
            self.k += 1
            self.start = None
            what = describe_run(self.lines, self.run_lines[self.k]) if 0 <= self.k < len(self.run_lines) else "MD"
            self.phase = f"{what}（{self.k + 1}/{len(self.runs)}）" if self.runs else "MD"
            self.columns, self.values = None, {}
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
                if self.k >= 0 and self.start is None:
                    self.start = self.values.get("Step", 0.0)
                return True
        return False

    @property
    def progress(self):
        if self.k < 0 or "Step" not in self.values or self.start is None:
            return 0.0
        done = sum(self.runs[:self.k]) + max(0.0, self.values["Step"] - self.start)
        return min(done / self.total, 1.0)


TAIL = "# ここから下に MD（fix、compute、run など）を追加します（＋ fix / ＋ compute / ＋ その他）"


def after_minimize(lines):
    """the command lines after the last minimize that basic_part would remove (blank lines, comments and walls kept)."""
    last = max((i for i, l in enumerate(lines) if split(l)[0] == "minimize"), default=None)
    if last is None:
        return []
    return [l for l in lines[last + 1:] if split(l)[0] and not (split(l)[0] == "fix" and len(split(l)[1]) > 2
                                                                 and split(l)[1][2].startswith("wall"))]


def basic_part(lines):
    """the basic input: everything up to and including the last minimize (force field, neighbour list, thermo output,
    minimisation); the MD (fix, compute, run ...) is added by the user from the catalog."""
    last = max((i for i, l in enumerate(lines) if split(l)[0] == "minimize"), default=None)
    if last is None:
        return list(lines)
    # walls of non-periodic faces belong to the boundary (system.settings), keep them
    walls = [l for l in lines[last + 1:] if split(l)[0] == "fix" and len(split(l)[1]) > 2 and split(l)[1][2].startswith("wall")]
    keep = list(lines[:last + 1])
    if walls:
        keep += ["# walls of the non-periodic faces (boundary)"] + walls
    return keep + ["", TAIL]


def ids(lines, command):
    """IDs of the fixes / computes / dumps defined in the lines (minus the ones removed by unfix / uncompute / undump)."""
    out = []
    un = {"fix": "unfix", "compute": "uncompute", "dump": "undump"}[command]
    for line in lines:
        c, a, _ = split(line)
        if c == command and a:
            out.append(a[0])
        elif c == un and a and a[0] in out:
            out.remove(a[0])
    return out


def unique_id(lines, command, want):
    """want, or want2, want3 ... not yet used by a fix / compute / dump of that kind."""
    used = set()
    for line in lines:
        c, a, _ = split(line)
        if c == command and a:
            used.add(a[0])
    if want not in used:
        return want
    k = 2
    while f"{want}{k}" in used:
        k += 1
    return f"{want}{k}"
