"""The fix / compute styles of LAMMPS (and a few other commands) that soapff-GUI offers to insert after the basic input
(force field + minimisation): title, description, a template line, the names of its arguments and the manual page.

Placeholders in templates: {ID} (made unique in the input), {T} temperature (K), {P} pressure (atm), {DT} time step (fs),
{N} steps.  The manual: https://docs.lammps.org/<doc>."""

DOC = "https://docs.lammps.org/"
BASE = ["ID", "対象グループ", "スタイル"]

# kind, style, title, description, template, labels after BASE (None: fixed keyword, not a field), doc page
_FIX = [
    # --- time integration and thermostats / barostats
    ("時間積分・温度・圧力", "nvt", "NVT（Nosé–Hoover の温度一定）",
     "温度一定の MD。減衰は 100 fs 程度（時間刻みの 100 倍）が目安。",
     "fix {ID} all nvt temp {T} {T} 100.0", [None, "開始温度 (K)", "終了温度 (K)", "減衰 (fs)"], "fix_nh.html"),
    ("時間積分・温度・圧力", "npt", "NPT（等方）", "温度と圧力一定。箱の 3 辺を同じ割合で変える。",
     "fix {ID} all npt temp {T} {T} 100.0 iso {P} {P} 1000.0",
     [None, "開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "方向", "開始圧力 (atm)", "終了圧力 (atm)", "減衰 (fs)"], "fix_nh.html"),
    ("時間積分・温度・圧力", "npt", "NPT（面に垂直な方向だけ：層間・界面）",
     "z 方向の箱の長さだけを変える。粘土の層間や界面に。",
     "fix {ID} all npt temp {T} {T} 100.0 z {P} {P} 1000.0",
     [None, "開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "方向", "開始圧力 (atm)", "終了圧力 (atm)", "減衰 (fs)"], "fix_nh.html"),
    ("時間積分・温度・圧力", "npt", "NPT（3 辺を別々に：バルク結晶）", "x、y、z を別々に変える。",
     "fix {ID} all npt temp {T} {T} 100.0 aniso {P} {P} 1000.0",
     [None, "開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "方向", "開始圧力 (atm)", "終了圧力 (atm)", "減衰 (fs)"], "fix_nh.html"),
    ("時間積分・温度・圧力", "npt", "NPT（箱の傾きも含めて全部）",
     "傾いた箱の全成分を変える（直方体の箱では先に change_box all triclinic が要る）。",
     "fix {ID} all npt temp {T} {T} 100.0 tri {P} {P} 1000.0",
     [None, "開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "方向", "開始圧力 (atm)", "終了圧力 (atm)", "減衰 (fs)"], "fix_nh.html"),
    ("時間積分・温度・圧力", "nphug", "NPHug（Hugoniostat：衝撃圧縮のユゴニオ状態）",
     "衝撃波のユゴニオ曲線上の状態（圧力・体積・エネルギーがユゴニオ条件を満たす）に系を導く。圧力は衝撃の圧力（atm）、"
     "温度の開始・終了値は使われず減衰だけが効く。基準状態（E0、V0、P0）は最初のステップの値（キーワードで指定も可）。"
     "z は一軸の衝撃（z 方向だけ箱が変わる）。",
     "fix {ID} all nphug temp {T} {T} 100.0 z 10000.0 10000.0 1000.0",
     [None, "開始温度（未使用）", "終了温度（未使用）", "減衰 (fs)", "方向", "開始圧力 (atm)", "終了圧力 (atm)", "減衰 (fs)"], "fix_nphug.html"),
    ("時間積分・温度・圧力", "nph", "NPH（圧力一定、温度制御なし）", "",
     "fix {ID} all nph iso {P} {P} 1000.0", ["方向", "開始圧力 (atm)", "終了圧力 (atm)", "減衰 (fs)"], "fix_nh.html"),
    ("時間積分・温度・圧力", "nve", "NVE（エネルギー一定）", "時間積分だけ。温度制御の fix（langevin など）と組み合わせることもある。",
     "fix {ID} all nve", [], "fix_nve.html"),
    ("時間積分・温度・圧力", "langevin", "Langevin 熱浴", "fix nve と一緒に使う。乱数の種は正の整数。",
     "fix {ID} all langevin {T} {T} 100.0 48279", ["開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "乱数の種"], "fix_langevin.html"),
    ("時間積分・温度・圧力", "temp/csvr", "CSVR 熱浴（Bussi）", "fix nve と一緒に使う。",
     "fix {ID} all temp/csvr {T} {T} 100.0 54324", ["開始温度 (K)", "終了温度 (K)", "減衰 (fs)", "乱数の種"], "fix_temp_csvr.html"),
    ("時間積分・温度・圧力", "temp/berendsen", "Berendsen 熱浴", "fix nve と一緒に使う。平衡化向き。",
     "fix {ID} all temp/berendsen {T} {T} 100.0", ["開始温度 (K)", "終了温度 (K)", "減衰 (fs)"], "fix_temp_berendsen.html"),
    ("時間積分・温度・圧力", "temp/rescale", "速度のスケーリング", "N ステップごとに温度を目標に戻す。fix nve と一緒に。",
     "fix {ID} all temp/rescale 100 {T} {T} 10.0 1.0",
     ["間隔 (ステップ)", "開始温度 (K)", "終了温度 (K)", "許容幅 (K)", "割合"], "fix_temp_rescale.html"),
    ("時間積分・温度・圧力", "press/berendsen", "Berendsen 圧力浴", "温度制御の fix と一緒に使う。",
     "fix {ID} all press/berendsen iso {P} {P} 1000.0", ["方向", "開始圧力 (atm)", "終了圧力 (atm)", "減衰 (fs)"],
     "fix_press_berendsen.html"),
    ("時間積分・温度・圧力", "box/relax", "箱の最小化（minimize の前に置く）", "エネルギー最小化で箱の大きさも緩める。",
     "fix {ID} all box/relax iso 0.0 vmax 0.001", ["方向", "目標圧力 (atm)", None, "1 回の最大変化"], "fix_box_relax.html"),
    ("時間積分・温度・圧力", "rigid", "剛体（分子ごと）", "分子を剛体として動かす（温度制御は rigid/nvt など）。",
     "fix {ID} all rigid molecule", ["単位"], "fix_rigid.html"),
    # --- constraints and restraints
    ("拘束・束縛", "shake", "SHAKE（結合と角度の固定）",
     "b と a の後は結合・角度の型番号（system.settings の bond_coeff / angle_coeff）。柔軟な SPC 水を剛体にするとき。",
     "fix {ID} water shake 1.0e-4 20 0 b 2 a 1", ["許容誤差", "最大反復", "出力間隔", None, "結合の型", None, "角度の型"],
     "fix_shake.html"),
    ("拘束・束縛", "rattle", "RATTLE（SHAKE と速度の拘束）", "", "fix {ID} water rattle 1.0e-4 20 0 b 2 a 1",
     ["許容誤差", "最大反復", "出力間隔", None, "結合の型", None, "角度の型"], "fix_shake.html"),
    ("拘束・束縛", "setforce", "力を固定（0 なら原子を止める）", "NULL の成分はそのまま。",
     "fix {ID} all setforce 0.0 0.0 0.0", ["Fx", "Fy", "Fz"], "fix_setforce.html"),
    ("拘束・束縛", "spring/self", "初めの位置へのばね", "原子をその場所につなぐ。",
     "fix {ID} all spring/self 10.0", ["ばね定数 (kcal/mol/Å²)"], "fix_spring_self.html"),
    ("拘束・束縛", "momentum", "全体の運動量を 0 に", "系の並進（と回転）を止める。",
     "fix {ID} all momentum 100 linear 1 1 1", ["間隔 (ステップ)", None, "x", "y", "z"], "fix_momentum.html"),
    ("拘束・束縛", "recenter", "重心を固定", "", "fix {ID} all recenter INIT INIT INIT", ["x", "y", "z"], "fix_recenter.html"),
    ("拘束・束縛", "restrain", "結合・角度などの束縛", "原子番号の組と力の定数。",
     "fix {ID} all restrain bond 1 2 1000.0 1000.0 1.0", [None, "原子 1", "原子 2", "開始の強さ", "終了の強さ", "距離 (Å)"],
     "fix_restrain.html"),
    ("拘束・束縛", "viscous", "粘性抵抗", "", "fix {ID} all viscous 0.1", ["係数"], "fix_viscous.html"),
    # --- external fields and walls
    ("外場・壁", "efield", "電場", "real 単位では V/Å。", "fix {ID} all efield 0.0 0.0 0.1", ["Ex", "Ey", "Ez"], "fix_efield.html"),
    ("外場・壁", "addforce", "力を加える", "", "fix {ID} all addforce 0.0 0.0 -0.1", ["Fx", "Fy", "Fz"], "fix_addforce.html"),
    ("外場・壁", "aveforce", "力をグループで平均", "", "fix {ID} all aveforce NULL NULL 0.0", ["Fx", "Fy", "Fz"], "fix_aveforce.html"),
    ("外場・壁", "wall/reflect", "反射する壁", "周期でない面で原子が失われないように。", "fix {ID} all wall/reflect zlo EDGE zhi EDGE",
     ["面", "位置", "面", "位置"], "fix_wall_reflect.html"),
    ("外場・壁", "wall/lj93", "LJ 9-3 の壁", "", "fix {ID} all wall/lj93 zlo EDGE 1.0 1.0 2.5",
     ["面", "位置", "ε", "σ", "カットオフ"], "fix_wall.html"),
    ("外場・壁", "deform", "箱を変形させる", "一定のひずみ速度で引っ張るなど。", "fix {ID} all deform 1 z erate 1.0e-6 remap x",
     ["間隔", "方向", "方法", "値", None, "remap"], "fix_deform.html"),
    ("外場・壁", "heat", "熱を加える・取る", "", "fix {ID} all heat 10 0.1", ["間隔", "エネルギー (kcal/mol/fs)"], "fix_heat.html"),
    # --- output / averaging
    ("出力・平均", "ave/time", "時間平均してファイルに出す", "compute の値（c_ID）を Nfreq ステップごとに書く。",
     "fix {ID} all ave/time 10 100 1000 c_thermo_temp file ave_time.dat", ["Nevery", "Nrepeat", "Nfreq", "値", None, "ファイル"],
     "fix_ave_time.html"),
    ("出力・平均", "ave/chunk", "chunk ごとの平均（密度分布など）",
     "先に compute chunk/atom（ID cz）を置く。z 方向の質量密度の分布など。",
     "fix {ID} all ave/chunk 10 100 1000 cz density/mass file density_z.dat",
     ["Nevery", "Nrepeat", "Nfreq", "chunk の ID", "値", None, "ファイル"], "fix_ave_chunk.html"),
    ("出力・平均", "ave/correlate", "時間相関（自己相関など）", "", "fix {ID} all ave/correlate 10 100 1000 c_thermo_press[4] type auto file corr.dat",
     ["Nevery", "Nrepeat", "Nfreq", "値", None, "相関の種類", None, "ファイル"], "fix_ave_correlate.html"),
    ("出力・平均", "ave/atom", "原子ごとの時間平均", "", "fix {ID} all ave/atom 10 100 1000 x y z",
     ["Nevery", "Nrepeat", "Nfreq", "値", "値", "値"], "fix_ave_atom.html"),
    ("出力・平均", "dt/reset", "時間刻みを自動で変える", "", "fix {ID} all dt/reset 10 NULL 1.0 0.1 units box",
     ["間隔", "最小 (fs)", "最大 (fs)", "最大移動量 (Å)"], "fix_dt_reset.html"),
]

_COMPUTE = [
    ("拡散・構造", "msd", "平均二乗変位（MSD）", "fix ave/time で書き出す（c_ID[4] が全体）。",
     "compute {ID} all msd", [], "compute_msd.html"),
    ("拡散・構造", "msd/chunk", "chunk ごとの MSD", "先に compute chunk/atom を置く。", "compute {ID} all msd/chunk cmol",
     ["chunk の ID"], "compute_msd_chunk.html"),
    ("拡散・構造", "vacf", "速度自己相関", "", "compute {ID} all vacf", [], "compute_vacf.html"),
    ("拡散・構造", "rdf", "動径分布関数（RDF）", "fix ave/time ... mode vector で書き出す。",
     "compute {ID} all rdf 100", ["ビンの数"], "compute_rdf.html"),
    ("拡散・構造", "displace/atom", "原子ごとの変位", "", "compute {ID} all displace/atom", [], "compute_displace_atom.html"),
    ("拡散・構造", "coord/atom", "配位数", "", "compute {ID} all coord/atom cutoff 3.0", [None, "カットオフ (Å)"],
     "compute_coord_atom.html"),
    ("拡散・構造", "cluster/atom", "クラスター（つながり）", "", "compute {ID} all cluster/atom 3.5", ["カットオフ (Å)"],
     "compute_cluster_atom.html"),
    ("拡散・構造", "centro/atom", "中心対称性パラメータ", "", "compute {ID} all centro/atom 12", ["近接の数"],
     "compute_centro_atom.html"),
    ("拡散・構造", "orientorder/atom", "配向秩序（Steinhardt）", "", "compute {ID} all orientorder/atom", [],
     "compute_orientorder_atom.html"),
    ("拡散・構造", "gyration", "回転半径", "", "compute {ID} all gyration", [], "compute_gyration.html"),
    ("拡散・構造", "com", "重心", "", "compute {ID} all com", [], "compute_com.html"),
    ("拡散・構造", "chunk/atom", "原子を chunk に分ける（z のビンなど）", "fix ave/chunk や msd/chunk と組み合わせる。",
     "compute {ID} all chunk/atom bin/1d z lower 0.5 units box", ["方法", "方向", "原点", "幅 (Å)", None, "単位"],
     "compute_chunk_atom.html"),
    ("温度・エネルギー・圧力", "temp", "温度", "", "compute {ID} all temp", [], "compute_temp.html"),
    ("温度・エネルギー・圧力", "temp/com", "重心の運動を除いた温度", "", "compute {ID} all temp/com", [], "compute_temp_com.html"),
    ("温度・エネルギー・圧力", "temp/partial", "一部の方向だけの温度", "", "compute {ID} all temp/partial 1 1 0", ["x", "y", "z"],
     "compute_temp_partial.html"),
    ("温度・エネルギー・圧力", "pe", "ポテンシャルエネルギー", "", "compute {ID} all pe", [], "compute_pe.html"),
    ("温度・エネルギー・圧力", "pe/atom", "原子ごとのポテンシャルエネルギー", "", "compute {ID} all pe/atom", [], "compute_pe_atom.html"),
    ("温度・エネルギー・圧力", "ke", "運動エネルギー", "", "compute {ID} all ke", [], "compute_ke.html"),
    ("温度・エネルギー・圧力", "ke/atom", "原子ごとの運動エネルギー", "", "compute {ID} all ke/atom", [], "compute_ke_atom.html"),
    ("温度・エネルギー・圧力", "pressure", "圧力テンソル", "", "compute {ID} all pressure thermo_temp", ["温度の compute"],
     "compute_pressure.html"),
    ("温度・エネルギー・圧力", "stress/atom", "原子ごとの応力", "", "compute {ID} all stress/atom NULL", ["温度の compute"],
     "compute_stress_atom.html"),
    ("温度・エネルギー・圧力", "heat/flux", "熱流束", "先に ke/atom、pe/atom、stress/atom を置く。",
     "compute {ID} all heat/flux ke pe stress", ["ke/atom の ID", "pe/atom の ID", "stress/atom の ID"], "compute_heat_flux.html"),
    ("温度・エネルギー・圧力", "group/group", "グループ間の相互作用", "", "compute {ID} water group/group all", ["相手のグループ"],
     "compute_group_group.html"),
    ("温度・エネルギー・圧力", "pair", "pair の成分ごとのエネルギー", "", "compute {ID} all pair lj/cut/coul/long", ["pair のスタイル"],
     "compute_pair.html"),
    ("原子・結合の量", "property/atom", "原子の量（電荷、速度など）", "", "compute {ID} all property/atom q", ["量"],
     "compute_property_atom.html"),
    ("原子・結合の量", "reduce", "原子ごとの量の和・平均", "", "compute {ID} all reduce sum c_pea", ["演算", "値"], "compute_reduce.html"),
    ("原子・結合の量", "count/type", "型ごとの原子数", "", "compute {ID} all count/type atom", ["対象"], "compute_count_type.html"),
    ("原子・結合の量", "bond", "結合のエネルギー", "", "compute {ID} all bond", [], "compute_bond.html"),
    ("原子・結合の量", "angle", "角度のエネルギー", "", "compute {ID} all angle", [], "compute_angle.html"),
    ("原子・結合の量", "bond/local", "結合ごとの長さ", "", "compute {ID} all bond/local dist", ["量"], "compute_bond_local.html"),
    ("原子・結合の量", "angle/local", "角度ごとの値", "", "compute {ID} all angle/local theta", ["量"], "compute_angle_local.html"),
    ("原子・結合の量", "dipole", "双極子モーメント", "", "compute {ID} all dipole", [], "compute_dipole.html"),
]

_OTHER = [
    ("MD の基本", "run", "MD を回す", "{N} ステップ。時間 = ステップ × 時間刻み。", "run {N}", ["ステップ数"], "run.html"),
    ("MD の基本", "timestep", "時間刻み", "柔軟な SPC 水（ClayFF）では 0.5 fs。", "timestep {DT}", ["時間刻み (fs)"], "timestep.html"),
    ("MD の基本", "velocity", "初速度を与える", "温度 {T} K のマクスウェル分布。", "velocity all create {T} 4928459 mom yes rot no dist gaussian",
     ["対象グループ", "方法", "温度 (K)", "乱数の種"], "velocity.html"),
    ("MD の基本", "unfix", "fix を外す", "", "unfix {LASTFIX}", ["ID"], "unfix.html"),
    ("MD の基本", "uncompute", "compute を外す", "", "uncompute {LASTCOMPUTE}", ["ID"], "uncompute.html"),
    ("MD の基本", "reset_timestep", "ステップを 0 に戻す", "", "reset_timestep 0", ["ステップ"], "reset_timestep.html"),
    ("MD の基本", "minimize", "エネルギー最小化", "", "minimize 1.0e-6 1.0e-6 5000 50000",
     ["エネルギーの許容値", "力の許容値", "最大反復", "最大評価回数"], "minimize.html"),
    ("出力", "thermo", "画面出力の間隔", "", "thermo 1000", ["間隔 (ステップ)"], "thermo.html"),
    ("出力", "thermo_style", "画面出力の項目", "", "thermo_style custom step temp press pe ke etotal vol density lx ly lz",
     ["スタイル"], "thermo_style.html"),
    ("出力", "dump", "トラジェクトリ（型の名前と電荷つき、OVITO 用）", "element 列に ClayFF の型の名前を書く。",
     "dump {ID} all custom 2000 traj.lammpstrj id type element x y z ix iy iz q", ["ID", "対象グループ", "スタイル", "間隔 (ステップ)", "ファイル"],
     "dump.html"),
    ("出力", "write_data", "最後の構造を書き出す", "", "write_data final.data", ["ファイル"], "write_data.html"),
    ("出力", "write_restart", "リスタートファイル", "", "write_restart restart.bin", ["ファイル"], "write_restart.html"),
    ("グループ", "group", "型でグループを作る", "型の番号は system.data の Masses。", "group {ID} type 1", ["グループ名", "方法", "型"],
     "group.html"),
]

# special entries built by the application (need the type numbers of the system)
SPECIAL = [
    ("soapff の手順", "soapff:diffusion", "拡散係数の MSD（soapff の解析用）",
     "水の O と各イオンのグループを作り、compute msd を fix ave/time で msd_<名前>.dat に 0.1 ps ごとに書きます。"
     "「解析」はこのファイルから拡散係数を求めます。run の前に置いてください。"),
    ("soapff の手順", "soapff:protocol", "soapff の緩和（初速度 → 緩和 NVT）",
     "初速度、時間刻み、緩和 NVT（記録しない）、ステップのリセットまでをまとめて入れます。本計算（fix nvt / npt、dump、"
     "拡散係数の MSD、run、write_data）は、そのあとに一覧から選んで足します。長さと温度は計算設定の MD の条件です。"),
    ("soapff の手順", "soapff:swelling", "膨潤曲線（層間に水を加えながら NPT、d001 を記録）",
     "Cygan ら（2004）の方法：層間に水分子を 1 個ずつ（その後 2 個ずつ）加え、各状態で最小化と NPT（平衡化のあと平均）を"
     "回して d001（= lz）を swelling.dat に書きます。最後に 2 個ずつ減らす脱水も選べます。構造に水の型がなければ、"
     "system.data に水の型を足します（water.mol も書きます）。粘土は箱の中央に保ち、層間は周期境界をまたぐ側です。"
     "swelling.dat の列は step、加えた水の数、段階（1 = 記録）、lx ly lz xy xz yz、PE、T です（例の analyze_swelling.py で"
     "d001 と水の量の図にできます）。in.lammps の最小化のすぐ後に置いてください。"),
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
