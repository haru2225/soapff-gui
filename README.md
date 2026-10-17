# soapff-GUI：粘土鉱物の初期構造を選ぶだけで LAMMPS を回す

[soapff-cui](../soapff-cui) のウィンドウ版です。初期構造のファイルを選んで「実行」を押すと、

1. soapff-cui が出力フォルダに LAMMPS の入力一式（`system.data`、`system.settings`、`in.lammps`、`charges.csv`、`report.json`）を書き出し、
2. そのまま LAMMPS（`mpirun -np N lmp -in in.lammps`）を起動します。

入力ファイルを開いたり書き換えたりする必要はありません。SOAP のライブラリは、隣の soapff-cui の `examples/minerals/library_clayff.npz`（いちばん新しいもの）を自動で使うので、指定は不要です。それが見つからないときは、GUI に同梱したコピー（`soapff_gui/data/library_clayff.npz`、soapff-cui と同じ 5 つの系、5,604 原子）を使います。

## 計算できる鉱物と界面構造

いまのライブラリ（`library_clayff.npz`、5 つの系、5,604 原子）で確かめたのは次の鉱物と構造です。どれも全原子が SOAP で決まり（fallback 0）、電荷と型が参照（または元の入力、結合グラフから決めた ClayFF の型）と一致し、生成した入力で MD が安定に回りました。

| 鉱物 | 構造 | 確かめたこと |
|---|---|---|
| Na モンモリロナイト | Mg 置換の 2:1 層 2 枚、層間に Na⁺ 30 個と水 2 層（2 水和層）、2,502 原子（JAEA の入力） | 電荷・型・結合・LJ が元の入力と一致。水が粘土に重なった構造でも可 |
| 石英 | (101) と水の界面（完全に水酸基化、シラノール 5.9 個/nm²、水の層 30 Å）、819 原子 / 水のない結晶 | 界面：参照と一致。結晶：型が一致 |
| アルバイト | (001) と水の界面（Si–OH と Al(IV)–OH で終端、Na⁺ を含む）、820 原子 / 水のない結晶 | 同上 |
| ブルーサイト | (001) と水の界面（層の間で切った Mg(OH)₂ 3 層）、624 原子 / 水のない結晶 | 同上 |
| パイロフィライト | (001) と水の界面（層の間で切った 2:1 層 2 枚）、839 原子 / 水のない結晶 | 同上 |

界面構造は soapff-cui の `examples/minerals/build_ase.py` で、次の手順で作ったものです。
1. 結晶：CIF を ASE で読む（石英・アルバイト・パイロフィライトは COD 1011097、9000584、9000809。パイロフィライトは OH の H を付ける）。ブルーサイトは文献の構造（Catti et al. 1995）から ASE で作る
2. スラブ：`ase.build.surface` で面を xy にし、厚さ 12 Å 以上の整数周期を積む。切れる陽イオン–O 結合がいちばん少ない高さで切る（O–H は切らない。層状鉱物は層の間）
3. 終端：切れた結合を、解離した水として OH（陽イオン側）と H（O 側）で終端する（完全に水酸基化した表面。電荷は自動で中性）。面内は各辺 10 Å 以上に繰り返す
4. 水：スラブの上に 30 Å の空間を設け、SPC/E の水を密度 0.997 g/cm³（30 − 2 Å の厚さに対して）でランダムに置く（ほかの分子の原子と 2.0 Å 以上。Packmol と同じ許容距離）
5. 緩和：LAMMPS（ClayFF + SPC/E）で、鉱物の H 以外を固定し、水と H だけを最小化 → 20 ps NVT（300 K）→ 最小化

真空に面したスラブ（境界条件 `p p f` など）：
- パイロフィライト (001)（2:1 層 2 枚、320 原子）は全原子が SOAP で決まり、LAMMPS も通りました。
- ブルーサイト (001) は、真空に面した表面がライブラリにありません（ライブラリは水との界面だけ）。表面の OH の H が類似度 0.964 で水の H と判定され、LAMMPS が止まります。ライブラリにない環境では、類似度 0.95 以上でも型を取り違えることがあります。

確かめていないもの：
- 同じ鉱物・同じ面で、面内の大きさ、スラブの厚さ、水の層の厚さを変えた構造（原子のまわり 5 Å の環境は同じなので、扱える見込み）
- ほかの面、エッジ、ほかの層間陽イオン（K⁺、Cs⁺、Ca²⁺）、水の量の違うモンモリロナイト、ほかの鉱物
- ClayFF に型のない元素（Ti、S など）を含む鉱物（力場を作れません）

こうした構造を渡しても力場は作られますが、類似度が 0.95 未満の原子は ClayFF の型とその電荷に戻され、0.95 以上でもライブラリにない環境の可能性があります。`charges.csv`（soapff-gui では「状態」の最小類似度と警告）で確かめ、必要ならその系の参照をライブラリに加えてください。

## インストール（GitHub から）

```
# 1. soapff-cui（力場を作る本体）
git clone git@github.com:haru2225/soapff-cui.git && pip install -e soapff-cui
# 2. soapff-gui
git clone git@github.com:haru2225/soapff-gui.git && pip install -e soapff-gui
soapff-gui              # = python -m soapff_gui
```

Python には tkinter が必要です（多くの配布に含まれます）。LAMMPS（`lmp`）と `mpirun` は別に必要です。

## 起動

```
# A. macOS：Finder で soapff-gui.command をダブルクリック
#    （Python は既定で ~/ScoreMD/.pixi/envs/default/bin/python。PY で変えられます）

# B. ターミナルから
PYTHONPATH=~/clay/soapff-gui ~/ScoreMD/.pixi/envs/default/bin/python -m soapff_gui
```

soapff-cui がインストールされていないときは、`$SOAPFF_HOME`、このフォルダの隣の `soapff-cui/`（`~/clay/soapff-cui`）、`soapff/`（古い名前）の順に探します。見つからないときは、起動したときに知らせます。
LAMMPS（`lmp`、`lmp_mpi`、`lmp_serial`）と `mpirun` は、Python と同じ環境の `bin` を先に、次に PATH を探します。見つからないときは「詳細設定」で場所を指定してください。

## 画面

| 項目 | 内容 |
|---|---|
| 初期構造 | extxyz（Lattice 付き）、cif、POSCAR、LAMMPS データ。元素・座標・箱だけで足ります。選ぶと原子数・元素・箱の大きさを表示します |
| 出力フォルダ | 既定は構造ファイルの隣の `<名前>_soapff/`。前の結果があるときは上書きを確認します |
| MD の条件 | 温度、NVT の長さ、NPT の長さ（0 でなし）、圧力、圧力の方向（界面・層間は `z`、バルク結晶は `aniso` / `tri`）、並列数（MPI） |
| 境界条件 | x・y・z それぞれ `p`（周期）、`f`（固定）、`m`、`s`、面ごとに違う `fm` など（LAMMPS の boundary と同じ）。既定は `p p p`。周期でない方向があるときは NPT は使えません。長距離クーロン力の計算法は soapff-cui が選びます（soapff-cui の README の「境界条件」） |
| 詳細設定 | SOAP ライブラリ、lmp と mpirun の場所、類似度のしきい値（既定 0.95）、カットオフ（既定 10 Å） |
| ▶ 実行 / ■ 停止 | 力場の作成 → 最小化 → MD。停止は LAMMPS を止めます |
| 入力だけ作る | LAMMPS を回さず、入力一式だけ書き出します（スパコンに持っていくときなど） |
| OVITO で開く | 出力フォルダの `traj.lammpstrj`（なければ `final.data`、`system.data`）を OVITO で開きます。計算中でも、それまでのフレームを開けます。RDF、配位数、密度分布などの解析は OVITO の中で行います。macOS は Ovito.app、それ以外は `ovito` コマンドが必要です |
| 状態 | 段階（力場作成・最小化・MD）、MD の進み具合、ステップ・温度・圧力・ポテンシャルエネルギー・密度 |

力場を作ったあと、原子数、SOAP で電荷を決めた原子と fallback の原子の数、最小の類似度を表示します。fallback の原子や重なりの警告があるときは橙色で知らせます（ライブラリにない環境が含まれている可能性があります。`charges.csv` を確認してください）。

設定は `~/.config/soapff-gui/settings.json` に保存し、次回の起動で使います。

## 出力（出力フォルダ）

soapff-cui の出力（`system.data`、`system.settings`、`in.lammps`、`charges.csv`、`report.json`）に加えて、LAMMPS の `log.lammps`、`traj.lammpstrj`（1 ps ごと）、`final.data`（最後の構造）。

## 構成

- `soapff_gui/runner.py`：GUI を使わない部分（soapff と LAMMPS を探す、soapff-cui の `make` を呼ぶ、LAMMPS を起動する、画面出力から段階と thermo の値と進み具合を読む）。
- `soapff_gui/app.py`：tkinter のウィンドウ。計算は別スレッドで行い、画面はキューで更新します。

## 限界

- 力場とライブラリの限界は soapff-cui と同じです（上の「計算できる鉱物と界面構造」）。
- 1 台の計算機で回すためのものです。PBS などのジョブ投入はしません（「入力だけ作る」で入力を作って持っていってください）。
